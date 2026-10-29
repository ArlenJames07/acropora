#!/usr/bin/env python3
"""Step 10 - sequence conservation and species tree from single-copy orthologs.

Uses families with exactly one intact locus in each genome (grades G1/G2), so
paralogy and annotation differences cannot affect the result.

  * CDS of each locus extracted with gffread from the uniform loci (step 05)
  * protein alignment with MAFFT, codon alignment by back-translation
  * yn00 (PAML) pairwise dN, dS and dN/dS for all six species pairs
  * per-gene quartet topology (four-point condition on protein p-distances)
  * concatenated protein alignment -> IQ-TREE species tree with UFBoot and
    gene concordance factors

Outputs (10_sequence_evolution/): ortholog_evolution.tsv, pair_summary.tsv,
quartet_support.tsv, species_tree.* (IQ-TREE)
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from statistics import median

from config import CODES, ENV_BIN, INPUTS, OUT, ensure, read_fasta, write_fasta

LOCI = OUT / "05_loci"
GC = OUT / "09_gene_conservation"
EVO = ensure(OUT / "10_sequence_evolution")
THREADS = int(os.environ.get("THREADS", "16"))
PAIRS = [(a, b) for i, a in enumerate(CODES) for b in CODES[i + 1:]]
CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def translate(cds):
    return "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def extract_cds():
    cds = {}
    for code in CODES:
        out = EVO / f"{code}.loci.cds.fna"
        if not out.exists():
            subprocess.run([str(ENV_BIN / "gffread"), "-x", str(out), "-g",
                            str(INPUTS / "genomes" / f"{code}.fa"), str(LOCI / f"{code}.loci.gff3")], check=True)
        cds.update({name: seq.upper() for name, seq in read_fasta(out)})
    return cds


def single_copy_families():
    loci_of = defaultdict(dict)
    status = {}
    for code in CODES:
        for row in csv.DictReader(open(LOCI / f"{code}.loci.tsv"), delimiter="\t"):
            status[row["locus_id"]] = row["status"]
            loci_of[row["orthogroup"]].setdefault(code, []).append(row["locus_id"])
    fams = {}
    for row in csv.DictReader(open(GC / "family_table.tsv"), delimiter="\t"):
        if row["grade"] not in ("G1", "G2"):
            continue
        members = {c: [l for l in loci_of[row["orthogroup"]].get(c, []) if status[l] != "disrupted"] for c in CODES}
        if all(len(m) == 1 and status[m[0]] == "intact" for m in members.values()):
            fams[row["orthogroup"]] = ({c: members[c][0] for c in CODES}, row["grade"])
    return fams


def p_distance(a, b):
    sites = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    return sum(x != y for x, y in sites) / len(sites) if sites else float("nan")


def run_family(args):
    og, members, cds = args
    seqs = {}
    for code, lid in members.items():
        s = cds.get(lid, "")
        s = s[:len(s) - len(s) % 3]
        prot = translate(s)
        if prot.endswith("*"):
            prot, s = prot[:-1], s[:-3]
        if not prot or "*" in prot or len(prot) < 50:
            return None
        seqs[code] = (s, prot)
    tmp = tempfile.mkdtemp(prefix=f"{og}_", dir=EVO / "tmp")
    try:
        with open(f"{tmp}/p.faa", "w") as fa:
            for code in CODES:
                write_fasta(fa, code, seqs[code][1])
        aln = subprocess.run([str(ENV_BIN / "mafft"), "--quiet", "--auto", f"{tmp}/p.faa"],
                             capture_output=True, text=True, check=True).stdout
        paln = {}
        for block in aln.split(">")[1:]:
            name, *lines = block.split("\n")
            paln[name.strip()] = "".join(lines).upper()
        codon = {}
        for code in CODES:
            nuc, i, out = seqs[code][0], 0, []
            for aa in paln[code]:
                if aa == "-":
                    out.append("---")
                else:
                    out.append(nuc[i:i + 3])
                    i += 3
            codon[code] = "".join(out)
        length = len(codon[CODES[0]])
        with open(f"{tmp}/c.phy", "w") as ph:
            ph.write(f" 4 {length}\n")
            for code in CODES:
                ph.write(f"{code}  {codon[code]}\n")
        with open(f"{tmp}/yn00.ctl", "w") as ctl:
            ctl.write("seqfile = c.phy\noutfile = yn\nverbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n")
        subprocess.run([str(ENV_BIN / "yn00"), "yn00.ctl"], cwd=tmp, capture_output=True, check=False)
        result = {"orthogroup": og, "aln_codons": length // 3}
        text = open(f"{tmp}/yn").read() if os.path.exists(f"{tmp}/yn") else ""
        # Yang & Nielsen (2000) table: "seq. seq. S N t kappa omega dN +- SE dS +- SE"
        section = text.split("(B) Yang & Nielsen (2000) method")[-1].split("(C)")[0]
        for line in section.splitlines():
            f = line.split()
            if len(f) >= 11 and f[0].isdigit() and f[1].isdigit():
                a, b = CODES[int(f[0]) - 1], CODES[int(f[1]) - 1]
                key = f"{min(a, b, key=CODES.index)}-{max(a, b, key=CODES.index)}"
                result[f"dN_{key}"], result[f"dS_{key}"], result[f"w_{key}"] = float(f[7]), float(f[10]), float(f[6])
        for a, b in PAIRS:
            result[f"ident_{a}-{b}"] = 1 - p_distance(paln[a], paln[b])
        d = {(a, b): p_distance(paln[a], paln[b]) for a, b in PAIRS}
        d.update({(b, a): v for (a, b), v in list(d.items())})
        A, B, C, D = CODES
        quartets = {f"{A}+{B}|{C}+{D}": d[A, B] + d[C, D], f"{A}+{C}|{B}+{D}": d[A, C] + d[B, D],
                    f"{A}+{D}|{B}+{C}": d[A, D] + d[B, C]}
        ordered = sorted(quartets.values())
        result["quartet"] = min(quartets, key=quartets.get) if ordered[0] < ordered[1] else "unresolved"
        result["protein_alignment"] = {c: paln[c] for c in CODES}
        return result
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ensure(EVO / "tmp")
    cds = extract_cds()
    fams = single_copy_families()
    print(f"[10] {len(fams)} single-copy intact families", flush=True)
    jobs = [(og, members, {l: cds.get(l, "") for l in members.values()}) for og, (members, _) in fams.items()]
    with ProcessPoolExecutor(THREADS) as pool:
        results = [r for r in pool.map(run_family, jobs, chunksize=20) if r]

    # Concatenated protein alignment for the species tree.
    concat = {c: [] for c in CODES}
    for r in results:
        for c in CODES:
            concat[c].append(r["protein_alignment"][c])
    with open(EVO / "concatenated_proteins.faa", "w") as fa:
        for c in CODES:
            write_fasta(fa, c, "".join(concat[c]))
    gdir = ensure(EVO / "gene_alignments")
    for r in results:
        with open(gdir / f"{r['orthogroup']}.faa", "w") as fa:
            for c in CODES:
                write_fasta(fa, c, r["protein_alignment"][c])

    if not (EVO / "species_tree.treefile").exists():
        subprocess.run([str(ENV_BIN / "iqtree"), "-s", str(EVO / "concatenated_proteins.faa"), "-m", "LG+F+R4",
                        "-B", "1000", "-T", str(min(THREADS, 8)), "--prefix", str(EVO / "species_tree"), "-redo"],
                       check=True, capture_output=True)

    keys = [k for k in results[0] if k not in ("protein_alignment",)]
    for r in results:
        for k in keys:
            r.setdefault(k, "")
    keys = ["orthogroup", "aln_codons", "quartet"] + sorted(k for k in set().union(*results) if k not in
                                                            ("orthogroup", "aln_codons", "quartet", "protein_alignment"))
    with open(EVO / "ortholog_evolution.tsv", "w") as out:
        out.write("\t".join(keys + ["grade", "mean_w", "conservation_class"]) + "\n")
        for r in results:
            ws = [r.get(f"w_{a}-{b}") for a, b in PAIRS]
            ds = [r.get(f"dS_{a}-{b}") for a, b in PAIRS]
            valid = [w for w, s in zip(ws, ds) if isinstance(w, float) and isinstance(s, float) and 0.005 < s < 2]
            mean_w = sum(valid) / len(valid) if valid else float("nan")
            r["mean_w"] = mean_w
            ident = min(r[f"ident_{a}-{b}"] for a, b in PAIRS)
            if not valid:
                cls = "undetermined"
            elif any(w > 1 for w in valid):
                cls = "dN/dS>1_in_some_pair"
            elif mean_w < 0.1 and ident >= 0.95:
                cls = "strongly_constrained"
            elif mean_w < 0.3:
                cls = "constrained"
            else:
                cls = "relaxed_or_fast"
            vals = [f"{r.get(k):.4f}" if isinstance(r.get(k), float) else str(r.get(k, "")) for k in keys]
            out.write("\t".join(vals + [fams[r["orthogroup"]][1], f"{mean_w:.4f}", cls]) + "\n")

    with open(EVO / "pair_summary.tsv", "w") as out:
        out.write("pair\tgenes\tmedian_protein_identity\tmedian_dN\tmedian_dS\tmedian_dN/dS\n")
        for a, b in PAIRS:
            k = f"{a}-{b}"
            ok = [r for r in results if isinstance(r.get(f"dS_{k}"), float) and 0 < r[f"dS_{k}"] < 2]
            out.write(f"{k}\t{len(ok)}\t{median(r[f'ident_{k}'] for r in results):.4f}\t"
                      f"{median(r[f'dN_{k}'] for r in ok):.4f}\t{median(r[f'dS_{k}'] for r in ok):.4f}\t"
                      f"{median(r[f'w_{k}'] for r in ok):.4f}\n")
    with open(EVO / "quartet_support.tsv", "w") as out:
        out.write("quartet\tgenes\n")
        for q, n in Counter(r["quartet"] for r in results).most_common():
            out.write(f"{q}\t{n}\n")
    shutil.rmtree(EVO / "tmp", ignore_errors=True)
    print(open(EVO / "pair_summary.tsv").read(), open(EVO / "quartet_support.tsv").read())


if __name__ == "__main__":
    main()
