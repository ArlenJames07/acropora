#!/usr/bin/env python3
"""Step 06 - gene-order synteny on the uniform loci (annotation independent).

Homology = two intact/partial loci assigned to the same orthogroup. Families
with more than MAX_COPIES loci in any genome (mostly transposon-derived) are
left out because they cannot anchor positions. MCScanX then finds collinear
blocks both between species (orthologous synteny) and within each genome
(segmental duplications).

Outputs (06_synteny/):
  acropora.gff / acropora.blast       MCScanX inputs (homology pairs in BLAST m8 layout)
  acropora.collinearity               MCScanX output
  blocks.tsv                          one row per collinear block
  anchors.tsv                         one row per anchor pair
  single_copy_anchors.tsv             positions of 1:1:1:1 families in all genomes
"""

from __future__ import annotations

import re
import subprocess
from collections import defaultdict
from itertools import combinations

from config import CODES, ENV_BIN, OUT, ensure

LOCI = OUT / "05_loci"
SYN = ensure(OUT / "06_synteny")
MAX_COPIES = 10
PREFIX = {"Adig": "ad", "Amil": "am", "Aten": "at", "Acer": "ac"}


def load_loci():
    loci = {}
    for code in CODES:
        with open(LOCI / f"{code}.loci.tsv") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                if row["status"] == "disrupted":
                    continue
                row["species"] = code
                row["start"], row["end"] = int(row["start"]), int(row["end"])
                loci[row["locus_id"]] = row
    return loci


def mcscanx_chrom(code, chrom, scaffolds):
    if chrom.startswith("chr"):
        return f"{PREFIX[code]}{chrom[3:]}"
    # Unplaced scaffolds get numbers above 100 so they never clash with chr1-14.
    key = (code, chrom)
    scaffolds.setdefault(key, 100 + len([k for k in scaffolds if k[0] == code]))
    return f"{PREFIX[code]}{scaffolds[key]}"


def run_mcscanx():
    loci = load_loci()
    by_og = defaultdict(list)
    for lid, row in loci.items():
        by_og[row["orthogroup"]].append(lid)

    scaffolds = {}
    with open(SYN / "acropora.gff", "w") as gff:
        for lid, row in loci.items():
            gff.write(f"{mcscanx_chrom(row['species'], row['chrom'], scaffolds)}\t{lid}\t"
                      f"{row['start']}\t{row['end']}\n")

    with open(SYN / "acropora.blast", "w") as out:
        for og, members in by_og.items():
            per_sp = defaultdict(int)
            for lid in members:
                per_sp[loci[lid]["species"]] += 1
            if max(per_sp.values()) > MAX_COPIES:
                continue
            for a, b in combinations(members, 2):
                # MCScanX expects BLAST m8 columns; only the pair and e-value matter.
                out.write(f"{a}\t{b}\t100\t100\t0\t0\t1\t100\t1\t100\t1e-100\t500\n")

    subprocess.run([str(ENV_BIN / "MCScanX"), "-b", "0", str(SYN / "acropora")],
                   check=True, stdout=open(SYN / "mcscanx.log", "w"), stderr=subprocess.STDOUT)


def parse_collinearity():
    loci = load_loci()
    blocks, anchors, block = [], [], None
    pattern = re.compile(r"## Alignment (\d+): score=([\d.]+) e_value=(\S+) N=(\d+) (\S+)&(\S+) (plus|minus)")
    with open(SYN / "acropora.collinearity") as handle:
        for line in handle:
            match = pattern.match(line)
            if match:
                block = {"block": int(match.group(1)), "n": int(match.group(4)),
                         "orientation": match.group(7), "pairs": []}
                blocks.append(block)
            elif not line.startswith("#") and line.strip() and block is not None:
                fields = line.split("\t")
                block["pairs"].append((fields[1].strip(), fields[2].strip()))

    with open(SYN / "blocks.tsv", "w") as bt, open(SYN / "anchors.tsv", "w") as at:
        bt.write("block\tspecies_1\tchrom_1\tstart_1\tend_1\tspecies_2\tchrom_2\tstart_2\tend_2\t"
                 "anchors\torientation\ttype\n")
        at.write("block\tlocus_1\tlocus_2\tspecies_1\tspecies_2\torthogroup\n")
        for b in blocks:
            a_loci = [loci[p[0]] for p in b["pairs"]]
            b_loci = [loci[p[1]] for p in b["pairs"]]
            s1, s2 = a_loci[0]["species"], b_loci[0]["species"]
            kind = "intra" if s1 == s2 else "inter"
            bt.write(f"{b['block']}\t{s1}\t{a_loci[0]['chrom']}\t{min(l['start'] for l in a_loci)}\t"
                     f"{max(l['end'] for l in a_loci)}\t{s2}\t{b_loci[0]['chrom']}\t"
                     f"{min(l['start'] for l in b_loci)}\t{max(l['end'] for l in b_loci)}\t"
                     f"{b['n']}\t{b['orientation']}\t{kind}\n")
            for x, y in zip(a_loci, b_loci):
                at.write(f"{b['block']}\t{x['locus_id']}\t{y['locus_id']}\t{s1}\t{s2}\t{x['orthogroup']}\n")

    # Families with exactly one intact/partial locus in every genome: the
    # position anchors used to test rearrangements and gene absences.
    by_og = defaultdict(dict)
    multi = set()
    for lid, row in loci.items():
        og, sp = row["orthogroup"], row["species"]
        if sp in by_og[og]:
            multi.add(og)
        by_og[og][sp] = row
    with open(SYN / "single_copy_anchors.tsv", "w") as out:
        out.write("orthogroup\t" + "\t".join(f"{c}_locus\t{c}_chrom\t{c}_start\t{c}_end\t{c}_strand"
                                             for c in CODES) + "\n")
        for og, rows in sorted(by_og.items()):
            if og in multi or len(rows) != len(CODES):
                continue
            out.write(og + "\t" + "\t".join(
                f"{rows[c]['locus_id']}\t{rows[c]['chrom']}\t{rows[c]['start']}\t{rows[c]['end']}\t"
                f"{rows[c]['strand']}" for c in CODES) + "\n")


def main():
    if not (SYN / "acropora.collinearity").exists():
        run_mcscanx()
    parse_collinearity()


if __name__ == "__main__":
    main()
