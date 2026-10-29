#!/usr/bin/env python3
"""Step 08 - quantify methodological differences between the assemblies.

These numbers are the baseline against which biological differences are judged:
  * BUSCO (cnidaria_odb12) completeness and duplication, genome and protein mode
  * uniform repeat content (Red), so genome-size differences can be split into
    repeat and non-repeat sequence
  * highly similar intra-genome duplications (>=97% identity, >=10 kb, from 20 kb
    genome fragments mapped back to their own genome). The threshold admits
    retained haplotypes, which differ by the ~1-3% heterozygosity of corals: the signature of retained haplotypes / false
    duplications rather than biological segmental duplications
  * which duplicated BUSCO genes fall in those near-identical segments
  * fraction of projected gene loci disrupted by frameshifts or stops: in
    conserved single-copy genes this mostly measures base-level indel errors
  * agreement between each annotation and the uniform projection

Output: 08_assembly_summary/assembly_methodology.tsv (+ helper tables)
"""

from __future__ import annotations

import csv
import glob
import gzip
import json
from collections import defaultdict

from config import ASSEMBLY_METHOD, CODES, INPUTS, OUT, ensure

QC = OUT / "02_assembly_qc"
LOCI = OUT / "05_loci"
SUM = ensure(OUT / "08_assembly_summary")
MIN_LEN, MIN_IDENTITY = 10_000, 0.97


def busco(mode, code):
    files = glob.glob(str(QC / "busco" / f"{mode}_{code}" / "short_summary*.json"))
    if not files:
        return {}
    res = json.load(open(files[0]))["results"]
    return {"complete": res.get("Complete percentage"), "single": res.get("Single copy percentage"),
            "duplicated": res.get("Multi copy percentage"), "fragmented": res.get("Fragmented percentage"),
            "missing": res.get("Missing percentage"), "n": res.get("n_markers")}


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def near_identical_segments(code):
    """Merged near-identical duplicated segments (both copies) per sequence."""
    segs = defaultdict(list)
    pairs = []
    path = QC / "self" / f"{code}.self.paf.gz"
    if not path.exists():
        return {}, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            f = line.split("\t")
            # Query = 20 kb fragment "<seq>_sliding:<start>-<end>" (1-based); convert
            # its alignment to coordinates on the sequence it was cut from.
            name, window = f[0].rsplit("_sliding:", 1)
            offset = int(window.split("-")[0]) - 1
            qn, qs, qe = name, offset + int(f[2]), offset + int(f[3])
            tn, ts, te = f[5], int(f[7]), int(f[8])
            matches, length = int(f[9]), int(f[10])
            if length < MIN_LEN or matches / length < MIN_IDENTITY:
                continue
            if qn == tn and qs < te and ts < qe:
                continue  # the fragment aligned back to its own origin
            segs[qn].append((qs, qe))
            segs[tn].append((ts, te))
            pairs.append((qn, qs, qe, tn, ts, te, matches / length))
    merged = {name: merge(v) for name, v in segs.items()}
    with open(SUM / f"{code}.near_identical_segments.bed", "w") as out:
        for name, ivs in sorted(merged.items()):
            for s, e in ivs:
                out.write(f"{name}\t{s}\t{e}\n")
    return merged, pairs


def overlaps(merged, name, start, end):
    return any(s < end and start < e for s, e in merged.get(name, ()))


def repeat_bp(code):
    total = 0
    path = QC / "red" / f"{code}.bed"
    if path.exists():
        with open(path) as handle:
            for line in handle:
                _, s, e = line.split()[:3]
                total += int(e) - int(s)
    return total


def main():
    stats = {row["species"]: row for row in csv.DictReader(open(INPUTS / "assembly_stats.tsv"), delimiter="\t")}
    ann = {}
    if (LOCI / "annotation_vs_projection.tsv").exists():
        ann = {row["species"]: row for row in csv.DictReader(open(LOCI / "annotation_vs_projection.tsv"), delimiter="\t")}

    rows = []
    for code in CODES:
        s = stats[code]
        total = int(s["total_bp"])
        merged, pairs = near_identical_segments(code)
        dup_bp = sum(e - st for ivs in merged.values() for st, e in ivs)
        dup_unplaced = sum(e - st for n, ivs in merged.items() if not n.startswith("chr") for st, e in ivs)
        rep = repeat_bp(code)

        g, p = busco("genome", code), busco("protein", code)
        dup_in_segments, dup_total = 0, 0
        table = glob.glob(str(QC / "busco" / f"genome_{code}" / "run_*" / "full_table.tsv"))
        if table:
            for line in open(table[0]):
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[1] == "Duplicated":
                    dup_total += 1
                    dup_in_segments += overlaps(merged, f[2], int(f[3]), int(f[4]))

        status = defaultdict(int)
        sc_disrupted = sc_total = 0
        lpath = LOCI / f"{code}.loci.tsv"
        if lpath.exists():
            for row in csv.DictReader(open(lpath), delimiter="\t"):
                status[row["status"]] += 1
        # Disruption among genes that are intact in the other three genomes is
        # the cleanest proxy for base-level (indel) assembly errors.
        if (LOCI / "orthogroup_counts.tsv").exists():
            for row in csv.DictReader(open(LOCI / "orthogroup_counts.tsv"), delimiter="\t"):
                others = [c for c in CODES if c != code]
                if all(row[f"{c}_intact"] == "1" and row[f"{c}_partial"] == "0" and row[f"{c}_disrupted"] == "0"
                       for c in others):
                    n_here = int(row[f"{code}_intact"]) + int(row[f"{code}_partial"]) + int(row[f"{code}_disrupted"])
                    if n_here == 1:
                        sc_total += 1
                        sc_disrupted += row[f"{code}_disrupted"] == "1"

        rows.append({
            "species": code, "method": ASSEMBLY_METHOD[code],
            "total_Mb": total / 1e6, "pct_in_chromosomes": float(s["pct_in_chromosomes"]),
            "unplaced_Mb": int(s["unplaced_bp"]) / 1e6, "contig_N50_Mb": int(s["contig_N50"]) / 1e6,
            "gaps": int(s["gaps"]),
            "repeat_Mb_Red": rep / 1e6, "repeat_pct_Red": 100 * rep / total,
            "nonrepeat_Mb": (total - rep) / 1e6,
            "near_identical_dup_Mb": dup_bp / 1e6, "near_identical_dup_on_unplaced_Mb": dup_unplaced / 1e6,
            "near_identical_pairs": len(pairs),
            "BUSCO_genome_complete": g.get("complete"), "BUSCO_genome_duplicated": g.get("duplicated"),
            "BUSCO_genome_fragmented": g.get("fragmented"), "BUSCO_genome_missing": g.get("missing"),
            "BUSCO_dup_in_near_identical_segments": f"{dup_in_segments}/{dup_total}",
            "BUSCO_protein_complete": p.get("complete"), "BUSCO_protein_duplicated": p.get("duplicated"),
            "annotated_genes": s["annotated_genes"],
            "loci_intact": status["intact"], "loci_partial": status["partial"], "loci_disrupted": status["disrupted"],
            "pct_disrupted_in_conserved_single_copy": 100 * sc_disrupted / sc_total if sc_total else None,
            "annotation_genes_overlapping_projection": ann.get(code, {}).get("annotated_genes_overlapping_a_locus", ""),
            "projected_loci_missing_from_annotation": ann.get(code, {}).get("loci_not_covered_by_annotation", ""),
        })

    with open(SUM / "assembly_methodology.tsv", "w") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.2f}" if isinstance(v, float) else v) for k, v in row.items()})
    print(open(SUM / "assembly_methodology.tsv").read())


if __name__ == "__main__":
    main()
