#!/usr/bin/env python3
"""Step 09 - gene conservation grades and duplications (annotation independent).

Copy numbers come from the uniform miniprot loci (step 05), not from the four
different annotations. For every gene family (orthogroup):

Copy number
  raw   = intact + partial loci
  corr  = raw minus copies lying in highly similar (>=97%, >=10 kb) duplicated
          assembly segments, which are candidate false duplications from
          uncollapsed haplotypes (at least one copy is kept)

Absence verification (family present elsewhere, corr = 0 in genome G)
  disrupted      a frameshifted/stop-containing copy exists (pseudogene or indel error)
  fragment       only a 30-50% coverage hit exists
  absent_gap     the syntenic region in G (between single-copy flanking genes)
                 contains an assembly gap: absence may be an assembly hole
  absent_synteny_intact  the flanking genes are collinear in G and no gap lies
                 between them: confident gene loss
  absent_unresolved      flanks not found in G (rearranged or unassembled)

Conservation grades
  G1 core single-copy, syntenic           1 copy in all 4, collinear in >=5 of 6 pairs
  G2 core single-copy, relocated          1 copy in all 4, not collinear in >=2 pairs
  G3 core multi-copy, stable copy number
  G4 core, copy-number variable
  G5 present in 2-3 genomes (lineage loss; see absence class)
  G6 single-genome family
  R  repeat-like family (>50 loci in a genome; transposon-derived, not graded)

Duplication mode of each extra copy, per genome (MCScanX precedence):
  segmental (in an intra-genome collinear block) > tandem (<=2 loci apart,
  <=100 kb) > proximal (<=20 loci apart) > dispersed; copies on unplaced
  scaffolds are "unplaced"; copies in near-identical segments are flagged
  assembly_dup_suspect.

Outputs (09_gene_conservation/): family_table.tsv, locus_duplication.tsv,
grade_summary.tsv, duplication_summary.tsv, copy_number_changes.tsv,
gene_grades_<species>.tsv (per annotated gene).
"""

from __future__ import annotations

import bisect
import csv
from collections import Counter, defaultdict

from config import ANNOTATED, CODES, INPUTS, OUT, ensure, orthogroups_dir

LOCI = OUT / "05_loci"
SYN = OUT / "06_synteny"
QC = OUT / "02_assembly_qc"
SUMDIR = OUT / "08_assembly_summary"
GC = ensure(OUT / "09_gene_conservation")
REPEAT_FAMILY = 50
PAIRS = [(a, b) for i, a in enumerate(CODES) for b in CODES[i + 1:]]


def read_intervals(path):
    data = defaultdict(list)
    if path.exists():
        for line in open(path):
            c, s, e = line.split()[:3]
            data[c].append((int(s), int(e)))
    for c in data:
        data[c].sort()
    return data


def overlap_frac(ivs, start, end):
    total = 0
    for s, e in ivs:
        if s >= end:
            break
        total += max(0, min(e, end) - max(s, start))
    return total / max(1, end - start)


def load():
    loci, by_og = {}, defaultdict(lambda: defaultdict(list))
    disrupted = defaultdict(Counter)
    for code in CODES:
        for row in csv.DictReader(open(LOCI / f"{code}.loci.tsv"), delimiter="\t"):
            if row["status"] == "disrupted":
                disrupted[row["orthogroup"]][code] += 1
                continue
            row["species"] = code
            row["start"], row["end"] = int(row["start"]), int(row["end"])
            loci[row["locus_id"]] = row
            by_og[row["orthogroup"]][code].append(row["locus_id"])
    families = {row["orthogroup"]: row for row in csv.DictReader(open(LOCI / "orthogroup_counts.tsv"), delimiter="\t")}
    return loci, by_og, disrupted, families


def synteny():
    """Inter-species anchor pairs and loci inside intra-genome blocks."""
    inter, intra = set(), set()
    for row in csv.DictReader(open(SYN / "anchors.tsv"), delimiter="\t"):
        if row["species_1"] == row["species_2"]:
            intra.update((row["locus_1"], row["locus_2"]))
        else:
            inter.add(frozenset((row["locus_1"], row["locus_2"])))
    return inter, intra


def gene_order(loci):
    """Rank of each locus along its sequence, per genome."""
    per_seq = defaultdict(list)
    for lid, row in loci.items():
        per_seq[(row["species"], row["chrom"])].append((row["start"], lid))
    rank = {}
    for key, items in per_seq.items():
        items.sort()
        for i, (_, lid) in enumerate(items):
            rank[lid] = i
    return rank


def duplication_modes(loci, by_og, intra, near_ident, rank):
    modes = {}
    for og, per in by_og.items():
        for code, members in per.items():
            if len(members) < 2:
                continue
            for lid in members:
                row = loci[lid]
                others = [loci[o] for o in members if o != lid]
                same = [o for o in others if o["chrom"] == row["chrom"]]
                gaps = [abs(rank[o["locus_id"]] - rank[lid]) for o in same]
                dist = [abs(o["start"] - row["start"]) for o in same]
                if lid in intra:
                    mode = "segmental"
                elif any(g <= 2 and d <= 100_000 for g, d in zip(gaps, dist)):
                    mode = "tandem"
                elif any(g <= 20 for g in gaps):
                    mode = "proximal"
                elif not row["chrom"].startswith("chr"):
                    mode = "unplaced"
                else:
                    mode = "dispersed"
                suspect = overlap_frac(near_ident[code].get(row["chrom"], []), row["start"], row["end"]) >= 0.5
                modes[lid] = (mode, suspect)
    return modes


class AnchorIndex:
    def __init__(self):
        self.pos = {c: {} for c in CODES}
        self.order = {c: defaultdict(list) for c in CODES}
        for row in csv.DictReader(open(SYN / "single_copy_anchors.tsv"), delimiter="\t"):
            for c in CODES:
                mid = (int(row[f"{c}_start"]) + int(row[f"{c}_end"])) // 2
                self.pos[c][row["orthogroup"]] = (row[f"{c}_chrom"], mid)
                self.order[c][row[f"{c}_chrom"]].append((mid, row["orthogroup"]))
        for c in CODES:
            for chrom in self.order[c]:
                self.order[c][chrom].sort()

    def flanks(self, code, chrom, pos, n=2, max_dist=1_000_000):
        items = self.order[code].get(chrom, [])
        i = bisect.bisect_left(items, (pos, ""))
        left = [og for p, og in items[max(0, i - n):i] if pos - p <= max_dist]
        right = [og for p, og in items[i:i + n] if p - pos <= max_dist]
        return left, right


def verify_absence(og, missing, by_og, loci, anchors, gaps, families, disrupted):
    if disrupted[og].get(missing):
        return "disrupted"
    if float(families[og][f"{missing}_weak_cov"]) >= 0.3:
        return "fragment"
    for code in CODES:
        members = by_og[og].get(code, [])
        if len(members) != 1 or code == missing:
            continue
        row = loci[members[0]]
        if not row["chrom"].startswith("chr"):
            continue
        left, right = anchors.flanks(code, row["chrom"], (row["start"] + row["end"]) // 2)
        if not left or not right:
            continue
        lpos = [anchors.pos[missing][og_] for og_ in left]
        rpos = [anchors.pos[missing][og_] for og_ in right]
        chroms = {p[0] for p in lpos + rpos}
        if len(chroms) != 1:
            return "absent_unresolved"
        chrom = chroms.pop()
        # The expected position is between the nearest flanking genes; only a
        # gap inside that interval can hide the gene.
        lo, hi = sorted((lpos[-1][1], rpos[0][1]))
        if hi - lo > 1_000_000:
            return "absent_unresolved"
        if any(s < hi and lo < e for s, e in gaps[missing].get(chrom, [])):
            return "absent_gap"
        return "absent_synteny_intact"
    return "absent_unresolved"


def main():
    loci, by_og, disrupted, families = load()
    inter, intra = synteny()
    rank = gene_order(loci)
    near_ident = {c: read_intervals(SUMDIR / f"{c}.near_identical_segments.bed") for c in CODES}
    gaps = {c: read_intervals(INPUTS / f"{c}.gaps.bed") for c in CODES}
    repeats = {c: read_intervals(QC / "red" / f"{c}.bed") for c in CODES}
    modes = duplication_modes(loci, by_og, intra, near_ident, rank)
    anchors = AnchorIndex()

    with open(GC / "locus_duplication.tsv", "w") as out:
        out.write("locus_id\tspecies\tchrom\tstart\tend\torthogroup\tmode\tassembly_dup_suspect\n")
        for lid, (mode, suspect) in sorted(modes.items()):
            r = loci[lid]
            out.write(f"{lid}\t{r['species']}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{r['orthogroup']}\t"
                      f"{mode}\t{int(suspect)}\n")

    grade_of, rows = {}, []
    for og, fam in families.items():
        per = by_og.get(og, {})
        raw = {c: len(per.get(c, [])) for c in CODES}
        if max(raw.values()) == 0:
            continue
        suspect = {c: sum(modes.get(l, ("", False))[1] for l in per.get(c, [])) for c in CODES}
        corr = {c: max(raw[c] - suspect[c], 1) if raw[c] else 0 for c in CODES}
        present = [c for c in CODES if corr[c] > 0]
        absence = {c: verify_absence(og, c, by_og, loci, anchors, gaps, families, disrupted)
                   for c in CODES if c not in present}
        syntenic_pairs = ""
        if max(raw.values()) > REPEAT_FAMILY:
            grade = "R"
        elif len(present) == 4 and all(corr[c] == 1 for c in CODES):
            single = {c: per[c][0] for c in CODES} if all(raw[c] == 1 for c in CODES) else None
            n_syn = sum(frozenset((single[a], single[b])) in inter for a, b in PAIRS) if single else 0
            syntenic_pairs = n_syn
            grade = "G1" if n_syn >= 5 else "G2"
        elif len(present) == 4:
            grade = "G3" if len(set(corr.values())) == 1 else "G4"
        elif len(present) >= 2:
            grade = "G5"
        else:
            grade = "G6"
        # Single-genome families from one annotation that look like transposon
        # ORFs or short spurious models are flagged rather than trusted.
        flag = ""
        if grade == "G6":
            sp = present[0]
            lrow = loci[per[sp][0]]
            rep = overlap_frac(repeats[sp].get(lrow["chrom"], []), lrow["start"], lrow["end"])
            if rep >= 0.5:
                flag = "repeat_overlap>=50%"
        grade_of[og] = grade
        rows.append({"orthogroup": og, "grade": grade, **{f"{c}_raw": raw[c] for c in CODES},
                     **{f"{c}_corr": corr[c] for c in CODES},
                     **{f"{c}_disrupted": disrupted[og].get(c, 0) for c in CODES},
                     **{f"{c}_absence": absence.get(c, "") for c in CODES},
                     "syntenic_pairs_of_6": syntenic_pairs, "flag": flag,
                     **{f"{a}_annotated": fam[f"{a}_annotated"] for a in ANNOTATED},
                     "product": fam["product"]})

    with open(GC / "family_table.tsv", "w") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Grade summary, and losses split by how well the absence is supported.
    with open(GC / "grade_summary.tsv", "w") as out:
        out.write("grade\tfamilies\n")
        for grade, n in sorted(Counter(r["grade"] for r in rows).items()):
            out.write(f"{grade}\t{n}\n")
        out.write("\nspecies\tabsence_class\tfamilies\n")
        for c in CODES:
            for cls, n in sorted(Counter(r[f"{c}_absence"] for r in rows
                                         if r[f"{c}_absence"] and r["grade"] in ("G5", "G6")).items()):
                out.write(f"{c}\t{cls}\t{n}\n")

    with open(GC / "duplication_summary.tsv", "w") as out:
        out.write("species\tmode\tloci\tassembly_dup_suspect\n")
        agg = defaultdict(lambda: [0, 0])
        for lid, (mode, suspect) in modes.items():
            if grade_of.get(loci[lid]["orthogroup"]) == "R":
                continue
            agg[(loci[lid]["species"], mode)][0] += 1
            agg[(loci[lid]["species"], mode)][1] += suspect
        for (sp, mode), (n, s) in sorted(agg.items()):
            out.write(f"{sp}\t{mode}\t{n}\t{s}\n")

    # Lineage-specific copy-number changes on corrected counts.
    with open(GC / "copy_number_changes.tsv", "w") as out:
        out.write("orthogroup\tspecies\tchange\t" + "\t".join(f"{c}_corr" for c in CODES) +
                  "\tdominant_mode\tproduct\n")
        for r in rows:
            if r["grade"] not in ("G4", "G5"):
                continue
            corr = {c: r[f"{c}_corr"] for c in CODES}
            for c in CODES:
                others = [corr[o] for o in CODES if o != c]
                change = ""
                if corr[c] >= 2 * max(others) and corr[c] - max(others) >= 2 and min(others) >= 1:
                    change = "expansion"
                elif corr[c] >= 1 and 2 * corr[c] <= min(others) and min(others) - corr[c] >= 2:
                    change = "contraction"
                if change:
                    m = Counter(modes[l][0] for l in by_og[r["orthogroup"]].get(c, []) if l in modes)
                    dom = m.most_common(1)[0][0] if m else ""
                    out.write(f"{r['orthogroup']}\t{c}\t{change}\t" + "\t".join(str(corr[x]) for x in CODES) +
                              f"\t{dom}\t{r['product']}\n")

    # Per annotated gene: its family's grade (so users can look up any gene).
    for code in ANNOTATED:
        gene_og = {}
        results = orthogroups_dir()
        for name in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
            reader = csv.DictReader(open(results / name), delimiter="\t")
            for row in reader:
                for g in filter(None, (x.strip() for x in row.get(code, "").split(","))):
                    gene_og[g] = row["Orthogroup"]
        with open(GC / f"gene_grades_{code}.tsv", "w") as out:
            out.write("gene_id\torthogroup\tgrade\tproduct\n")
            for row in csv.DictReader(open(INPUTS / "annotations" / f"{code}.genes.tsv"), delimiter="\t"):
                og = gene_og.get(row["gene_id"], "")
                out.write(f"{row['gene_id']}\t{og}\t{grade_of.get(og, 'not_projected')}\t{row['product']}\n")
    print(open(GC / "grade_summary.tsv").read())


if __name__ == "__main__":
    main()
