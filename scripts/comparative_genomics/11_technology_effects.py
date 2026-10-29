#!/usr/bin/env python3
"""Step 11 - estimate how much of each between-genome difference comes from
sequencing / assembly technology rather than biology.

No reads are available, so technology effects are estimated from the assemblies
themselves, with signals that biology cannot (or can only rarely) produce:

1. Genome-wide indel decomposition (base-level error; weak at this divergence).
   Long-read assemblies that were not HiFi-based keep small indels in
   homopolymers. In each pairwise SyRI comparison, the rate of 1-bp indels
   inside homopolymers (>=4 identical bases) per aligned Mb is modelled as
       rate(a,b) = k * SNP_density(a,b) + e_a + e_b
   where k * SNP_density is the biological, divergence-proportional part and
   e_x is the excess contributed by assembly x. With 4 genomes there are 6
   pairs and 5 unknowns; they are solved by non-negative least squares. Indels
   outside homopolymers are fitted the same way. At 3-8% divergence real
   indels (~800 per aligned Mb) swamp assembly errors (typically tens per Mb),
   so this decomposition is reported but is not the primary estimate.
2. Frameshifts in conserved single-copy genes, query-controlled. The uniform
   projection used proteins from Adig, Amil and Aten, so an annotated genome
   receives near-identical "self" queries while Acer does not. Here only
   A. tenuis proteins of single-copy families are used and only Adig, Amil and
   Acer are scored; if A. tenuis is the outgroup (see species tree), it is
   equally distant from the three, so differences in frameshift/stop rates
   reflect the assemblies, not divergence. Rates are compared within
   alignment-identity strata (frameshift_by_identity.tsv).
3. False duplication: BUSCO duplicated %, highly similar (>=97%) duplicated segments.
4. Scaffolding: lineage-specific rearrangements (step 07) per unit of
   terminal branch length of the species tree (step 10). Biological
   rearrangements accumulate roughly with time, so an assembly with far more
   events per branch length than the others carries assembly errors; the share
   at assembly gaps shows how many are scaffolding joins.
5. Annotation pipeline: annotated genes vs uniform projection (step 05).

Outputs (11_technology_effects/):
  indel_by_pair.tsv            per-pair SNP and indel rates
  indel_decomposition.tsv      per-assembly excess (e_x) and fitted k
  technology_scorecard.tsv     all technology-sensitive metrics per assembly
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict

import numpy as np
from scipy.optimize import nnls

from config import CODES, OUT, ensure, read_fasta

WGA = OUT / "03_wga"
TECH = ensure(OUT / "11_technology_effects")
PAIRS = [(a, b) for i, a in enumerate(CODES) for b in CODES[i + 1:]]
HOMOPOLYMER = 4

# Sequencing / assembly strategy as documented by the data providers.
TECHNOLOGY = {
    "Adig": "haplotype-phased chromosome-scale assembly, reference haplotype (Takeuchi et al. 2025, OIST); platforms not stated in the source files",
    "Amil": "PacBio CLR + Illumina/10X (Canu/Falcon); chromosomes ordered with linkage maps (Fuller et al. 2020)",
    "Aten": "haplotype-phased chromosome-scale assembly, reference haplotype (Takeuchi et al. 2025, OIST); platforms not stated in the source files",
    "Acer": "PacBio HiFi 235x (hifiasm) + purge_dups + Arima Hi-C (YaHS) + manual curation (Darwin Tree of Life / ASG)",
}


def homopolymer_run(seq, pos, base):
    """Length of the run of `base` touching 0-based position pos."""
    left = pos
    while left > 0 and seq[left - 1] == base:
        left -= 1
    right = pos
    while right < len(seq) and seq[right] == base:
        right += 1
    return right - left


def count_pair(ref, qry):
    pair = f"{ref}_vs_{qry}"
    genome = dict(read_fasta(WGA / "chrom" / f"{ref}.fa"))
    snp = aligned = hp = non_hp = 0
    with open(WGA / pair / f"{pair}.syri.out") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            kind = f[10]
            # SNPs/indels are reported inside every aligned block (syntenic,
            # inverted, translocated, duplicated), so all of them form the denominator.
            if kind in ("SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDPAL"):
                aligned += int(f[2]) - int(f[1]) + 1
            elif kind == "SNP":
                snp += 1
            elif kind in ("INS", "DEL"):
                ref_allele, alt_allele = f[3], f[4]
                if abs(len(ref_allele) - len(alt_allele)) != 1:
                    continue
                seq = genome[f[0]]
                pos = int(f[1]) - 1  # SyRI: 1-based anchor base before the event
                if kind == "DEL":
                    base = ref_allele[-1].upper()
                    run = homopolymer_run(seq, pos + 1, base)
                else:
                    base = alt_allele[-1].upper()
                    # inserted base extends a run that starts after the anchor
                    run = 1 + (homopolymer_run(seq, pos + 1, base) if pos + 1 < len(seq) and seq[pos + 1] == base
                               else (homopolymer_run(seq, pos, base) if seq[pos] == base else 0))
                if run >= HOMOPOLYMER:
                    hp += 1
                else:
                    non_hp += 1
    return {"pair": pair, "aligned_Mb": aligned / 1e6, "snps": snp, "indel1_homopolymer": hp,
            "indel1_other": non_hp}


def decompose(rows, column):
    """Fit rate(a,b) = k*snp_density + e_a + e_b by non-negative least squares."""
    A, y = [], []
    for row in rows:
        a, b = row["pair"].split("_vs_")
        x = [row["snps"] / row["aligned_Mb"]] + [1.0 if c in (a, b) else 0.0 for c in CODES]
        A.append(x)
        y.append(row[column] / row["aligned_Mb"])
    coef, residual = nnls(np.array(A), np.array(y))
    fitted = np.array(A) @ coef
    return coef[0], dict(zip(CODES, coef[1:])), residual, fitted, y


STRATA = ((0.0, 0.90), (0.90, 0.95), (0.95, 1.01))


def controlled_frameshifts():
    """Frameshift/stop rate of A. tenuis single-copy proteins in the other genomes.

    Rates are computed within alignment-identity strata (>=90% query coverage):
    miniprot calls spurious frameshifts more often in divergent alignments, so
    genomes are compared at equal identity; the >=95% stratum is the headline.
    """
    from config import orthogroups_dir
    families = read_tsv(OUT / "09_gene_conservation" / "family_table.tsv")
    single = {r["orthogroup"] for r in families if r["grade"] in ("G1", "G2")}
    queries = set()
    for row in csv.DictReader(open(orthogroups_dir() / "Orthogroups.tsv"), delimiter="\t"):
        if row["Orthogroup"] in single:
            queries.update(g.strip() for g in row["Aten"].split(",") if g.strip())
    result, rows = {}, []
    for code in ("Adig", "Amil", "Acer"):
        best = {}
        with gzip.open(OUT / "04_orthology" / "miniprot" / f"{code}.miniprot.gff.gz", "rt") as handle:
            for line in handle:
                if not line.startswith("##PAF"):
                    continue
                f = line.rstrip("\n").split("\t")  # f[0] is "##PAF", so PAF column i is f[i+1]
                if f[1] not in queries:
                    continue
                tags = dict(t.split(":", 2)[::2] for t in f[13:] if t[:3] in ("ms:", "fs:", "st:"))
                hit = (int(tags["ms"]), int(tags["fs"]) > 0, int(tags["st"]) > 0,
                       int(f[10]) / int(f[11]), (int(f[4]) - int(f[3])) / int(f[2]))
                if f[1] not in best or hit[0] > best[f[1]][0]:
                    best[f[1]] = hit
        for lo, hi in STRATA:
            sel = [b for b in best.values() if lo <= b[3] < hi and b[4] >= 0.9]
            n = len(sel)
            row = {"species": code, "identity": f"{lo:.2f}-{min(hi, 1):.2f}", "genes": n,
                   "frameshift_pct": 100 * sum(b[1] for b in sel) / max(1, n),
                   "stop_pct": 100 * sum(b[2] for b in sel) / max(1, n)}
            rows.append(row)
            if lo == 0.95:
                result[code] = row
    with open(TECH / "frameshift_by_identity.tsv", "w") as out:
        out.write("species\tidentity\tgenes\tframeshift_pct\tstop_pct\n")
        for r in rows:
            out.write(f"{r['species']}\t{r['identity']}\t{r['genes']}\t{r['frameshift_pct']:.2f}\t{r['stop_pct']:.2f}\n")
    return result


def terminal_branches():
    """Terminal branch lengths of the (unrooted) IQ-TREE species tree."""
    path = OUT / "10_sequence_evolution" / "species_tree.treefile"
    if not path.exists():
        return {}
    text = open(path).read()
    return {m.group(1): float(m.group(2)) for m in re.finditer(r"(Adig|Amil|Aten|Acer):([0-9.eE-]+)", text)}


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if path.exists() else []


def main():
    rows = []
    for a, b in PAIRS:
        print(f"[11] counting indels {a} vs {b}", flush=True)
        rows.append(count_pair(a, b))
    with open(TECH / "indel_by_pair.tsv", "w") as out:
        out.write("pair\taligned_Mb\tSNPs_per_Mb\tindel1_homopolymer_per_Mb\tindel1_other_per_Mb\t"
                  "homopolymer_share_of_1bp_indels\n")
        for r in rows:
            total = r["indel1_homopolymer"] + r["indel1_other"]
            out.write(f"{r['pair']}\t{r['aligned_Mb']:.1f}\t{r['snps'] / r['aligned_Mb']:.1f}\t"
                      f"{r['indel1_homopolymer'] / r['aligned_Mb']:.2f}\t{r['indel1_other'] / r['aligned_Mb']:.2f}\t"
                      f"{r['indel1_homopolymer'] / max(1, total):.3f}\n")

    fits = {}
    with open(TECH / "indel_decomposition.tsv", "w") as out:
        out.write("indel_class\tk_per_SNP\t" + "\t".join(f"excess_{c}_per_Mb" for c in CODES) +
                  "\tresidual_norm\tobserved\tfitted\n")
        for column in ("indel1_homopolymer", "indel1_other"):
            k, excess, residual, fitted, observed = decompose(rows, column)
            fits[column] = excess
            out.write(f"{column}\t{k:.5f}\t" + "\t".join(f"{excess[c]:.2f}" for c in CODES) +
                      f"\t{residual:.3f}\t{','.join(f'{v:.2f}' for v in observed)}\t"
                      f"{','.join(f'{v:.2f}' for v in fitted)}\n")

    method = {r["species"]: r for r in read_tsv(OUT / "08_assembly_summary" / "assembly_methodology.tsv")}
    # Unique events (merged across comparisons) that only one genome carries.
    events = read_tsv(OUT / "07_structural_variants" / "sv_events_merged.tsv")
    lineage = defaultdict(lambda: [0, 0])
    for ev in events:
        derived = ev["derived_pattern"]
        if derived in CODES and ev["class"] in ("biological_lineage_specific", "artifact_scaffolding",
                                                 "uncertain_one_gap", "uncertain_error_prone_assembly"):
            lineage[derived][0] += 1
            lineage[derived][1] += ev["class"] == "artifact_scaffolding"

    frameshift = controlled_frameshifts()
    branches = terminal_branches()
    absences = defaultdict(dict)
    for row in read_tsv(OUT / "09_gene_conservation" / "family_table.tsv"):
        if row["grade"] in ("G5", "G6"):
            for c in CODES:
                if row[f"{c}_absence"]:
                    absences[c][row[f"{c}_absence"]] = absences[c].get(row[f"{c}_absence"], 0) + 1

    with open(TECH / "technology_scorecard.tsv", "w") as out:
        cols = ["species", "technology", "contig_N50_Mb", "gaps", "pct_in_chromosomes",
                "BUSCO_genome_duplicated", "near_identical_dup_Mb", "pct_disrupted_in_conserved_single_copy",
                "frameshift_pct_ident95_Aten_queries", "stop_pct_ident95_Aten_queries",
                "homopolymer_indel_excess_per_Mb", "other_indel_excess_per_Mb",
                "lineage_specific_rearrangements", "of_which_gap_bounded_artifacts",
                "terminal_branch_length", "lineage_rearrangements_per_branch_length",
                "gene_absences_at_assembly_gap", "gene_absences_syntenic_region_intact",
                "annotated_genes", "projected_loci_missing_from_annotation"]
        out.write("\t".join(cols) + "\n")
        for c in CODES:
            m = method.get(c, {})
            vals = [c, TECHNOLOGY[c], m.get("contig_N50_Mb", ""), m.get("gaps", ""), m.get("pct_in_chromosomes", ""),
                    m.get("BUSCO_genome_duplicated", ""), m.get("near_identical_dup_Mb", ""),
                    m.get("pct_disrupted_in_conserved_single_copy", ""),
                    f"{frameshift[c]['frameshift_pct']:.2f}" if c in frameshift else "NA (query species)",
                    f"{frameshift[c]['stop_pct']:.2f}" if c in frameshift else "NA (query species)",
                    f"{fits['indel1_homopolymer'][c]:.2f}", f"{fits['indel1_other'][c]:.2f}",
                    str(lineage[c][0]), str(lineage[c][1]),
                    f"{branches[c]:.4f}" if c in branches else "",
                    f"{lineage[c][0] / branches[c]:.0f}" if branches.get(c) else "",
                    str(absences[c].get("absent_gap", 0)), str(absences[c].get("absent_synteny_intact", 0)),
                    m.get("annotated_genes", ""),
                    m.get("projected_loci_missing_from_annotation", "")]
            out.write("\t".join(vals) + "\n")
    for name in ("indel_by_pair.tsv", "indel_decomposition.tsv", "technology_scorecard.tsv"):
        print(open(TECH / name).read())


if __name__ == "__main__":
    main()
