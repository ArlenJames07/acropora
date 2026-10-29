#!/usr/bin/env python3
"""Step 05 - resolve miniprot projections into one uniform locus set per genome.

Every annotated protein (A. digitifera, A. millepora, A. tenuis) was aligned to
all four genomes with identical miniprot settings (step 04). Here, per genome:

  * hits covering >=50% of the query protein are processed from best to worst
    score; a hit whose CDS overlaps an already accepted locus becomes supporting
    evidence for that locus, otherwise it founds a new locus
  * each locus inherits the OrthoFinder orthogroup (gene family) of its best hit
  * status: intact (>=80% query coverage, no frameshift, no in-frame stop),
    partial (50-80% coverage, no frameshift/stop) or disrupted (frameshift or
    premature stop: pseudogene-like, or an indel error in the assembly)
  * weaker hits (30-50% coverage) are kept only as evidence that a family is
    not completely absent from a genome

Outputs (05_loci/):
  <sp>.loci.tsv           one row per locus
  <sp>.loci.gff3          CDS models of the loci (for CDS extraction)
  <sp>.loci.faa           translated locus proteins
  orthogroup_counts.tsv   copies per family per genome + annotation counts
  weak_evidence.tsv       best sub-threshold hit per family and genome
  annotation_vs_projection.tsv
"""

from __future__ import annotations

import gzip
from collections import Counter, defaultdict

from config import ANNOTATED, CODES, INPUTS, OUT, ensure, orthogroups_dir, write_fasta

ORTHO = OUT / "04_orthology"
LOCI = ensure(OUT / "05_loci")
MIN_COV, INTACT_COV, WEAK_COV = 0.5, 0.8, 0.3
BIN = 50_000


def load_orthogroups():
    gene_og, og_members = {}, defaultdict(list)
    results = orthogroups_dir()
    for name in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
        with open(f"{results}/{name}") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                og = fields[0]
                for column, cell in zip(header[1:], fields[1:]):
                    for gene in filter(None, (g.strip() for g in cell.split(","))):
                        gene_og[gene] = og
                        og_members[og].append(gene)
    return gene_og, og_members


def load_products():
    products = {}
    for code in ANNOTATED:
        with open(INPUTS / "annotations" / f"{code}.genes.tsv") as handle:
            next(handle)
            for line in handle:
                f = line.rstrip("\n").split("\t")
                products[f[0]] = f[7]
    return products


def parse_miniprot(path):
    """Yield one dict per alignment (PAF line + translated protein + GFF lines)."""
    hit = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##PAF"):
                if hit and "chrom" in hit:
                    yield hit
                f = line.rstrip("\n").split("\t")
                tags = dict(t.split(":", 2)[::2] for t in f[13:] if t[:3] in ("ms:", "fs:", "st:"))
                hit = {"query": f[1], "qlen": int(f[2]), "qcov": (int(f[4]) - int(f[3])) / int(f[2]),
                       "score": int(tags["ms"]), "fs": int(tags["fs"]), "st": int(tags["st"]),
                       "cds": [], "protein": ""}
            elif line.startswith("##STA"):
                hit["protein"] = line.rstrip("\n").split("\t")[1]
            elif not line.startswith("#") and hit is not None:
                f = line.rstrip("\n").split("\t")
                if f[2] == "mRNA":
                    attrs = dict(a.split("=", 1) for a in f[8].split(";"))
                    hit.update(chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                               identity=float(attrs["Identity"]))
                elif f[2] == "CDS":
                    hit["cds"].append((int(f[3]), int(f[4]), f[7]))
    if hit and "chrom" in hit:
        yield hit


def overlaps(cds_a, cds_b):
    return any(a0 <= b1 and b0 <= a1 for a0, a1, _ in cds_a for b0, b1, _ in cds_b)


def resolve(code, gene_og):
    hits, weak = [], defaultdict(float)
    for hit in parse_miniprot(ORTHO / "miniprot" / f"{code}.miniprot.gff.gz"):
        og = gene_og.get(hit["query"])
        if og is None:
            continue
        hit["og"] = og
        if hit["qcov"] >= MIN_COV:
            hits.append(hit)
        elif hit["qcov"] >= WEAK_COV:
            weak[og] = max(weak[og], hit["qcov"])
    hits.sort(key=lambda h: (-h["score"], h["query"]))

    loci, index = [], defaultdict(list)  # (chrom, bin) -> locus numbers
    for hit in hits:
        bins = range(hit["start"] // BIN, hit["end"] // BIN + 1)
        candidates = {i for b in bins for i in index[(hit["chrom"], b)]}
        owner = next((i for i in sorted(candidates) if overlaps(hit["cds"], loci[i]["cds"])), None)
        if owner is None:
            hit["support"] = Counter({hit["og"]: 1})
            hit["n_hits"] = 1
            for b in bins:
                index[(hit["chrom"], b)].append(len(loci))
            loci.append(hit)
        else:
            loci[owner]["support"][hit["og"]] += 1
            loci[owner]["n_hits"] += 1
    return loci, weak


def status(locus):
    if locus["fs"] or locus["st"]:
        return "disrupted"
    return "intact" if locus["qcov"] >= INTACT_COV else "partial"


def compare_annotation(code, loci):
    """How the species' own annotation and the uniform projection agree."""
    genes = defaultdict(list)
    with open(INPUTS / "annotations" / f"{code}.genes.tsv") as handle:
        next(handle)
        for line in handle:
            f = line.split("\t")
            genes[f[1]].append((int(f[2]), int(f[3])))
    good = [l for l in loci if l["status"] != "disrupted"]
    spans = defaultdict(list)
    for l in good:
        spans[l["chrom"]].append((l["start"], l["end"]))

    def covered(intervals, chrom, start, end):
        return any(s <= end and start <= e for s, e in intervals.get(chrom, ()))

    # Bin intervals so the overlap test stays fast.
    def binned(table):
        out = defaultdict(list)
        for chrom, ivs in table.items():
            for s, e in ivs:
                for b in range(s // BIN, e // BIN + 1):
                    out[(chrom, b)].append((s, e))
        return out

    gene_bins, locus_bins = binned(genes), binned(spans)
    hit = lambda bins, c, s, e: any(covered({c: bins[(c, b)]}, c, s, e) for b in range(s // BIN, e // BIN + 1))
    n_genes = sum(len(v) for v in genes.values())
    recovered = sum(hit(locus_bins, c, s, e) for c, ivs in genes.items() for s, e in ivs)
    missed = sum(not hit(gene_bins, l["chrom"], l["start"], l["end"]) for l in good)
    return f"{code}\t{n_genes}\t{recovered}\t{len(good)}\t{missed}\n"


def chrom_key(chrom):
    return (0, int(chrom[3:]), "") if chrom.startswith("chr") else (1, 0, chrom)


def main():
    gene_og, og_members = load_orthogroups()
    products = load_products()
    counts = defaultdict(lambda: defaultdict(int))
    weak_all = {}
    with open(LOCI / "annotation_vs_projection.tsv", "w") as cmp:
        cmp.write("species\tannotated_genes\tannotated_genes_overlapping_a_locus\t"
                  "loci_intact_or_partial\tloci_not_covered_by_annotation\n")
        for code in CODES:
            print(f"[05] {code}", flush=True)
            loci, weak = resolve(code, gene_og)
            weak_all[code] = weak
            loci.sort(key=lambda h: (chrom_key(h["chrom"]), h["start"]))
            with open(LOCI / f"{code}.loci.tsv", "w") as table, \
                    open(LOCI / f"{code}.loci.gff3", "w") as gff, \
                    open(LOCI / f"{code}.loci.faa", "w") as faa:
                table.write("locus_id\tchrom\tstart\tend\tstrand\torthogroup\tstatus\tbest_query\t"
                            "query_coverage\tidentity\tframeshifts\tstops\tscore\tn_hits\t"
                            "orthogroups_supporting\n")
                gff.write("##gff-version 3\n")
                for n, locus in enumerate(loci, 1):
                    lid = f"{code}L{n:06d}"
                    locus["id"], locus["status"] = lid, status(locus)
                    counts[locus["og"]][(code, locus["status"])] += 1
                    support = ",".join(f"{og}:{c}" for og, c in locus["support"].most_common(3))
                    table.write(f"{lid}\t{locus['chrom']}\t{locus['start']}\t{locus['end']}\t"
                                f"{locus['strand']}\t{locus['og']}\t{locus['status']}\t{locus['query']}\t"
                                f"{locus['qcov']:.3f}\t{locus['identity']:.4f}\t{locus['fs']}\t"
                                f"{locus['st']}\t{locus['score']}\t{locus['n_hits']}\t{support}\n")
                    gff.write(f"{locus['chrom']}\tminiprot\tmRNA\t{locus['start']}\t{locus['end']}\t.\t"
                              f"{locus['strand']}\t.\tID={lid};orthogroup={locus['og']}\n")
                    for start, end, phase in locus["cds"]:
                        gff.write(f"{locus['chrom']}\tminiprot\tCDS\t{start}\t{end}\t.\t"
                                  f"{locus['strand']}\t{phase}\tParent={lid}\n")
                    write_fasta(faa, lid, locus["protein"].rstrip("*"))
            if code in ANNOTATED:
                cmp.write(compare_annotation(code, loci))

    with open(LOCI / "orthogroup_counts.tsv", "w") as out:
        cols = [f"{c}_{s}" for c in CODES for s in ("intact", "partial", "disrupted")]
        out.write("orthogroup\t" + "\t".join(f"{c}_annotated" for c in ANNOTATED) + "\t" +
                  "\t".join(cols) + "\t" + "\t".join(f"{c}_weak_cov" for c in CODES) + "\tproduct\n")
        for og in sorted(og_members):
            members = og_members[og]
            ann = [sum(m.startswith(c + "_") for m in members) for c in ANNOTATED]
            names = Counter(products.get(m, "") for m in members if products.get(m, ""))
            informative = [n for n, _ in names.most_common() if not n.startswith("uncharacterized")]
            product = informative[0] if informative else (names.most_common(1)[0][0] if names else "")
            row = [str(counts[og][(c, s)]) for c in CODES for s in ("intact", "partial", "disrupted")]
            weak = [f"{weak_all[c].get(og, 0):.2f}" for c in CODES]
            out.write(f"{og}\t" + "\t".join(map(str, ann)) + "\t" + "\t".join(row) + "\t" +
                      "\t".join(weak) + f"\t{product}\n")


if __name__ == "__main__":
    main()
