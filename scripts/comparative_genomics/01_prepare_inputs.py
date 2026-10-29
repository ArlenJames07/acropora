#!/usr/bin/env python3
"""Step 01 - build uniform inputs for every downstream comparison.

For each assembly:
  * genome FASTA with homologous chromosomes renamed chr1..chr14 in A. digitifera
    numbering and orientation (table from scripts/prepare_chromosomes.py);
    unplaced scaffolds are kept under their accession so no sequence is lost;
    soft-masking is removed because each assembly was masked differently
  * BED files of assembly gaps (N runs) and contigs: breakpoints that fall on
    these are candidate scaffolding artifacts
  * for annotated assemblies: longest-isoform proteins and a gene table whose
    coordinates (CDS span, so UTR-calling differences are ignored) are lifted
    into the normalized chromosome frame
"""

from __future__ import annotations

import re
from urllib.parse import unquote

from config import (ANNOTATED, CHROMOSOME_MAPS, INPUTS, RESOURCES, GENOMES,
                    REVERSE_COMPLEMENT, SPECIES, ensure, read_fasta, write_fasta,
                    open_text)
from prepare_chromosomes import COMPLEMENT, chromosome_number

GAP = re.compile(r"[Nn]+")


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    half, total = sum(lengths) / 2, 0
    for length in lengths:
        total += length
        if total >= half:
            return length
    return 0


def prepare_genome(code, key, genome_name):
    """Write normalized genome; return {source_seqid: (new_name, reverse, length)}."""
    out_dir = ensure(INPUTS / "genomes")
    fasta = out_dir / f"{code}.fa"
    seq_map, stats = {}, {"lengths": [], "chrom_len": 0, "contigs": [], "gaps": 0, "gap_bp": 0}
    records = {}
    for header, sequence in read_fasta(GENOMES / genome_name, full_header=True):
        seqid = header.split()[0]
        number = chromosome_number(header)
        sequence = sequence.upper()
        if number is not None and 1 <= number <= 14:
            target = CHROMOSOME_MAPS[key][number]
            reverse = number in REVERSE_COMPLEMENT[key]
            if reverse:
                sequence = sequence.translate(COMPLEMENT)[::-1]
            name = f"chr{target}"
            stats["chrom_len"] += len(sequence)
        else:
            name, reverse = seqid, False
        seq_map[seqid] = (name, reverse, len(sequence))
        records[name] = sequence

    chroms = sorted((n for n in records if n.startswith("chr")), key=lambda n: int(n[3:]))
    order = chroms + sorted(n for n in records if not n.startswith("chr"))
    if len(chroms) != 14:
        raise SystemExit(f"{code}: expected 14 chromosomes, found {len(chroms)}")

    with open(fasta, "w") as fa, open(INPUTS / f"{code}.gaps.bed", "w") as gaps, \
            open(INPUTS / f"{code}.contigs.bed", "w") as contigs, \
            open(INPUTS / f"{code}.chrom.sizes", "w") as sizes:
        for name in order:
            sequence = records[name]
            write_fasta(fa, name, sequence)
            sizes.write(f"{name}\t{len(sequence)}\n")
            stats["lengths"].append(len(sequence))
            previous = 0
            for match in GAP.finditer(sequence):
                gaps.write(f"{name}\t{match.start()}\t{match.end()}\n")
                stats["gaps"] += 1
                stats["gap_bp"] += match.end() - match.start()
                if match.start() > previous:
                    contigs.write(f"{name}\t{previous}\t{match.start()}\n")
                    stats["contigs"].append(match.start() - previous)
                previous = match.end()
            if len(sequence) > previous:
                contigs.write(f"{name}\t{previous}\t{len(sequence)}\n")
                stats["contigs"].append(len(sequence) - previous)
    return seq_map, stats


def attributes(text):
    out = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = unquote(v)
    return out


def lift(seq_map, seqid, start, end, strand):
    name, reverse, length = seq_map[seqid]
    if reverse:
        start, end = length - end + 1, length - start + 1
        strand = {"+": "-", "-": "+"}.get(strand, strand)
    return name, start, end, strand


def prepare_annotation(code, gff_name, protein_name, style, seq_map):
    """Write longest-isoform proteins and a normalized gene table."""
    proteins = {name: seq.rstrip("*.") for name, seq in read_fasta(RESOURCES / protein_name)}
    products = {}
    if style == "refseq":
        for header, _ in read_fasta(RESOURCES / protein_name, full_header=True):
            pid, _, rest = header.partition(" ")
            products[pid] = re.sub(r"\s*\[Acropora millepora\]$", "", rest)

    cds = {}  # protein -> [gene, seqid, start, end, strand]
    with open_text(RESOURCES / gff_name) as gff:
        for line in gff:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS" or fields[0] not in seq_map:
                continue
            attrs = attributes(fields[8])
            if style == "augustus":
                transcript = attrs["Parent"]
                protein = f"{fields[0]}.{transcript}"
                gene = f"{fields[0]}.{transcript.rsplit('.', 1)[0]}"
            else:
                protein = attrs.get("protein_id", attrs.get("Name", ""))
                gene = attrs.get("gene", protein)
            if protein not in proteins:
                continue
            start, end = int(fields[3]), int(fields[4])
            if protein in cds:
                cds[protein][2] = min(cds[protein][2], start)
                cds[protein][3] = max(cds[protein][3], end)
            else:
                cds[protein] = [gene, fields[0], start, end, fields[6]]

    best = {}
    for protein, (gene, *_rest) in cds.items():
        key = (len(proteins[protein]), protein)
        if gene not in best or key > best[gene][0]:
            best[gene] = (key, protein)

    rows = []
    for gene, (_, protein) in best.items():
        _, seqid, start, end, strand = cds[protein]
        chrom, start, end, strand = lift(seq_map, seqid, start, end, strand)
        gene_id = f"{code}_{gene}"
        rows.append((chrom, start, end, strand, gene_id, protein, products.get(protein, "")))
    rows.sort(key=lambda r: (0 if r[0].startswith("chr") else 1,
                             int(r[0][3:]) if r[0].startswith("chr") else r[0], r[1]))

    ann_dir = ensure(INPUTS / "annotations")
    with open(ann_dir / f"{code}.genes.tsv", "w") as table, \
            open(ann_dir / f"{code}.proteins.faa", "w") as fasta:
        table.write("gene_id\tchrom\tstart\tend\tstrand\tprotein_id\tprotein_length\tproduct\n")
        for chrom, start, end, strand, gene_id, protein, product in rows:
            table.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\t{protein}\t"
                        f"{len(proteins[protein])}\t{product}\n")
            write_fasta(fasta, gene_id, proteins[protein])
    on_chrom = sum(r[0].startswith("chr") for r in rows)
    return len(rows), on_chrom, len(cds)


def main():
    ensure(INPUTS)
    with open(INPUTS / "assembly_stats.tsv", "w") as out:
        out.write("species\ttotal_bp\tsequences\tchromosome_bp\tpct_in_chromosomes\t"
                  "unplaced_bp\tscaffold_N50\tcontigs\tcontig_N50\tgaps\tgap_bp\t"
                  "annotated_genes\tgenes_on_chromosomes\tannotated_isoforms\n")
        for code, _name, key, genome, gff, prot, style in SPECIES:
            print(f"[01] {code}: genome", flush=True)
            seq_map, s = prepare_genome(code, key, genome)
            genes = (0, 0, 0)
            if code in ANNOTATED:
                print(f"[01] {code}: annotation", flush=True)
                genes = prepare_annotation(code, gff, prot, style, seq_map)
            total = sum(s["lengths"])
            out.write(f"{code}\t{total}\t{len(s['lengths'])}\t{s['chrom_len']}\t"
                      f"{100 * s['chrom_len'] / total:.2f}\t{total - s['chrom_len']}\t"
                      f"{n50(s['lengths'])}\t{len(s['contigs'])}\t{n50(s['contigs'])}\t"
                      f"{s['gaps']}\t{s['gap_bp']}\t{genes[0]}\t{genes[1]}\t{genes[2]}\n")


if __name__ == "__main__":
    main()
