#!/usr/bin/env python3
"""Compare three Acropora genomes with DIAMOND and MCScanX.

Run normally:
    ./scripts/compare_genomes_mcscanx.py --threads 14

Rebuild existing results:
    ./scripts/compare_genomes_mcscanx.py --threads 14 --force

The code is arranged in the same order as the analysis:
    Phase 1 - prepare proteins and gene positions
    Phase 2 - build a DIAMOND database
    Phase 3 - find similar proteins
    Phase 4 - find syntenic blocks with MCScanX
    Phase 5 - classify duplicate genes in each genome
    Phase 6 - make dual, circle, dot, bar, and duplication plots
    Phase 7 - write a Markdown report
"""

import argparse
import gzip
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote


# All paths are fixed relative to this script, so it can be run from anywhere.
PROJECT = Path(__file__).resolve().parent.parent
RESOURCES = PROJECT / "resources"
GENOMES = PROJECT / "genomes"
RESULTS = PROJECT / "results" / "mcscanx"
MCSCANX_DIR = PROJECT / "programs" / "MCScanX"
PLOTTER_DIR = MCSCANX_DIR / "downstream_analyses"

# name, two-letter MCScanX prefix, genome, GFF, proteins, annotation format
SPECIES = [
    ("A_digitifera", "Ad", "AdiV3.2Ref_genome.fa.gz",
     "AdiV3.2Ref.gff.gz", "AdiV3.2Ref_aa.fa.gz", "augustus"),
    ("A_millepora", "Am", "A_millepora_v2.1.fasta",
     "A_millepora_v2.1.gff", "A_millepora_v2.1_protein.faa", "refseq"),
    ("A_tenuis", "At", "AteV2.2Ref_genome.fa.gz",
     "AteV2.2Ref.gff.gz", "AteV2.2Ref_aa.fa.gz", "augustus"),
]


def options():
    """Read the few command-line options supported by this workflow."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=min(14, os.cpu_count() or 1))
    parser.add_argument("--prepare-only", action="store_true",
                        help="stop after creating the MCScanX input files")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="repeat DIAMOND and MCScanX even if results exist")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    return args


def open_text(path):
    """Open a normal or gzip-compressed text file."""
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_fasta(path):
    """Read a FASTA file into {sequence_name: sequence}."""
    records, name, pieces = {}, None, []
    with open_text(path) as fasta:
        for line in fasta:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(pieces).rstrip("*")
                name, pieces = line[1:].split()[0], []
            elif name is not None:
                pieces.append(line.strip())
    if name is not None:
        records[name] = "".join(pieces).rstrip("*")
    return records


def gff_attributes(text):
    """Turn GFF column 9 (ID=x;Parent=y) into a dictionary."""
    answer = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            answer[key] = unquote(value)
    return answer


def chromosome(style, sequence_name):
    """Return chromosome number 1--14, or None for a scaffold."""
    if style == "augustus":
        match = re.fullmatch(r"chr([1-9]|1[0-4])Ref", sequence_name)
        return int(match.group(1)) if match else None

    # Millepora chromosome accessions run from NC_058066.1 to NC_058079.1.
    match = re.fullmatch(r"NC_(\d+)\.1", sequence_name)
    number = int(match.group(1)) - 58065 if match else 0
    return number if 1 <= number <= 14 else None


def read_gene_positions(prefix, style, gff_path, proteins):
    """Match protein IDs to their genomic coordinates.

    AUGUSTUS has one useful 'transcript' row per protein. RefSeq splits each
    protein over several CDS rows, so those rows must first be merged.
    """
    protein_positions = {}

    with open_text(gff_path) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            seqid, feature = fields[0], fields[2]
            number = chromosome(style, seqid)
            attrs = gff_attributes(fields[8])
            if number is None:
                continue

            if style == "augustus" and feature == "transcript":
                transcript = attrs.get("ID", "")
                protein = f"{seqid}.{transcript}"
                if protein in proteins:
                    protein_positions[protein] = {
                        "gene": f"{seqid}:{attrs.get('Parent', transcript)}",
                        "chrom": f"{prefix}{number}",
                        "start": int(fields[3]), "end": int(fields[4]),
                    }

            elif style == "refseq" and feature == "CDS":
                protein = attrs.get("protein_id", "")
                if protein not in proteins:
                    continue
                start, end = int(fields[3]), int(fields[4])
                gene = f"{seqid}:{attrs.get('gene', attrs.get('Parent', protein))}"
                if protein not in protein_positions:
                    protein_positions[protein] = {
                        "gene": gene, "chrom": f"{prefix}{number}",
                        "start": start, "end": end,
                    }
                else:
                    # Expand the interval as additional CDS exons are read.
                    old = protein_positions[protein]
                    old["start"] = min(old["start"], start)
                    old["end"] = max(old["end"], end)

    # Several proteins may be isoforms of one gene. Keep the longest isoform.
    longest = {}
    for protein, position in protein_positions.items():
        gene = position["gene"]
        old_protein = longest.get(gene)
        # Protein ID breaks exact length ties, making repeated runs identical.
        if old_protein is None or (len(proteins[protein]), protein) > (
                len(proteins[old_protein]), old_protein):
            longest[gene] = protein

    selected = [(protein, protein_positions[protein]) for protein in longest.values()]
    selected.sort(key=lambda item: (
        int(item[1]["chrom"][2:]), item[1]["start"], item[1]["end"]
    ))
    return selected


# PHASE 1 --------------------------------------------------------------------
def prepare_inputs():
    """Create the combined protein FASTA, MCScanX GFF, and ID mapping table."""
    proteins_out = RESULTS / "acropora.proteins.fasta"
    gff_out = RESULTS / "acropora.gff"
    map_out = RESULTS / "gene_id_map.tsv"
    total = 0

    with proteins_out.open("w") as fasta, gff_out.open("w") as gff, map_out.open("w") as table:
        table.write("species\toriginal_id\tmcscanx_id\tchromosome\tstart\tend\n")

        for name, prefix, genome, gff_name, protein_name, style in SPECIES:
            proteins = read_fasta(RESOURCES / protein_name)
            genes = read_gene_positions(prefix, style, RESOURCES / gff_name, proteins)
            print(f"  {name}: selected {len(genes):,} genes", flush=True)

            for old_id, position in genes:
                # Species prefixes prevent identical IDs in different genomes.
                new_id = f"{prefix}_{old_id}"
                sequence = proteins[old_id]
                fasta.write(f">{new_id}\n")
                for start in range(0, len(sequence), 80):
                    fasta.write(sequence[start:start + 80] + "\n")
                gff.write(
                    f"{position['chrom']}\t{new_id}\t"
                    f"{position['start']}\t{position['end']}\n"
                )
                table.write(
                    f"{name}\t{old_id}\t{new_id}\t{position['chrom']}\t"
                    f"{position['start']}\t{position['end']}\n"
                )
            total += len(genes)

    print(f"  Total: {total:,} genes", flush=True)
    return proteins_out, gff_out


def run(command, log_name):
    """Show an exact command, run it, and save its messages in a log."""
    command = [str(part) for part in command]
    print("  $ " + shlex.join(command), flush=True)
    with (RESULTS / log_name).open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


def summarize(collinearity):
    """Write a short table of block and collinear-gene-pair counts."""
    names = {prefix: name for name, prefix, *_ in SPECIES}
    counts = {}
    pattern = re.compile(r"^## Alignment .*?N=(\d+)\s+(\S+)&(\S+)\s")
    with collinearity.open() as source:
        for line in source:
            match = pattern.search(line)
            if match:
                pair = tuple(sorted((names[match.group(2)[:2]], names[match.group(3)[:2]])))
                blocks, genes = counts.get(pair, (0, 0))
                counts[pair] = blocks + 1, genes + int(match.group(1))
    with (RESULTS / "comparison_summary.tsv").open("w") as table:
        table.write("species_1\tspecies_2\tsyntenic_blocks\tcollinear_gene_pairs\n")
        for pair, (blocks, genes) in sorted(counts.items()):
            table.write(f"{pair[0]}\t{pair[1]}\t{blocks}\t{genes}\n")


# PHASE 5 --------------------------------------------------------------------
def classify_duplicates(combined_fasta, combined_gff, threads, force):
    """Run Duplicate_gene_classifier separately for each species."""
    duplicate_dir = RESULTS / "duplicate_genes"
    duplicate_dir.mkdir(exist_ok=True)
    all_proteins = read_fasta(combined_fasta)
    classifier = MCSCANX_DIR / "duplicate_gene_classifier"

    # MCScanX requires a one-genome GFF and self-vs-self BLAST for this tool.
    if not os.access(classifier, os.X_OK):
        run(["make", "-C", MCSCANX_DIR, "mcscanx"], "duplicate_genes/build.log")

    for name, prefix, *_ in SPECIES:
        analysis = duplicate_dir / name
        protein_file = duplicate_dir / f"{name}.proteins.fasta"
        gene_type = duplicate_dir / f"{name}.gene_type"
        if gene_type.exists() and gene_type.stat().st_size and not force:
            print(f"  {name}: reuse {gene_type.name}")
            continue

        selected = {key: value for key, value in all_proteins.items()
                    if key.startswith(prefix + "_")}
        with protein_file.open("w") as fasta:
            for protein, sequence in selected.items():
                fasta.write(f">{protein}\n")
                for start in range(0, len(sequence), 80):
                    fasta.write(sequence[start:start + 80] + "\n")
        with combined_gff.open() as source, Path(str(analysis) + ".gff").open("w") as target:
            for line in source:
                if line.split("\t", 2)[1].startswith(prefix + "_"):
                    target.write(line)

        # Build and search a separate protein database for this species.
        # diamond makedb --in SPECIES.proteins.fasta --db SPECIES
        run(["diamond", "makedb", "--in", protein_file, "--db", analysis],
            f"duplicate_genes/{name}.makedb.log")
        # diamond blastp --query PROTEINS --db SPECIES --out SPECIES.blast ...
        run([
            "diamond", "blastp", "--query", protein_file, "--db", analysis,
            "--out", str(analysis) + ".blast", "--outfmt", "6",
            "--evalue", "1e-10", "--max-target-seqs", "5",
            "--threads", threads, "--sensitive",
        ], f"duplicate_genes/{name}.blastp.log")
        # programs/MCScanX/duplicate_gene_classifier SPECIES
        run([classifier, analysis], f"duplicate_genes/{name}.classifier.log")

    return summarize_duplicates(duplicate_dir)


def summarize_duplicates(duplicate_dir):
    """Count duplicate classes and draw a proportional stacked-bar SVG."""
    categories = ["Singleton", "Dispersed", "Proximal", "Tandem", "WGD/segmental"]
    colors = ["#9ca3af", "#4e79a7", "#f28e2b", "#e15759", "#59a14f"]
    totals = {}
    with (duplicate_dir / "duplicate_gene_summary.tsv").open("w") as table:
        table.write("species\tcategory\tcode\tgenes\tpercent\n")
        for name, *_ in SPECIES:
            counts = [0] * 5
            with (duplicate_dir / f"{name}.gene_type").open() as source:
                for line in source:
                    fields = line.split()
                    if len(fields) == 2:
                        counts[int(fields[1])] += 1
            totals[name] = counts
            gene_total = sum(counts)
            for code, category in enumerate(categories):
                percent = 100 * counts[code] / gene_total
                table.write(f"{name}\t{category}\t{code}\t{counts[code]}\t{percent:.2f}\n")

    # Write an SVG directly, avoiding an additional plotting dependency.
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="520">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="550" y="38" text-anchor="middle" font-family="sans-serif" '
        'font-size="24" font-weight="bold">MCScanX duplicate-gene classes</text>',
    ]
    for row, (name, counts) in enumerate(totals.items()):
        y, x, gene_total = 90 + row * 105, 190, sum(counts)
        svg.append(f'<text x="175" y="{y + 35}" text-anchor="end" '
                   f'font-family="sans-serif" font-size="18">{name}</text>')
        for color, count in zip(colors, counts):
            segment = 820 * count / gene_total
            svg.append(f'<rect x="{x:.2f}" y="{y}" width="{segment:.2f}" '
                       f'height="55" fill="{color}"/>')
            if segment > 65:
                svg.append(f'<text x="{x + segment / 2:.2f}" y="{y + 34}" '
                           f'text-anchor="middle" font-family="sans-serif" '
                           f'font-size="14" fill="white">{100 * count / gene_total:.1f}%</text>')
            x += segment
    for index, (category, color) in enumerate(zip(categories, colors)):
        x = 110 + index * 195
        svg.append(f'<rect x="{x}" y="445" width="20" height="20" fill="{color}"/>')
        svg.append(f'<text x="{x + 27}" y="461" font-family="sans-serif" '
                   f'font-size="14">{category}</text>')
    svg.append("</svg>\n")
    (RESULTS / "plots").mkdir(exist_ok=True)
    svg_file = RESULTS / "plots" / "E_duplicate_gene_classes.svg"
    svg_file.write_text("\n".join(svg))
    # Also provide a PNG when ImageMagick is available.
    if shutil.which("convert"):
        run(["convert", svg_file, RESULTS / "plots" / "E_duplicate_gene_classes.png"],
            "plots/E_duplicate_gene_classes.log")
    return totals


# PHASE 6 --------------------------------------------------------------------
def make_plots(gff, collinearity):
    """Make dual-synteny, circle, dot, and bar plots with MCScanX tools."""
    plot_dir = RESULTS / "plots"
    plot_dir.mkdir(exist_ok=True)
    classpath = str(PLOTTER_DIR)

    # Compile MCScanX's four Java plotters only when necessary.
    plotters = ("dual_synteny_plotter", "circle_plotter", "dot_plotter", "bar_plotter")
    if not all((PLOTTER_DIR / f"{name}.class").exists() for name in plotters):
        run(["make", "-C", PLOTTER_DIR], "plotter_build.log")

    chromosome_lists = {
        prefix: ",".join(f"{prefix}{number}" for number in range(1, 15))
        for _, prefix, *_ in SPECIES
    }

    # A, C, and D are pairwise plots. Their control-file layouts are identical.
    for left_index in range(len(SPECIES) - 1):
        for right_index in range(left_index + 1, len(SPECIES)):
            left_name, left_prefix, *_ = SPECIES[left_index]
            right_name, right_prefix, *_ = SPECIES[right_index]
            stem = f"{left_name}_vs_{right_name}"
            for letter, plotter, x_size, y_size in (
                ("A", "dual_synteny_plotter", 900, 1600),
                ("C", "dot_plotter", 1600, 1600),
                ("D", "bar_plotter", 1600, 1600),
            ):
                control = plot_dir / f"{letter}_{stem}.{plotter}.ctl"
                image = plot_dir / f"{letter}_{stem}.{plotter}.png"
                control.write_text(
                    f"{x_size}\n{y_size}\n{chromosome_lists[left_prefix]}\n"
                    f"{chromosome_lists[right_prefix]}\n"
                )
                # java -cp PLOTTER_DIR PLOTTER -g GFF -s COLLINEARITY
                #   -c CONTROL -o PNG
                run([
                    "java", "-Djava.awt.headless=true", "-cp", classpath,
                    plotter, "-g", gff, "-s", collinearity,
                    "-c", control, "-o", image,
                ], f"plots/{letter}_{stem}.{plotter}.log")

    # B is one circular overview containing all 42 chromosomes.
    # java -cp PLOTTER_DIR circle_plotter -g GFF -s COLLINEARITY -c CTL -o PNG
    circle_control = plot_dir / "B_all_species.circle_plotter.ctl"
    all_chromosomes = ",".join(chromosome_lists[prefix] for _, prefix, *_ in SPECIES)
    circle_control.write_text("2000\n" + all_chromosomes + "\n")
    run([
        "java", "-Djava.awt.headless=true", "-cp", classpath,
        "circle_plotter", "-g", gff, "-s", collinearity,
        "-c", circle_control, "-o", plot_dir / "B_all_species.circle_plotter.png",
    ], "plots/B_all_species.circle_plotter.log")


# PHASE 7 --------------------------------------------------------------------
def write_report(duplicate_counts):
    """Create a Markdown report whose tables and insights come from results."""
    synteny_rows = []
    with (RESULTS / "comparison_summary.tsv").open() as table:
        next(table)
        for line in table:
            left, right, blocks, genes = line.rstrip().split("\t")
            synteny_rows.append((left, right, int(blocks), int(genes)))

    strongest = max(synteny_rows, key=lambda row: row[3])
    categories = ["Singleton", "Dispersed", "Proximal", "Tandem", "WGD/segmental"]
    lines = [
        "# Acropora MCScanX analysis report", "",
        "## Analysis overview", "",
        "Three chromosome-level annotated genomes were compared: *A. digitifera*, "
        "*A. millepora*, and *A. tenuis*. The longest protein isoform per gene was "
        "used. DIAMOND supplied protein similarities and MCScanX identified "
        "inter-species collinear blocks.", "",
        "## Main insights", "",
        f"- **{strongest[0]} vs {strongest[1]}** had the largest number of "
        f"collinear gene pairs ({strongest[3]:,} in {strongest[2]:,} blocks).",
    ]
    for species, counts in duplicate_counts.items():
        total = sum(counts)
        duplicate_only = max(range(1, 5), key=lambda code: counts[code])
        lines.append(
            f"- In **{species}**, the largest duplicate category was "
            f"**{categories[duplicate_only]}** ({counts[duplicate_only]:,} genes; "
            f"{100 * counts[duplicate_only] / total:.1f}%)."
        )
    longest = max(synteny_rows, key=lambda row: row[3] / row[2])
    lines.append(
        f"- **{longest[0]} vs {longest[1]}** had the longest blocks on average "
        f"({longest[3] / longest[2]:.1f} collinear gene pairs per block), suggesting "
        "less fragmented chromosome-scale collinearity in this comparison."
    )
    most_wgd = max(duplicate_counts, key=lambda name: duplicate_counts[name][4] /
                   sum(duplicate_counts[name]))
    wgd_counts = duplicate_counts[most_wgd]
    lines.append(
        f"- **{most_wgd}** had the highest WGD/segmental fraction "
        f"({100 * wgd_counts[4] / sum(wgd_counts):.1f}%; {wgd_counts[4]:,} genes). "
        "This is a computational classification and is not by itself evidence of a "
        "lineage-specific whole-genome duplication."
    )

    lines += [
        "", "## Synteny summary", "",
        "| Species 1 | Species 2 | Blocks | Collinear gene pairs |",
        "|---|---|---:|---:|",
    ]
    for left, right, blocks, genes in synteny_rows:
        lines.append(f"| {left} | {right} | {blocks:,} | {genes:,} |")

    lines += [
        "", "## Duplicate-gene classification", "",
        "| Species | Singleton | Dispersed | Proximal | Tandem | WGD/segmental |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for species, counts in duplicate_counts.items():
        lines.append(f"| {species} | " + " | ".join(f"{value:,}" for value in counts) + " |")

    lines += ["", "## Figures", "",
              "| Comparison | A: dual synteny | C: dot plot | D: bar plot |",
              "|---|---|---|---|"]
    for left_index in range(len(SPECIES) - 1):
        for right_index in range(left_index + 1, len(SPECIES)):
            left, right = SPECIES[left_index][0], SPECIES[right_index][0]
            stem = f"{left}_vs_{right}"
            lines.append(
                f"| {left} vs {right} | [dual](plots/A_{stem}.dual_synteny_plotter.png) "
                f"| [dot](plots/C_{stem}.dot_plotter.png) "
                f"| [bar](plots/D_{stem}.bar_plotter.png) |"
            )
    lines += [
        "", "- **B:** [Three-species circle plot](plots/B_all_species.circle_plotter.png).",
        "- **E:** [Duplicate-gene classes](plots/E_duplicate_gene_classes.png) "
        "([vector SVG](plots/E_duplicate_gene_classes.svg)).", "",
        "## Interpretation notes", "",
        "MCScanX duplicate codes are: 0 singleton, 1 dispersed, 2 proximal, "
        "3 tandem, and 4 WGD/segmental. Classification was performed separately "
        "for each genome using a species-specific self-vs-self DIAMOND search. "
        "The results describe computational homology and collinearity calls; they "
        "should be interpreted alongside annotation quality and phylogenetic evidence.", "",
    ]
    (RESULTS / "REPORT.md").write_text("\n".join(lines))


def check_files():
    """Give a clear error before starting if an input or program is missing."""
    files = [MCSCANX_DIR]
    for _, _, genome, gff, proteins, _ in SPECIES:
        files.extend((GENOMES / genome, RESOURCES / gff, RESOURCES / proteins))
    missing = [str(path) for path in files if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing:\n" + "\n".join(missing))


def main():
    args = options()
    RESULTS.mkdir(parents=True, exist_ok=True)
    check_files()

    proteins = RESULTS / "acropora.proteins.fasta"
    gff = RESULTS / "acropora.gff"
    blast = RESULTS / "acropora.blast"
    collinearity = RESULTS / "acropora.collinearity"

    # Reuse completed phases unless --force was supplied.
    if (collinearity.exists() and collinearity.stat().st_size
            and not args.force and not args.prepare_only):
        print("PHASES 1-4: reuse existing MCScanX results")
        summarize(collinearity)
    else:
        print("PHASE 1: prepare proteins and gene positions")
        proteins, gff = prepare_inputs()
        if args.prepare_only:
            return

        print("\nPHASE 2: build the DIAMOND database")
        # diamond makedb --in PROTEINS --db DATABASE
        run(["diamond", "makedb", "--in", proteins, "--db", RESULTS / "acropora"],
            "diamond_makedb.log")

        print("\nPHASE 3: run the all-against-all protein search")
        # diamond blastp --query PROTEINS --db DATABASE --out BLAST --outfmt 6
        #   --evalue 1e-10 --max-target-seqs 5 --threads 14 --sensitive
        run([
            "diamond", "blastp", "--query", proteins,
            "--db", RESULTS / "acropora", "--out", blast,
            "--outfmt", "6", "--evalue", "1e-10",
            "--max-target-seqs", "5", "--threads", args.threads, "--sensitive",
        ], "diamond_blastp.log")

        print("\nPHASE 4: find inter-species syntenic blocks")
        executable = MCSCANX_DIR / "MCScanX"
        if not os.access(executable, os.X_OK):
            # make -C programs/MCScanX mcscanx
            run(["make", "-C", MCSCANX_DIR, "mcscanx"], "mcscanx_build.log")
        # programs/MCScanX/MCScanX results/mcscanx/acropora -b 2
        run([executable, RESULTS / "acropora", "-b", "2"], "mcscanx.log")
        summarize(collinearity)

    print("\nPHASE 5: classify duplicate genes within each species")
    duplicate_counts = classify_duplicates(proteins, gff, args.threads, args.force)

    if not args.skip_plots:
        print("\nPHASE 6: make dual, circle, dot, bar, and duplication plots")
        make_plots(gff, collinearity)

    print("\nPHASE 7: write the Markdown report")
    write_report(duplicate_counts)
    print(f"\nFinished. Results: {RESULTS}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
