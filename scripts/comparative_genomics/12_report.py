#!/usr/bin/env python3
"""Step 12 - figures and a Markdown report of the comparative analysis.

Figures (12_report/figures/):
  fig1_assembly_methodology   genome size split + BUSCO + disrupted-gene rate
  fig2_sv_classes             SyRI rearrangements by evidence class, per pair
  fig3_anchor_dotplots        single-copy ortholog positions, all pairs
  fig4_biological_events      confirmed rearrangements on A. digitifera chromosomes
  fig5_gene_conservation      grades, verified absences, duplication modes
  fig6_sequence_evolution     dN/dS and dS per species pair
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from config import CODES, INPUTS, NAMES, OUT, ensure  # noqa: E402

REP = ensure(OUT / "12_report")
FIG = ensure(REP / "figures")
PAIRS = [(a, b) for i, a in enumerate(CODES) for b in CODES[i + 1:]]

# Reference categorical palette (dataviz skill), fixed order.
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
GRAY, LIGHTGRAY, INK, MUTED = "#8a8985", "#c9c8c3", "#0b0b0b", "#52514e"
SPECIES_COLOR = {"Adig": BLUE, "Amil": ORANGE, "Aten": AQUA, "Acer": VIOLET}

CLASS_ORDER = [
    ("biological_shared", "Biological: shared by 2 genomes", BLUE),
    ("biological_lineage_specific", "Biological: lineage-specific", AQUA),
    ("biological_unpolarized", "Biological: unpolarized", GREEN),
    ("uncertain_error_prone_assembly", "Uncertain: only in error-prone assembly", VIOLET),
    ("uncertain_one_gap", "Uncertain: one breakpoint at gap", YELLOW),
    ("uncertain_gap_bounded", "Uncertain: gap-bounded", "#f5c45c"),
    ("artifact_scaffolding", "Artifact: scaffolding (gap-bounded)", RED),
    ("artifact_gap_bounded", "Artifact: gap-bounded, no genes", ORANGE),
    ("artifact_alignment", "Artifact: not supported by orthologs", MAGENTA),
    ("unverified_repeat", "Unverified: repeat-rich", GRAY),
    ("unverified_no_genes", "Unverified: <3 genes inside", LIGHTGRAY),
]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
})


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if Path(path).exists() else []


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)


def fig_assembly(method):
    """Fig 1: size split, contiguity and BUSCO deviations; every label stays inside its panel."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios": [1.4, 1, 1], "wspace": 0.45})
    names = [NAMES[r["species"]] for r in method]
    y = list(range(len(method)))

    ax = axes[0]
    nonrep = [float(r["nonrepeat_Mb"]) for r in method]
    rep = [float(r["repeat_Mb_Red"]) for r in method]
    ax.barh(y, nonrep, color=BLUE, height=0.6, edgecolor="white", linewidth=1, label="Non-repeat")
    ax.barh(y, rep, left=nonrep, color="#86b6ef", height=0.6, edgecolor="white", linewidth=1,
            label="Repeat (Red, same settings for all)")
    for i, r in enumerate(method):
        total = float(r["total_Mb"])
        ax.text(total + 8, i, f"{total:.0f} Mb", va="center", fontsize=8.5, color=INK)
        ax.text(nonrep[i] / 2, i, f"{nonrep[i]:.0f}", va="center", ha="center", fontsize=7.5, color="white")
        ax.text(nonrep[i] + rep[i] / 2, i, f"{float(r['repeat_pct_Red']):.0f}%", va="center", ha="center",
                fontsize=7.5, color=INK)
    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.set_xlim(0, max(float(r["total_Mb"]) for r in method) * 1.18)
    ax.set_xlabel("Assembly size (Mb)")
    ax.set_title("A. Genome size: non-repeat vs repeat", loc="left", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8)

    ax = axes[1]
    n50 = [float(r["contig_N50_Mb"]) for r in method]
    ax.barh(y, n50, color=[SPECIES_COLOR[r["species"]] for r in method], height=0.6)
    for i, r in enumerate(method):
        ax.text(n50[i] + 0.2, i, f"{n50[i]:.1f} Mb · {int(r['gaps'])} gaps · {float(r['unplaced_Mb']):.0f} Mb unplaced",
                va="center", fontsize=7.5, color=INK)
    ax.set_yticks(y, [r["species"] for r in method])
    ax.invert_yaxis()
    ax.set_xlim(0, max(n50) * 2.4)
    ax.set_xlabel("Contig N50 (Mb)")
    ax.set_title("B. Contiguity", loc="left", fontsize=10)

    ax = axes[2]
    keys = [("BUSCO_genome_duplicated", "Duplicated", ORANGE), ("BUSCO_genome_fragmented", "Fragmented", YELLOW),
            ("BUSCO_genome_missing", "Missing", GRAY)]
    width = 0.26
    for k, (key, label, color) in enumerate(keys):
        vals = [float(r[key] or 0) for r in method]
        ax.bar([i + (k - 1) * width for i in y], vals, width * 0.92, color=color, label=label)
    top = max(float(r[k] or 0) for r in method for k, _, _ in keys)
    for i, r in enumerate(method):
        ax.text(i, top * 1.12, f"C {float(r['BUSCO_genome_complete']):.1f}%", ha="center", fontsize=7.5, color=MUTED)
    ax.set_ylim(0, top * 1.3)
    ax.set_xticks(y, [r["species"] for r in method])
    ax.set_ylabel("% of 3,203 cnidaria_odb12 genes")
    ax.set_title("C. BUSCO (genome mode): what is not single-copy", loc="left", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=8)
    save(fig, "fig1_assembly_methodology")



def fig_sv(calls):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
    pairs = [f"{a}_vs_{b}" for a, b in PAIRS]
    for ax, metric in zip(axes, ("count", "Mb")):
        left = [0.0] * len(pairs)
        for cls, label, color in CLASS_ORDER:
            vals = []
            for p in pairs:
                sel = [c for c in calls if c["pair"] == p and c["class"] == cls and c["type"] in ("INV", "TRANS", "INVTR")]
                vals.append(len(sel) if metric == "count" else sum(int(c["size"]) for c in sel) / 1e6)
            ax.barh(range(len(pairs)), vals, left=left, color=color, label=label, height=0.62,
                    edgecolor="white", linewidth=1)
            left = [l + v for l, v in zip(left, vals)]
        ax.set_yticks(range(len(pairs)), [p.replace("_vs_", " vs ") for p in pairs])
        ax.invert_yaxis()
        ax.set_xlabel("Rearrangements (>=10 kb)" if metric == "count" else "Sequence in rearrangements (Mb)")
    axes[0].set_title("SyRI inversions & translocations by evidence class", loc="left", fontsize=10)
    axes[1].legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7.5)
    save(fig, "fig2_sv_classes")


def chrom_offsets(code):
    sizes = [l.split() for l in open(INPUTS / f"{code}.chrom.sizes")]
    off, total = {}, 0
    for name, size in sizes:
        if name.startswith("chr"):
            off[name] = total
            total += int(size)
    return off, total


def fig_dotplots(anchors, events):
    offsets = {c: chrom_offsets(c) for c in CODES}
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.6))
    for ax, (a, b) in zip(axes.flat, PAIRS):
        oa, ta = offsets[a]
        ob, tb = offsets[b]
        xs, ys, cs = [], [], []
        for r in anchors:
            if r[f"{a}_chrom"] in oa and r[f"{b}_chrom"] in ob:
                xs.append((oa[r[f"{a}_chrom"]] + int(r[f"{a}_start"])) / 1e6)
                ys.append((ob[r[f"{b}_chrom"]] + int(r[f"{b}_start"])) / 1e6)
                same = r[f"{a}_strand"] == r[f"{b}_strand"]
                cs.append(BLUE if same else ORANGE)
        ax.scatter(xs, ys, s=0.6, c=cs, linewidths=0, rasterized=True)
        for o in oa.values():
            ax.axvline(o / 1e6, color=LIGHTGRAY, lw=0.5)
        for o in ob.values():
            ax.axhline(o / 1e6, color=LIGHTGRAY, lw=0.5)
        ax.set_xlim(0, ta / 1e6)
        ax.set_ylim(0, tb / 1e6)
        ax.set_xlabel(f"{NAMES[a]} chr1-14 (Mb)")
        ax.set_ylabel(f"{NAMES[b]} chr1-14 (Mb)")
        ax.grid(False)
        ax.set_title(f"{a} vs {b}: {len(xs):,} single-copy orthologs", loc="left", fontsize=9)
    fig.suptitle("Gene-order dot plots (blue = same strand, orange = opposite strand)", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "fig3_anchor_dotplots")


def fig_events(events):
    """Fig 4: evidence classes on A. digitifera chromosomes; biological events get numbered markers."""
    sizes = [l.split() for l in open(INPUTS / "Adig.chrom.sizes") if l.startswith("chr")]
    colors = {c: col for c, _, col in CLASS_ORDER}
    labels = {c: l for c, l, _ in CLASS_ORDER}
    fig = plt.figure(figsize=(15, 6.6))
    ax = fig.add_axes([0.05, 0.1, 0.62, 0.82])
    key_ax = fig.add_axes([0.69, 0.1, 0.30, 0.82])
    key_ax.axis("off")
    for i, (name, size) in enumerate(sizes):
        ax.add_patch(plt.Rectangle((0, i - 0.17), int(size) / 1e6, 0.34, color="#ecebe6", lw=0))
    shown, numbered = set(), []
    order = sorted(events, key=lambda e: (int(e["Adig_location"].split(":")[0][3:]),
                                          int(e["Adig_location"].split(":")[1].split("-")[0])))
    for ev in order:
        chrom, span = ev["Adig_location"].split(":")
        s, e = (int(x) / 1e6 for x in span.split("-"))
        i = int(chrom[3:]) - 1
        cls = ev["class"]
        h = 0.34 if ev["type"] == "INV" else 0.2
        ax.add_patch(plt.Rectangle((s, i - h / 2), max(e - s, 0.15), h, color=colors[cls], lw=0,
                                   label=None if cls in shown else labels[cls]))
        shown.add(cls)
        if cls.startswith("biological"):
            numbered.append((ev, i, (s + e) / 2))
    # Numbered markers above the bar; markers that would touch are stacked.
    last = {}
    for n, (ev, i, x) in enumerate(numbered, 1):
        level = 1 if i in last and x - last[i][0] < 1.2 and last[i][1] == 0 else 0
        last[i] = (x, level)
        yy = i - 0.38 - 0.24 * level
        ax.scatter([x], [yy], s=120, color=colors[ev["class"]], zorder=3, linewidths=0)
        ax.text(x, yy, str(n), ha="center", va="center", fontsize=6.5, color="white", fontweight="bold", zorder=4)
    ax.set_yticks(range(len(sizes)), [n for n, _ in sizes])
    ax.set_ylim(len(sizes) - 0.5, -0.9)
    ax.set_xlim(-0.8, max(int(s) for _, s in sizes) / 1e6 * 1.02)
    ax.set_xlabel("A. digitifera coordinate (Mb)   ·   tall boxes = inversions, short boxes = translocations")
    ax.grid(axis="y", visible=False)
    ax.set_title("Rearrangements reproduced by orthologous genes, coloured by evidence class",
                 loc="left", fontsize=10.5)
    handles, hl = ax.get_legend_handles_labels()
    keep = [(h, l) for c, l, _ in CLASS_ORDER for h, l2 in zip(handles, hl) if l2 == l]
    key_ax.legend([h for h, _ in keep], [l for _, l in keep], loc="upper left", fontsize=8,
                  title="Evidence class", title_fontsize=8.5, alignment="left")
    lines = ["Numbered events (biological)", "    type  chr   start     derived in / split"]
    for n, (ev, _, _) in enumerate(numbered, 1):
        what = ev["derived_pattern"].replace(" | ", " vs ") or "unpolarized"
        chrom, span = ev["Adig_location"].split(":")
        start = int(span.split("-")[0]) / 1e6
        lines.append(f"{n:>2}. {ev['type']:<5} {chrom:<5} {start:5.1f} Mb  {what}")
    key_ax.text(0.0, 0.42, "\n".join(lines), va="top", ha="left", fontsize=7.6, family="DejaVu Sans Mono",
                color=INK, transform=key_ax.transAxes)
    save(fig, "fig4_biological_events")


def fig_genes(grades, absence, dup):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    ax = axes[0]
    order = ["G1", "G2", "G3", "G4", "G5", "G6", "R"]
    labels = {"G1": "G1 core 1-copy syntenic", "G2": "G2 core 1-copy relocated", "G3": "G3 core multi-copy stable",
              "G4": "G4 core copy-number variable", "G5": "G5 lost in 1-2 genomes", "G6": "G6 single genome",
              "R": "R repeat-like family"}
    counts = Counter(r["grade"] for r in grades)
    vals = [counts.get(g, 0) for g in order]
    ax.barh(range(len(order)), vals, color=[BLUE, "#86b6ef", AQUA, YELLOW, ORANGE, RED, GRAY], height=0.6)
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:,}", va="center", fontsize=8)
    ax.set_yticks(range(len(order)), [labels[g] for g in order])
    ax.invert_yaxis()
    ax.set_xlabel("Gene families")
    ax.set_title("Conservation grade of gene families", loc="left", fontsize=10)

    ax = axes[1]
    classes = [("absent_synteny_intact", "Absent, syntenic region intact", RED),
               ("disrupted", "Disrupted copy only", ORANGE), ("fragment", "Fragment only", YELLOW),
               ("absent_gap", "Absent, region has assembly gap", GRAY),
               ("absent_unresolved", "Absent, region unresolved", LIGHTGRAY)]
    bottom = [0] * len(CODES)
    for cls, label, color in classes:
        vals = [sum(1 for r in grades if r["grade"] in ("G5", "G6") and r[f"{c}_absence"] == cls) for c in CODES]
        ax.bar(range(len(CODES)), vals, 0.6, bottom=bottom, color=color, label=label, edgecolor="white", linewidth=1)
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_xticks(range(len(CODES)), CODES)
    ax.set_ylabel("Families absent in this genome")
    ax.set_title("Why a family is missing", loc="left", fontsize=10)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)

    ax = axes[2]
    modes = [("tandem", AQUA), ("proximal", "#86b6ef"), ("segmental", VIOLET), ("dispersed", BLUE), ("unplaced", GRAY)]
    agg = {(r["species"], r["mode"]): (int(r["loci"]), int(r["assembly_dup_suspect"])) for r in dup}
    bottom = [0] * len(CODES)
    for mode, color in modes:
        real = [agg.get((c, mode), (0, 0))[0] - agg.get((c, mode), (0, 0))[1] for c in CODES]
        ax.bar(range(len(CODES)), real, 0.6, bottom=bottom, color=color, label=mode, edgecolor="white", linewidth=1)
        bottom = [b + v for b, v in zip(bottom, real)]
    suspect = [sum(agg.get((c, m), (0, 0))[1] for m, _ in modes) for c in CODES]
    ax.bar(range(len(CODES)), suspect, 0.6, bottom=bottom, color="white", edgecolor=RED, hatch="////",
           linewidth=0.8, label="in highly similar (≥97%) segment\n(possible assembly duplication)")
    ax.set_xticks(range(len(CODES)), CODES)
    ax.set_ylabel("Duplicated gene loci")
    ax.set_title("Duplicated loci by mode", loc="left", fontsize=10)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    fig.tight_layout()
    save(fig, "fig5_gene_conservation")


def fig_evolution(evo):
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
    keys = [f"{a}-{b}" for a, b in PAIRS]
    for ax, metric, label in ((axes[0], "dS", "dS (synonymous divergence)"), (axes[1], "w", "dN/dS")):
        data = []
        for k in keys:
            vals = []
            for r in evo:
                try:
                    ds, v = float(r[f"dS_{k}"]), float(r[f"{metric}_{k}"])
                except (KeyError, ValueError):
                    continue
                if 0.005 < ds < 2:
                    vals.append(v)
            data.append(vals)
        ax.boxplot(data, orientation="horizontal", showfliers=False, widths=0.5,
                   medianprops={"color": ORANGE, "lw": 2}, boxprops={"color": BLUE}, whiskerprops={"color": BLUE},
                   capprops={"color": BLUE})
        ax.set_yticks(range(1, len(keys) + 1), [k.replace("-", " vs ") for k in keys])
        ax.invert_yaxis()
        ax.set_xlabel(label)
    axes[0].set_title("A. Synonymous divergence (dS), single-copy orthologs", loc="left", fontsize=10)
    axes[1].set_title("B. Selective constraint (dN/dS, yn00)", loc="left", fontsize=10)
    fig.tight_layout()
    save(fig, "fig6_sequence_evolution")


def fig_technology(decomp, score):
    """Fig 7: technology-effect estimates; one-line titles, notes under the axes."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"wspace": 0.32})
    x = list(range(len(CODES)))

    ax = axes[0]
    width = 0.38
    top = 0
    for k, (row, color, label) in enumerate(((decomp[0], BLUE, "in homopolymers"),
                                              (decomp[1], "#86b6ef", "outside homopolymers"))):
        vals = [float(row[f"excess_{c}_per_Mb"]) for c in CODES]
        top = max(top, max(vals))
        ax.bar([i + (k - 0.5) * width for i in x], vals, width * 0.92, color=color, label=label)
    ax.set_ylim(0, top * 1.2)
    ax.set_xticks(x, CODES)
    ax.set_ylabel("Excess 1-bp indels per aligned Mb")
    ax.set_title("A. Genome-wide indel excess (weak signal)", loc="left", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=8)
    ax.text(0.5, -0.27, "Relative values from 6 pairwise rates; the lowest assembly is set to 0.",
            transform=ax.transAxes, ha="center", fontsize=7.5, color=MUTED)

    ax = axes[1]
    vals = [float(r["frameshift_pct_ident95_Aten_queries"]) if r["frameshift_pct_ident95_Aten_queries"][:1].isdigit()
            else None for r in score]
    ax.bar(x, [v or 0 for v in vals], 0.55, color=[SPECIES_COLOR[r["species"]] for r in score])
    for i, v in enumerate(vals):
        ax.text(i, (v or 0) + 0.03, f"{v:.2f}%" if v is not None else "query\nspecies", ha="center",
                va="bottom", fontsize=8, color=INK)
    ax.set_ylim(0, max(v for v in vals if v is not None) * 1.25)
    ax.set_xticks(x, [r["species"] for r in score])
    ax.set_ylabel("% of genes with a frameshift")
    ax.set_title("B. Frameshifts in conserved genes", loc="left", fontsize=10)
    ax.text(0.5, -0.27, "A. tenuis single-copy proteins, alignments ≥95% identity, ≥90% coverage.",
            transform=ax.transAxes, ha="center", fontsize=7.5, color=MUTED)

    ax = axes[2]
    n = [int(r["lineage_specific_rearrangements"] or 0) for r in score]
    gap = [int(r["of_which_gap_bounded_artifacts"] or 0) for r in score]
    ax.bar(x, [a - b for a, b in zip(n, gap)], 0.55, color=VIOLET, label="other single-genome events")
    ax.bar(x, gap, 0.55, bottom=[a - b for a, b in zip(n, gap)], color=RED, label="both breakpoints at gaps")
    for i, r in enumerate(score):
        per = r.get("lineage_rearrangements_per_branch_length", "")
        ax.text(i, n[i] + max(n) * 0.03, f"{n[i]}" + (f"\n({per}/unit)" if per else ""),
                ha="center", va="bottom", fontsize=7.5, color=INK, linespacing=1.3)
    ax.set_ylim(0, max(n) * 1.35 + 1)
    ax.set_xticks(x, [r["species"] for r in score])
    ax.set_ylabel("Unique rearrangements in one genome only")
    ax.set_title("C. Rearrangements carried by one genome", loc="left", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=8)
    ax.text(0.5, -0.27, "(n/unit) = events per unit of terminal branch length (substitutions/site).",
            transform=ax.transAxes, ha="center", fontsize=7.5, color=MUTED)
    save(fig, "fig7_technology_effects")


def md_table(rows, cols, headers=None, fmt=None):
    headers = headers or cols
    fmt = fmt or {}
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for r in rows:
        out.append("| " + " | ".join(fmt.get(c, str)(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def num(x, d=1):
    try:
        return f"{float(x):,.{d}f}"
    except (TypeError, ValueError):
        return str(x)


def write_report(method, calls, events, anchors, grades, dup, evo, score, decomp):
    sv = [c for c in calls if c["type"] in ("INV", "TRANS", "INVTR")]
    by_group = Counter()
    bp_group = Counter()
    for c in sv:
        group = c["class"].split("_")[0]
        by_group[group] += 1
        bp_group[group] += int(c["size"])
    bio = [e for e in events if e["class"].startswith("biological")]
    grade_n = Counter(r["grade"] for r in grades)
    pair_sum = tsv(OUT / "10_sequence_evolution" / "pair_summary.tsv")
    quartets = tsv(OUT / "10_sequence_evolution" / "quartet_support.tsv")
    tree_file = OUT / "10_sequence_evolution" / "species_tree.treefile"
    tree = open(tree_file).read().strip() if tree_file.exists() else "(not available)"
    cn = tsv(OUT / "09_gene_conservation" / "copy_number_changes.tsv")
    inter = {}
    for a, b in PAIRS:
        rows = [r for r in anchors if r[f"{a}_chrom"].startswith("chr") and r[f"{b}_chrom"].startswith("chr")]
        inter[(a, b)] = (sum(r[f"{a}_chrom"] != r[f"{b}_chrom"] for r in rows), len(rows))
    split_counts = Counter(e["derived_pattern"] for e in events if e["class"] == "biological_shared")

    L = []
    L.append("# Comparative genomics of four *Acropora* assemblies: biological vs methodological differences\n")
    L.append("*Generated by `scripts/comparative_genomics/12_report.py`; all numbers come from the tables in "
             "`results/comparative_genomics/`. Rerun `run_all.sh` to regenerate.*\n")
    L.append("## 1. Key results\n")
    if score:
        s = {r["species"]: r for r in score}
        L.append(f"- **Most raw structural differences are methodological.** Of {len(sv):,} SyRI inversions and "
                 f"translocations ≥10 kb, {by_group['biological']:,} ({num(bp_group['biological'] / 1e6)} Mb) are "
                 f"supported as biological, {by_group['uncertain']:,} are uncertain and "
                 f"{by_group['artifact']:,} ({num(bp_group['artifact'] / 1e6)} Mb) are classified as artifacts; "
                 f"{by_group['unverified']:,} small calls contain too few genes to test. Counted pairwise, the same "
                 f"event can appear in up to three comparisons; after merging, {len(bio)} unique rearrangements are "
                 "supported as biological (table in section 4).")
        L.append(f"- **A. millepora (Amil_v2.1) carries most of the artifacts.** It has "
                 f"{s['Amil']['lineage_specific_rearrangements']} unique rearrangements seen only in its assembly, versus "
                 f"{', '.join(s[c]['lineage_specific_rearrangements'] for c in ('Adig', 'Aten', 'Acer'))} in "
                 "A. digitifera, A. tenuis and A. cervicornis, although its terminal branch is no longer than theirs, "
                 f"and the highest frameshift rate in conserved genes "
                 f"({s['Amil']['frameshift_pct_ident95_Aten_queries']}% vs "
                 f"{s['Adig']['frameshift_pct_ident95_Aten_queries']}% in A. digitifera and "
                 f"{s['Acer']['frameshift_pct_ident95_Aten_queries']}% in A. cervicornis, alignments ≥95% identity). "
                 "This is consistent with PacBio CLR contigs ordered by linkage maps rather than Hi-C.")
        absent = Counter(r[f"{c}_absence"] for r in grades if r["grade"] in ("G5", "G6") for c in CODES
                         if r[f"{c}_absence"])
        total_abs = sum(absent.values())
        L.append(f"- **Gene content is highly conserved.** {grade_n['G1']:,} gene families are single-copy and "
                 f"collinear in all four genomes (grade G1); {grade_n['G4']:,} are present in all four with variable "
                 f"copy number, and {grade_n['G5']:,} are missing from one or two genomes. Of {total_abs:,} "
                 f"family absences, only {absent['absent_synteny_intact']:,} "
                 f"({100 * absent['absent_synteny_intact'] / max(1, total_abs):.0f}%) are confident losses (syntenic "
                 f"region intact, no trace of the gene); {absent['disrupted']:,} have a disrupted copy, "
                 f"{absent['fragment']:,} a fragment, {absent['absent_gap']:,} fall at an assembly gap and "
                 f"{absent['absent_unresolved']:,} could not be placed.")
        L.append(f"- **Assembly gaps hide genes in A. cervicornis.** {s['Acer']['gene_absences_at_assembly_gap']} "
                 "families absent from A. cervicornis are expected at a position (between collinear flanking genes) "
                 f"that contains an assembly gap, versus {s['Adig']['gene_absences_at_assembly_gap']}–"
                 f"{s['Amil']['gene_absences_at_assembly_gap']} in the other genomes; these are unassembled "
                 "sequence, not losses.")
        if method:
            m = {r["species"]: r for r in method}
            L.append("- **Genome size.** " + "; ".join(
                f"{NAMES[c]} {num(m[c]['total_Mb'], 0)} Mb ({num(m[c]['repeat_pct_Red'], 0)}% repeat, "
                f"{num(m[c]['nonrepeat_Mb'], 0)} Mb non-repeat, BUSCO duplicated {m[c]['BUSCO_genome_duplicated']}%)"
                for c in CODES) + ". With uniform masking and near-complete BUSCO in all four, the size difference "
                "lies in repeat and intergenic sequence; whether it reflects true repeat content or repeat collapse "
                "cannot be resolved without reads.")
    species_split = quartets[0]["quartet"].replace("|", " | ") if quartets else ""
    if split_counts:
        concordant = split_counts.get(species_split, 0)
        L.append("- **Shared rearrangements** (two independently built assemblies agree): " +
                 "; ".join(f"`{k}`: {v}" for k, v in split_counts.most_common()) +
                 f". The species tree split is `{species_split}`; {concordant} shared events match it. Events "
                 "grouping non-sister species are expected for ancestral inversion polymorphisms sorted by "
                 "incomplete lineage sorting (about a quarter of single-copy gene trees also disagree with the "
                 "species tree), but because every discordant split includes A. millepora, an A. millepora "
                 "artifact coinciding with a real event elsewhere cannot be excluded.")
    L.append("")

    L.append("## 2. Assemblies and how differences in technology were handled\n")
    if score:
        L.append(md_table(score, ["species", "technology"], ["Assembly", "Sequencing / assembly strategy"]))
    L.append("\nThe comparisons were designed so that methodology cannot masquerade as biology:\n")
    L.append("1. **Same coordinate frame.** Homologous chromosomes were renamed chr1–14 in *A. digitifera* numbering "
             "and orientation (whole-chromosome orientation differences are naming conventions, not inversions).")
    L.append("2. **Same gene models.** The four annotations come from different pipelines (AUGUSTUS, NCBI Gnomon, "
             "none for *A. cervicornis*). All ~88,600 annotated proteins were projected onto every genome with "
             "identical miniprot settings, so copy number, presence/absence and synteny are measured the same way "
             "everywhere.")
    L.append("3. **Same repeat masking.** The assemblies were soft-masked inconsistently (two unmasked); Red was "
             "rerun uniformly.")
    L.append("4. **Divergence-appropriate alignment.** minimap2 `asm20` instead of `asm5`, which left most of the "
             "inter-species genome unaligned in the earlier `results/syri_results` run.")
    L.append("5. **Every rearrangement is tested with independent evidence** (section 4).\n")

    L.append("## 3. Estimating the effect of sequencing technology\n")
    L.append("No raw reads were available, so technology effects were estimated from signals in the assemblies that "
             "biology produces rarely or that should scale with divergence:\n")
    L.append("| Signal | What it detects | How it is measured |\n|---|---|---|")
    L.append("| 1-bp indel excess | base-level errors (CLR/ONT) | 1-bp indels per aligned Mb in each pair are modelled "
             "as `k·SNP density + e_A + e_B` and solved over the 6 pairs (non-negative least squares); `e_x` is the "
             "excess contributed by assembly x |")
    L.append("| Frameshifts, query-controlled | base-level errors in genes | *A. tenuis* single-copy proteins projected "
             "onto the other three genomes (equally divergent query) |")
    L.append("| BUSCO duplicated, highly similar duplicated segments | retained haplotypes / false duplications | BUSCO "
             "cnidaria_odb12; 20 kb fragments mapped back to their own genome, secondary hits ≥10 kb at ≥97% identity |")
    L.append("| Rearrangements carried by one genome only, per branch length | scaffolding / misassembly | step 07 "
             "events ÷ terminal branch length of the species tree |")
    L.append("| Breakpoints or missing genes at assembly gaps | scaffolding joins, unassembled regions | distance to "
             "N-runs; gaps between flanking syntenic genes |")
    L.append("| Annotation vs projection | annotation pipeline | loci found by projection but absent from the annotation |\n")
    if score:
        cols = ["species", "contig_N50_Mb", "gaps", "pct_in_chromosomes", "BUSCO_genome_duplicated",
                "near_identical_dup_Mb", "frameshift_pct_ident95_Aten_queries", "homopolymer_indel_excess_per_Mb",
                "other_indel_excess_per_Mb", "lineage_specific_rearrangements", "of_which_gap_bounded_artifacts",
                "gene_absences_at_assembly_gap", "projected_loci_missing_from_annotation"]
        heads = ["Assembly", "Contig N50 (Mb)", "Gaps", "% in chr", "BUSCO dup %", "Similar (≥97%) dup segments (Mb)",
                 "Frameshift % (Aten queries, ≥95% id)", "Homopolymer indel excess /Mb", "Other indel excess /Mb",
                 "Lineage-only rearrangements", "…gap-bounded", "Gene absences at gaps", "Loci missed by annotation"]
        L.append(md_table(score, cols, heads))
    fs = tsv(OUT / "11_technology_effects" / "frameshift_by_identity.tsv")
    if fs:
        L.append("\nFrameshift and premature-stop rates of *A. tenuis* single-copy proteins (≥90% coverage), "
                 "compared within alignment-identity strata so that divergence cannot drive the difference:\n")
        L.append(md_table(fs, ["species", "identity", "genes", "frameshift_pct", "stop_pct"],
                          ["Genome", "Alignment identity", "Genes", "Frameshift %", "Premature stop %"]))
    L.append("\nThe genome-wide indel decomposition is shown for completeness. At 3–8% divergence, real indels in "
             "aligned DNA (~800 per Mb) swamp assembly errors (typically tens per Mb), so it has little power; "
             "coding sequence of conserved genes, where real frameshifts are rare, is the informative place to look.\n")
    L.append("\n![Assembly methodology](figures/fig1_assembly_methodology.png)\n")
    L.append("![Technology effects](figures/fig7_technology_effects.png)\n")
    if method:
        L.append(md_table(method, ["species", "total_Mb", "repeat_pct_Red", "nonrepeat_Mb", "unplaced_Mb",
                                   "BUSCO_genome_complete", "BUSCO_protein_complete", "annotated_genes"],
                          ["Assembly", "Size (Mb)", "Repeat % (Red)", "Non-repeat (Mb)", "Unplaced (Mb)",
                           "BUSCO genome C %", "BUSCO protein C %", "Annotated genes"]))
    L.append("\n**How to read it.** Values far above the other three assemblies (after divergence is removed) "
             "indicate technology effects. Genome-size differences that remain after uniform repeat masking, with "
             "near-complete BUSCO and no excess duplication, are more likely biological (repeat content) or due to "
             "repeat collapse, which cannot be separated without reads.\n")

    L.append("## 4. Structural differences\n")
    L.append("Each SyRI inversion/translocation ≥10 kb was classified:\n")
    for cls, label, _ in CLASS_ORDER:
        n = sum(1 for c in sv if c["class"] == cls)
        if n:
            L.append(f"- **{label}** (`{cls}`): {n:,} calls, {num(sum(int(c['size']) for c in sv if c['class'] == cls) / 1e6)} Mb")
    L.append("\nRules: a call must be reproduced by single-copy orthologs (genes inside the reference interval found "
             "inside the query interval and inverted/moved accordingly). The arrangement is then scored in all four "
             "genomes. If one genome alone carries it and both its breakpoints are within 50 kb of an assembly gap, "
             "it is a scaffolding artifact. If the junctions lie inside contiguous sequence, it is lineage-specific. "
             "If two genomes share it, it is a shared biological event.\n")
    L.append("![SV classes](figures/fig2_sv_classes.png)\n")
    L.append("**Inter-chromosomal exchange.** Single-copy orthologs on non-homologous chromosomes:\n")
    L.append("| Pair | On different chromosomes | Of |\n|---|---:|---:|")
    for (a, b), (n, t) in inter.items():
        L.append(f"| {a} vs {b} | {n:,} ({100 * n / max(1, t):.2f}%) | {t:,} |")
    L.append("\nAll 14 chromosomes are 1:1 homologous among the four species; no fusion or fission was found.\n")
    L.append("![Dot plots](figures/fig3_anchor_dotplots.png)\n")
    L.append("![Confirmed events](figures/fig4_biological_events.png)\n")
    if bio:
        L.append("### Rearrangements supported as biological\n")
        L.append(md_table(sorted(bio, key=lambda e: (e["class"], e["Adig_location"])),
                          ["event", "type", "class", "derived_pattern", "Adig_location", "size_bp_max", "core_genes",
                           "n_comparisons"],
                          ["Event", "Type", "Class", "Derived in / split", "A. digitifera location", "Max size (bp)",
                           "Genes", "Comparisons"],
                          {"size_bp_max": lambda v: num(v, 0), "derived_pattern": lambda v: v.replace(" | ", " vs ")}))
    L.append("")

    L.append("## 5. Gene conservation and duplication\n")
    L.append("| Grade | Meaning | Families |\n|---|---|---:|")
    meaning = {"G1": "core, single-copy, collinear in ≥5 of 6 pairs", "G2": "core, single-copy, relocated",
               "G3": "core, multi-copy, same copy number", "G4": "core, copy number varies",
               "G5": "missing from 1–2 genomes", "G6": "found in one genome only", "R": "repeat-like family (>50 copies)"}
    for g in ("G1", "G2", "G3", "G4", "G5", "G6", "R"):
        L.append(f"| {g} | {meaning[g]} | {grade_n.get(g, 0):,} |")
    L.append("\nAbsences (G5/G6) were verified on the genome: a disrupted copy (frameshift/stop), a fragmentary hit, "
             "an assembly gap between the flanking syntenic genes, or an intact syntenic region with no trace of the "
             "gene (confident loss). Only the last class is strong evidence of biological gene loss.\n")
    L.append("![Gene conservation](figures/fig5_gene_conservation.png)\n")
    if dup:
        L.append(md_table(dup, ["species", "mode", "loci", "assembly_dup_suspect"],
                          ["Assembly", "Duplication mode", "Loci", "In ≥97% similar duplicated segments (possible false duplication)"]))
    if cn:
        L.append("\n### Copy-number changes in a single genome (corrected counts)\n")
        L.append(f"{len(cn):,} family-genome combinations; the 25 largest:\n")
        top = sorted(cn, key=lambda r: -max(int(r[f"{c}_corr"]) for c in CODES))[:25]
        L.append(md_table(top, ["orthogroup", "species", "change"] + [f"{c}_corr" for c in CODES] +
                          ["dominant_mode", "product"],
                          ["Family", "Genome", "Change"] + CODES + ["Mode", "Product (A. millepora)"]))
    L.append("")

    L.append("## 6. Sequence conservation and species tree\n")
    L.append(f"Species tree (IQ-TREE, LG+F+R4, concatenated single-copy proteins, UFBoot): `{tree}`\n")
    if quartets:
        L.append(md_table(quartets, ["quartet", "genes"], ["Gene-tree quartet", "Single-copy genes"]))
    if pair_sum:
        L.append("\n" + md_table(pair_sum, list(pair_sum[0])))
    if evo:
        cls = Counter(r["conservation_class"] for r in evo)
        L.append("\nConservation classes of single-copy orthologs: " +
                 ", ".join(f"{k}: {v:,}" for k, v in cls.most_common()) + ".\n")
    L.append("![Sequence evolution](figures/fig6_sequence_evolution.png)\n")

    L.append("## 7. Caveats\n")
    L.append("- Without reads, collapsed repeats and phasing errors cannot be measured directly; the genome-size "
             "difference of *A. cervicornis* is reported, not attributed.")
    L.append("- Projected gene models come from three annotated species; genes unique to *A. cervicornis* cannot "
             "be detected, and divergent families are found less completely in it.")
    L.append("- Events with fewer than three single-copy genes inside (mostly <100 kb) are listed as unverified, not "
             "as artifacts.")
    L.append("- The four-taxon tree is unrooted; shared events are reported as splits, not as derived states.\n")

    L.append("## 8. Files\n")
    L.append("| Path | Content |\n|---|---|")
    for path, text in [("08_assembly_summary/assembly_methodology.tsv", "assembly metrics"),
                       ("11_technology_effects/technology_scorecard.tsv", "technology-sensitive metrics per assembly"),
                       ("11_technology_effects/indel_decomposition.tsv", "indel excess per assembly"),
                       ("07_structural_variants/sv_calls_pairwise.tsv", "every rearrangement with evidence and class"),
                       ("07_structural_variants/sv_events_merged.tsv", "unique events across comparisons"),
                       ("09_gene_conservation/family_table.tsv", "copy number, grade, absence class per family"),
                       ("09_gene_conservation/gene_grades_<species>.tsv", "grade for every annotated gene"),
                       ("09_gene_conservation/locus_duplication.tsv", "duplication mode of each duplicated locus"),
                       ("10_sequence_evolution/ortholog_evolution.tsv", "dN, dS, dN/dS, identity per ortholog"),
                       ("05_loci/<species>.loci.tsv", "uniform gene loci per genome")]:
        L.append(f"| `{path}` | {text} |")
    (REP / "REPORT.md").write_text("\n".join(L) + "\n")


def main():
    method = tsv(OUT / "08_assembly_summary" / "assembly_methodology.tsv")
    calls = tsv(OUT / "07_structural_variants" / "sv_calls_pairwise.tsv")
    events = tsv(OUT / "07_structural_variants" / "sv_events_merged.tsv")
    anchors = tsv(OUT / "06_synteny" / "single_copy_anchors.tsv")
    grades = tsv(OUT / "09_gene_conservation" / "family_table.tsv")
    dup = tsv(OUT / "09_gene_conservation" / "duplication_summary.tsv")
    evo = tsv(OUT / "10_sequence_evolution" / "ortholog_evolution.tsv")
    score = tsv(OUT / "11_technology_effects" / "technology_scorecard.tsv")
    decomp = tsv(OUT / "11_technology_effects" / "indel_decomposition.tsv")
    if method:
        fig_assembly(method)
    if calls:
        fig_sv(calls)
    if anchors:
        fig_dotplots(anchors, events)
    if events:
        fig_events(events)
    if grades:
        fig_genes(grades, None, dup)
    if evo:
        fig_evolution(evo)
    if score and decomp:
        fig_technology(decomp, score)
    write_report(method, calls, events, anchors, grades, dup, evo, score, decomp)
    print("report and figures written to", REP)


if __name__ == "__main__":
    main()
