#!/usr/bin/env python3
"""Step 07 - separate biological rearrangements from methodological artifacts.

Candidate rearrangements are the SyRI calls (step 03) between homologous
chromosomes: inversions (INV), translocations (TRANS), inverted translocations
(INVTR) and duplications (DUP/INVDP). Each call is evaluated with evidence that
does not come from the DNA alignment itself:

1. Ortholog support. Single-copy orthologs present once in all four genomes
   (step 06) that lie inside the reference interval must map inside the query
   interval. Calls that orthologous genes do not reproduce are alignment
   artifacts (repeats, paralogs, divergent sequence).
2. Phylogenetic pattern. With 4 genomes, the arrangement of the event's genes
   relative to their flanking genes is scored in every genome as reference-like
   or query-like. If one genome alone carries the derived arrangement, the event
   is lineage-specific (or an artifact of that one assembly). If two genomes
   share it, two independently produced assemblies agree, which is strong
   evidence that the rearrangement is real.
3. Breakpoints at assembly gaps. Hi-C or linkage-map scaffolding errors flip or
   misplace whole contigs, so their breakpoints coincide with gaps. A
   lineage-specific event whose breakpoints both lie at gaps of the derived
   genome is classified as a probable scaffolding artifact; if the breakpoints
   lie inside contigs, the novel junction is spanned by contiguous sequence.
4. Repeat content (uniform Red masking) and size.

Outputs (07_structural_variants/):
  sv_calls_pairwise.tsv   every SyRI rearrangement >= MIN_SIZE with evidence and class
  sv_events_merged.tsv    confirmed events merged across comparisons (unique events)
  sv_class_summary.tsv    counts and bp per class and comparison
"""

from __future__ import annotations

import bisect
import glob
import re
from collections import Counter, defaultdict
from pathlib import Path

from config import CODES, INPUTS, OUT, ensure

WGA = OUT / "03_wga"
RED = OUT / "02_assembly_qc" / "red"
SYN = OUT / "06_synteny"
SV = ensure(OUT / "07_structural_variants")

MIN_SIZE = 10_000          # rearrangements smaller than this are only counted
GAP_WINDOW = 50_000        # breakpoint within this distance of a gap = gap-associated
FLANK_ANCHORS = 5          # flanking single-copy genes used on each side
FLANK_MAX_DIST = 3_000_000
SV_TYPES = {"INV", "TRANS", "INVTR", "DUP", "INVDP"}
SMALL = Counter()  # (pair, type) -> calls below MIN_SIZE
# A genome with more than this many times the median lineage-specific
# rearrangement rate (per unit terminal branch length) of the other genomes is
# treated as error-prone: its "lineage-specific" events become uncertain.
ERROR_PRONE_FOLD = 10
TREE = OUT / "10_sequence_evolution" / "species_tree.treefile"


class Intervals:
    """Sorted, merged intervals per chromosome with fast overlap queries."""

    def __init__(self, path: Path):
        self.data = defaultdict(list)
        if path.exists():
            with open(path) as handle:
                for line in handle:
                    c, s, e = line.split()[:3]
                    self.data[c].append((int(s), int(e)))
        self.starts = {}
        for c in self.data:
            self.data[c].sort()
            self.starts[c] = [s for s, _ in self.data[c]]

    def covered(self, chrom, start, end):
        ivs, starts, total = self.data.get(chrom, []), self.starts.get(chrom, []), 0
        i = max(bisect.bisect_right(starts, start) - 1, 0)
        while i < len(ivs) and ivs[i][0] < end:
            s, e = ivs[i]
            total += max(0, min(e, end) - max(s, start))
            i += 1
        return total

    def distance(self, chrom, pos):
        """Distance from pos to the nearest interval edge (gap = contig boundary)."""
        ivs, starts = self.data.get(chrom, []), self.starts.get(chrom, [])
        if not ivs:
            return float("inf")
        i = bisect.bisect_right(starts, pos)
        best = float("inf")
        for j in (i - 1, i):
            if 0 <= j < len(ivs):
                s, e = ivs[j]
                best = min(best, 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e)))
        return best


def load_anchors():
    """Return {species: {og: (chrom, mid, strand)}} and per-chrom sorted lists."""
    pos = {c: {} for c in CODES}
    with open(SYN / "single_copy_anchors.tsv") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            for c in CODES:
                mid = (int(row[f"{c}_start"]) + int(row[f"{c}_end"])) // 2
                pos[c][row["orthogroup"]] = (row[f"{c}_chrom"], mid, row[f"{c}_strand"])
    order = {c: defaultdict(list) for c in CODES}
    for c in CODES:
        for og, (chrom, mid, _) in pos[c].items():
            order[c][chrom].append((mid, og))
        for chrom in order[c]:
            order[c][chrom].sort()
    return pos, order


def anchors_in(order, sp, chrom, start, end):
    items = order[sp].get(chrom, [])
    i = bisect.bisect_left(items, (start, ""))
    out = []
    while i < len(items) and items[i][0] <= end:
        out.append(items[i][1])
        i += 1
    return out


def flanks(order, sp, chrom, start, end):
    items = order[sp].get(chrom, [])
    i = bisect.bisect_left(items, (start, ""))
    j = bisect.bisect_right(items, (end, "~"))
    left = [og for p, og in items[max(0, i - FLANK_ANCHORS):i] if start - p <= FLANK_MAX_DIST]
    right = [og for p, og in items[j:j + FLANK_ANCHORS] if p - end <= FLANK_MAX_DIST]
    return left, right


def spearman_sign(xs, ys):
    """Sign of the rank correlation (+1 collinear, -1 inverted, 0 undefined)."""
    if len(xs) < 2:
        return 0
    rx = {v: i for i, v in enumerate(sorted(xs))}
    ry = {v: i for i, v in enumerate(sorted(ys))}
    n = len(xs)
    mean = (n - 1) / 2
    cov = sum((rx[x] - mean) * (ry[y] - mean) for x, y in zip(xs, ys))
    return (cov > 0) - (cov < 0)


def majority_chrom(pos, sp, ogs):
    chroms = Counter(pos[sp][og][0] for og in ogs)
    return chroms.most_common(1)[0][0] if chroms else None


def placed_between(pos, sp, core, left, right):
    """Do the core genes sit next to their flanking genes in sp?

    With flanks on both sides the core must lie between them. Events that
    reach a chromosome end have flanks on one side only; the core must then lie
    on the flank's chromosome within FLANK_MAX_DIST of it.
    """
    if not core or not (left or right):
        return None
    chrom_core = majority_chrom(pos, sp, core)
    lp = [pos[sp][og][1] for og in left if pos[sp][og][0] == chrom_core]
    rp = [pos[sp][og][1] for og in right if pos[sp][og][0] == chrom_core]
    inside = sorted(pos[sp][og][1] for og in core if pos[sp][og][0] == chrom_core)
    mid = inside[len(inside) // 2]
    if left and right:
        if not lp or not rp:
            return False
        lo, hi = sorted((sorted(lp)[len(lp) // 2], sorted(rp)[len(rp) // 2]))
        return lo <= mid <= hi
    side = lp or rp
    if len(side) < 2:
        return False
    return min(min(abs(p - q) for q in (inside[0], inside[-1])) for p in side) <= FLANK_MAX_DIST


def state(pos, ref, sp, kind, core, rflank, qflank):
    """'R' (reference-like), 'Q' (query-like) or '?' for one genome.

    Inversions: all chromosomes were put in a common orientation (step 01), so
    the order of the core genes in sp is compared with the reference directly;
    flanking genes are not used because in fragmented assemblies they are often
    rearranged themselves. Translocations: are the core genes next to the
    reference flanks or next to the query flanks?
    """
    if kind == "INV":
        chrom = majority_chrom(pos, sp, core)
        on_chrom = [og for og in core if pos[sp][og][0] == chrom]
        if len(on_chrom) < max(3, 0.7 * len(core)):
            return "?"
        sign = spearman_sign([pos[ref][og][1] for og in on_chrom], [pos[sp][og][1] for og in on_chrom])
        return {1: "R", -1: "Q"}.get(sign, "?")
    in_r = placed_between(pos, sp, core, *rflank)
    in_q = placed_between(pos, sp, core, *qflank)
    if in_r and not in_q:
        return "R"
    if in_q and not in_r:
        return "Q"
    return "?"


def parse_syri(path):
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[10] not in SV_TYPES:
                continue
            yield {"rchr": f[0], "rs": int(f[1]), "re": int(f[2]), "qchr": f[5],
                   "qs": min(int(f[6]), int(f[7])), "qe": max(int(f[6]), int(f[7])),
                   "id": f[8], "type": f[10], "copy": f[11]}


def classify(ev, pos, order, gaps, repeats):
    R, Q = ev["ref"], ev["qry"]
    kind = "INV" if ev["type"] == "INV" else ev["type"]
    ev["size"] = max(ev["re"] - ev["rs"], ev["qe"] - ev["qs"]) + 1
    for sp, c, s, e, tag in ((R, ev["rchr"], ev["rs"], ev["re"], "ref"), (Q, ev["qchr"], ev["qs"], ev["qe"], "qry")):
        ev[f"{tag}_repeat_frac"] = repeats[sp].covered(c, s, e) / max(1, e - s)
        d1, d2 = gaps[sp].distance(c, s), gaps[sp].distance(c, e)
        ev[f"{tag}_gap_dist_kb"] = f"{min(d1, 1e9) / 1000:.0f},{min(d2, 1e9) / 1000:.0f}"
        ev[f"{tag}_gap_breakpoints"] = (d1 <= GAP_WINDOW) + (d2 <= GAP_WINDOW)

    inside = anchors_in(order, R, ev["rchr"], ev["rs"], ev["re"])
    qinside = set(anchors_in(order, Q, ev["qchr"], ev["qs"], ev["qe"]))
    core = [og for og in inside if og in qinside]
    ev["anchors_ref"], ev["anchors_confirmed"] = len(inside), len(core)
    ev["core"] = core
    ev["states"] = {}

    gap_bounded = max(ev["ref_gap_breakpoints"], ev["qry_gap_breakpoints"]) == 2
    repeat_rich = min(ev["ref_repeat_frac"], ev["qry_repeat_frac"]) >= 0.5

    if kind in ("DUP", "INVDP"):
        side = "qry" if ev["copy"] == "copygain" else "ref"
        if ev[f"{side}_gap_breakpoints"] == 2:
            return "duplication_gap_bounded_suspect", "both ends of the extra copy at assembly gaps"
        if repeat_rich:
            return "duplication_repeat", "duplicated segment is mostly repeat"
        return "duplication_candidate", "see gene copy-number analysis (step 08)"

    if len(inside) < 3:
        if gap_bounded:
            return "artifact_gap_bounded", "no gene anchors; both breakpoints at assembly gaps"
        if repeat_rich:
            return "unverified_repeat", "no gene anchors; region mostly repeat"
        return "unverified_no_genes", "fewer than 3 single-copy genes inside"
    if len(core) < 0.5 * len(inside) or len(core) < 3:
        return "artifact_alignment", "orthologous genes do not reproduce the rearrangement"

    rflank = flanks(order, R, ev["rchr"], ev["rs"], ev["re"])
    qflank = flanks(order, Q, ev["qchr"], ev["qs"], ev["qe"])
    states = {sp: state(pos, R, sp, kind, core, rflank, qflank) for sp in CODES}
    ev["states"] = states
    if states[R] != "R" or states[Q] != "Q":
        return "artifact_alignment", "gene order does not show the rearrangement between ref and query"

    alt = sorted(sp for sp, s in states.items() if s == "Q")
    ref_like = sorted(sp for sp, s in states.items() if s == "R")
    unknown = [sp for sp, s in states.items() if s == "?"]
    if unknown:
        ev["derived"] = ""
        if gap_bounded:
            return "uncertain_gap_bounded", "confirmed by genes but a third genome is uninformative and breakpoints at gaps"
        return "biological_unpolarized", "confirmed by orthologous genes; third-genome state unresolved"
    if len(alt) == 1 or len(ref_like) == 1:
        derived = alt[0] if len(alt) == 1 else ref_like[0]
        ev["derived"] = derived
        tag = "qry" if derived == Q else "ref"
        n_gap = ev[f"{tag}_gap_breakpoints"]
        if n_gap == 2:
            return "artifact_scaffolding", f"only {derived} differs and both its breakpoints lie at assembly gaps"
        if n_gap == 1:
            return "uncertain_one_gap", f"only {derived} differs; one of its breakpoints lies at an assembly gap"
        return "biological_lineage_specific", f"only {derived} differs; breakpoints inside contiguous sequence"
    # Canonical bipartition: the half containing the first genome of CODES first.
    halves = sorted((alt, ref_like), key=lambda h: min(CODES.index(c) for c in h))
    ev["derived"] = " | ".join("+".join(sorted(h, key=CODES.index)) for h in halves)
    return "biological_shared", f"arrangement splits the genomes {ev['derived']}; shared by two independently assembled genomes"


def main():
    pos, order = load_anchors()
    gaps = {c: Intervals(INPUTS / f"{c}.gaps.bed") for c in CODES}
    repeats = {c: Intervals(RED / f"{c}.bed") for c in CODES}

    calls, small = [], SMALL
    for path in sorted(glob.glob(str(WGA / "*_vs_*" / "*.syri.out"))):
        pair = Path(path).parent.name
        ref, qry = pair.split("_vs_")
        for ev in parse_syri(path):
            ev["ref"], ev["qry"], ev["pair"] = ref, qry, pair
            size = max(ev["re"] - ev["rs"], ev["qe"] - ev["qs"]) + 1
            if size < MIN_SIZE:
                small[(pair, ev["type"])] += 1
                continue
            ev["class"], ev["reason"] = classify(ev, pos, order, gaps, repeats)
            calls.append(ev)
        print(f"[07] {pair}: {sum(1 for c in calls if c['pair'] == pair)} calls >= {MIN_SIZE} bp", flush=True)

    flag_error_prone(calls, pos)


def write_calls(calls):
    small = SMALL
    columns = ["pair", "id", "type", "class", "reason", "size", "rchr", "rs", "re", "qchr", "qs", "qe",
               "anchors_ref", "anchors_confirmed", "ref_gap_breakpoints", "qry_gap_breakpoints",
               "ref_gap_dist_kb", "qry_gap_dist_kb", "ref_repeat_frac", "qry_repeat_frac"]
    with open(SV / "sv_calls_pairwise.tsv", "w") as out:
        out.write("\t".join(columns + ["derived"] + [f"state_{c}" for c in CODES]) + "\n")
        for ev in calls:
            vals = [f"{ev[c]:.3f}" if isinstance(ev[c], float) else str(ev[c]) for c in columns]
            vals.append(ev.get("derived", ""))
            vals += [ev["states"].get(c, "") for c in CODES]
            out.write("\t".join(vals) + "\n")

    with open(SV / "sv_class_summary.tsv", "w") as out:
        out.write("pair\ttype\tclass\tcount\tbp_ref\n")
        agg = defaultdict(lambda: [0, 0])
        for ev in calls:
            key = (ev["pair"], ev["type"], ev["class"])
            agg[key][0] += 1
            agg[key][1] += ev["re"] - ev["rs"] + 1
        for (pair, kind), n in small.items():
            agg[(pair, kind, f"small_<{MIN_SIZE // 1000}kb")][0] += n
        for key in sorted(agg):
            out.write("\t".join(key) + f"\t{agg[key][0]}\t{agg[key][1]}\n")


def flag_error_prone(calls, pos):
    """Downgrade lineage-specific events of assemblies with an implausible excess.

    Biological rearrangements accumulate roughly with time, so the number of
    unique events per genome is compared after dividing by its terminal branch
    length (species tree from step 10). Run step 07 again after step 10 to use it.
    """
    events = merge_events(calls, pos, write=False)
    branches = {}
    if TREE.exists():
        branches = {m.group(1): float(m.group(2))
                    for m in re.finditer(r"(Adig|Amil|Aten|Acer):([0-9.eE-]+)", TREE.read_text())}
    counts = Counter(e["derived"] for e in events if e["class"] == "biological_lineage_specific")
    rate = {c: (counts.get(c, 0) + 1) / branches[c] for c in CODES if c in branches}  # +1: avoid zero rates
    flagged = set()
    for c, r in rate.items():
        others = sorted(v for k, v in rate.items() if k != c)
        if others and r > ERROR_PRONE_FOLD * others[len(others) // 2]:
            flagged.add(c)
    with open(SV / "lineage_rate_check.tsv", "w") as out:
        out.write("species\tlineage_specific_events\tterminal_branch\tevents_per_branch_unit\terror_prone\n")
        for c in CODES:
            out.write(f"{c}\t{counts.get(c, 0)}\t{branches.get(c, '')}\t{rate.get(c, 0):.0f}\t{int(c in flagged)}\n")
    for ev in calls:
        if ev["class"] == "biological_lineage_specific" and ev.get("derived") in flagged:
            ev["class"] = "uncertain_error_prone_assembly"
            ev["reason"] = (f"only {ev['derived']} differs, but this assembly carries >{ERROR_PRONE_FOLD}x more "
                            "lineage-specific rearrangements per branch length than the others")
    write_calls(calls)
    merge_events(calls, pos)


def merge_events(calls, pos, write=True):
    """Collapse the same event seen in several pairwise comparisons.

    Two confirmed calls are the same event when they have the same type, the
    same phylogenetic pattern and share at least half of their core orthologs.
    """
    keep = [ev for ev in calls if ev["class"].startswith(("biological", "uncertain", "artifact_scaffolding"))]
    for ev in calls:
        ev.setdefault("derived", "")
    keep.sort(key=lambda ev: -len(ev["core"]))
    events = []
    for ev in keep:
        kind = "INV" if ev["type"] == "INV" else "TRANS"
        core = set(ev["core"])
        for event in events:
            if event["kind"] == kind and len(core & event["core"]) >= 0.5 * min(len(core), len(event["core"])):
                event["calls"].append(ev)
                event["core"] |= core
                break
        else:
            events.append({"kind": kind, "core": core, "calls": [ev]})

    priority = ["biological_shared", "biological_lineage_specific", "biological_unpolarized",
                "uncertain_error_prone_assembly", "uncertain_one_gap", "uncertain_gap_bounded",
                "artifact_scaffolding"]
    merged = []
    for event in sorted(events, key=lambda e: (-len(e["calls"]), -len(e["core"]))):
        classes = Counter(c["class"] for c in event["calls"])
        # An event counts as artifact only if every comparison says so.
        cls = min(classes, key=priority.index)
        best = next(c for c in event["calls"] if c["class"] == cls)
        merged.append({**event, "class": cls, "derived": best["derived"], "best": best, "classes": classes})
    if not write:
        return merged

    with open(SV / "sv_events_merged.tsv", "w") as out:
        out.write("event\ttype\tclass\tderived_pattern\tn_comparisons\tcomparisons\tcore_genes\t"
                  "Adig_location\tsize_bp_max\tclasses_seen\treason\n")
        for n, event in enumerate(merged, 1):
            adig = [pos["Adig"][og] for og in event["core"]]
            chrom = Counter(a[0] for a in adig).most_common(1)[0][0]
            ps = [a[1] for a in adig if a[0] == chrom]
            out.write(f"SV{n:04d}\t{event['kind']}\t{event['class']}\t{event['derived']}\t{len(event['calls'])}\t"
                      f"{','.join(sorted(c['pair'] for c in event['calls']))}\t{len(event['core'])}\t"
                      f"{chrom}:{min(ps)}-{max(ps)}\t{max(c['size'] for c in event['calls'])}\t"
                      f"{';'.join(f'{k}:{v}' for k, v in event['classes'].items())}\t{event['best']['reason']}\n")
    return merged


if __name__ == "__main__":
    main()
