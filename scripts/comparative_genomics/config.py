"""Shared paths, species table and helpers for the comparative_genomics workflow."""

from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
RESOURCES = PROJECT / "resources"
GENOMES = PROJECT / "genomes"
OUT = PROJECT / "results" / "comparative_genomics"
INPUTS = OUT / "01_inputs"
ENV_BIN = Path(os.environ.get("CONDA_ENV_BIN", "/home/rare/miniconda3/envs/acropora_cg/bin"))

# Reuse the chromosome homology/orientation table of the existing SyRI workflow,
# so every coordinate in this workflow is in A. digitifera chromosome numbering.
sys.path.insert(0, str(PROJECT / "scripts"))
from prepare_chromosomes import CHROMOSOME_MAPS, REVERSE_COMPLEMENT  # noqa: E402

# code, display name, key in CHROMOSOME_MAPS, genome, gff, proteins, annotation style
SPECIES = [
    ("Adig", "A. digitifera", "A_digitifera", "AdiV3.2Ref_genome.fa.gz",
     "AdiV3.2Ref.gff.gz", "AdiV3.2Ref_aa.fa.gz", "augustus"),
    ("Amil", "A. millepora", "A_millepora", "A_millepora_v2.1.fasta",
     "A_millepora_v2.1.gff", "A_millepora_v2.1_protein.faa", "refseq"),
    ("Aten", "A. tenuis", "A_tenius", "AteV2.2Ref_genome.fa.gz",
     "AteV2.2Ref.gff.gz", "AteV2.2Ref_aa.fa.gz", "augustus"),
    ("Acer", "A. cervicornis", "A_cervicornis", "GCA_964034985.1_jaAcrCerv1.1_genomic.fna.gz",
     None, None, None),
]
CODES = [s[0] for s in SPECIES]
ANNOTATED = [s[0] for s in SPECIES if s[4]]
NAMES = {s[0]: s[1] for s in SPECIES}

# How each assembly was produced; used only to explain methodological signal.
ASSEMBLY_METHOD = {
    "Adig": "AdiV3.2Ref (OIST); chromosome-scale, AUGUSTUS annotation",
    "Amil": "Amil_v2.1 (RefSeq GCF_013753865.1); chromosomes + 840 unplaced scaffolds, NCBI Gnomon annotation",
    "Aten": "AteV2.2Ref (OIST); chromosome-scale, AUGUSTUS annotation",
    "Acer": "jaAcrCerv1.1 (Darwin Tree of Life); PacBio HiFi + Arima Hi-C, no annotation",
}


def genome_fasta(code: str) -> Path:
    return INPUTS / "genomes" / f"{code}.fa"


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path: Path, full_header: bool = False):
    """Yield (name, sequence) records."""
    name, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].rstrip() if full_header else line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_fasta(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i + width] + "\n")


def ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def orthogroups_dir() -> Path:
    """Most recent OrthoFinder Orthogroups/ directory (a -b rerun nests it deeper)."""
    found = sorted((OUT / "04_orthology" / "orthofinder").rglob("Orthogroups/Orthogroups.tsv"),
                   key=lambda p: p.stat().st_mtime)
    if not found:
        raise SystemExit("OrthoFinder Orthogroups.tsv not found; run step 04 first")
    return found[-1].parent
