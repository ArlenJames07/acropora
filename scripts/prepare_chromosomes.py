#!/usr/bin/env python3
"""Extract Acropora chromosomes 1-14 and give them identical chrN names."""

from __future__ import annotations

import gzip
import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parent.parent
SOURCE_DIR = PROJECT / "genomes"
OUTPUT_DIR = SOURCE_DIR / "chromosomes_only"
ASSEMBLIES = ("A_digitifera", "A_millepora", "A_tenius")

# Mappings are source chromosome number -> A. digitifera homolog number.
# A. millepora homologs and orientations were inferred from dominant asm5
# whole-chromosome alignments. The other two assemblies already match.
CHROMOSOME_MAPS = {
    "A_digitifera": {number: number for number in range(1, 15)},
    "A_millepora": {
        1: 2,
        2: 3,
        3: 9,
        4: 4,
        5: 5,
        6: 12,
        7: 8,
        8: 1,
        9: 10,
        10: 6,
        11: 11,
        12: 13,
        13: 14,
        14: 7,
    },
    "A_tenius": {number: number for number in range(1, 15)},
}
REVERSE_COMPLEMENT = {
    "A_digitifera": set(),
    "A_millepora": {5, 6, 7, 10, 11, 12, 13, 14},
    "A_tenius": set(),
}
COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def chromosome_number(header: str) -> int | None:
    """Return chromosome number encoded by one of the two input header styles."""
    name = header.split()[0]
    match = re.fullmatch(r"chr(\d+)Ref", name)
    if match:
        return int(match.group(1))

    match = re.search(r"\bchromosome\s+(\d+)\b", header, flags=re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def extract(assembly: str, source: Path, destination: Path) -> None:
    records: dict[int, str] = {}
    source_number: int | None = None
    sequence: list[str] = []

    with gzip.open(source, "rt") as input_fasta:
        for line in input_fasta:
            if line.startswith(">"):
                if source_number is not None:
                    records[source_number] = "".join(sequence)
                number = chromosome_number(line[1:].strip())
                source_number = number if number is not None and 1 <= number <= 14 else None
                sequence = []
            elif source_number is not None:
                sequence.append(line.strip())
        if source_number is not None:
            records[source_number] = "".join(sequence)

    expected = set(range(1, 15))
    if set(records) != expected:
        destination.unlink(missing_ok=True)
        missing = sorted(expected - set(records))
        extra = sorted(set(records) - expected)
        raise ValueError(f"{source}: missing={missing}, unexpected={extra}")

    normalized: dict[int, str] = {}
    for original_number, target_number in CHROMOSOME_MAPS[assembly].items():
        sequence_text = records[original_number]
        if original_number in REVERSE_COMPLEMENT[assembly]:
            sequence_text = sequence_text.translate(COMPLEMENT)[::-1]
        normalized[target_number] = sequence_text

    with destination.open("w") as output_fasta:
        for target_number in range(1, 15):
            output_fasta.write(f">chr{target_number}\n")
            sequence_text = normalized[target_number]
            for start in range(0, len(sequence_text), 80):
                output_fasta.write(sequence_text[start : start + 80] + "\n")


def write_mapping_table() -> None:
    mapping_file = OUTPUT_DIR / "chromosome_mapping.tsv"
    with mapping_file.open("w") as table:
        table.write("assembly\tsource_chromosome\tnormalized_chromosome\torientation\n")
        for assembly in ASSEMBLIES:
            for source_number, target_number in CHROMOSOME_MAPS[assembly].items():
                orientation = (
                    "reverse_complemented"
                    if source_number in REVERSE_COMPLEMENT[assembly]
                    else "forward"
                )
                table.write(
                    f"{assembly}\tchr{source_number}\tchr{target_number}\t{orientation}\n"
                )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for assembly in ASSEMBLIES:
        source = SOURCE_DIR / f"{assembly}.fa.gz"
        destination = OUTPUT_DIR / f"{assembly}.chromosomes.fa"
        print(f"Preparing {destination.name} from {source.name}", flush=True)
        extract(assembly, source, destination)
    write_mapping_table()


if __name__ == "__main__":
    main()
