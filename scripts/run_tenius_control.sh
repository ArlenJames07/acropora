#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
genome_dir="${project_dir}/genomes/chromosomes_only"
result_dir="${project_dir}/results/syri_results"
threads="${THREADS:-14}"
name="A_tenius"
pair="A_tenius_vs_A_tenius_control"
pair_dir="${result_dir}/${pair}"
genome="${genome_dir}/${name}.chromosomes.fa"
bam="${pair_dir}/${pair}.sorted.bam"

mkdir -p "${pair_dir}"
printf 'Running %s\n' "${pair}"

minimap2 -ax asm5 --eqx -t "${threads}" "${genome}" "${genome}" \
    2> "${pair_dir}/minimap2.log" \
    | samtools sort -@ "${threads}" -o "${bam}" -
samtools index -@ "${threads}" "${bam}"

syri -c "${bam}" -F B --cigar \
    -r "${genome}" -q "${genome}" \
    --nc "${threads}" \
    --dir "${pair_dir}" \
    --prefix "${pair}." \
    --samplename "${name}_control" \
    --lf syri.log

test -s "${pair_dir}/${pair}.syri.out" || {
    printf 'ERROR: SyRI did not create %s\n' "${pair_dir}/${pair}.syri.out" >&2
    exit 1
}
