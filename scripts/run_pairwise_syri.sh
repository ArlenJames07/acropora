#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
genome_dir="${project_dir}/genomes/chromosomes_only"
result_dir="${project_dir}/results/syri_results"
threads="${THREADS:-14}"

mkdir -p "${result_dir}"

run_pair() {
    local reference_name="$1"
    local query_name="$2"
    local pair="${reference_name}_vs_${query_name}"
    local pair_dir="${result_dir}/${pair}"
    local reference="${genome_dir}/${reference_name}.chromosomes.fa"
    local query="${genome_dir}/${query_name}.chromosomes.fa"
    local bam="${pair_dir}/${pair}.sorted.bam"

    mkdir -p "${pair_dir}"
    printf 'Running %s\n' "${pair}"

    minimap2 -ax asm5 --eqx -t "${threads}" "${reference}" "${query}" \
        2> "${pair_dir}/minimap2.log" \
        | samtools sort -@ "${threads}" -o "${bam}" -
    samtools index -@ "${threads}" "${bam}"

    syri -c "${bam}" -F B --cigar \
        -r "${reference}" -q "${query}" \
        --nc "${threads}" \
        --dir "${pair_dir}" \
        --prefix "${pair}." \
        --samplename "${query_name}" \
        --lf "${pair}.syri.log"

    test -s "${pair_dir}/${pair}.syri.out" || {
        printf 'ERROR: SyRI did not create %s\n' "${pair_dir}/${pair}.syri.out" >&2
        return 1
    }
}

run_pair A_digitifera A_millepora
run_pair A_digitifera A_tenius
run_pair A_millepora A_tenius
