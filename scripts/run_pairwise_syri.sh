#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
genome_dir="${project_dir}/genomes/chromosomes_only"
result_dir="${project_dir}/results/syri_results"
threads="${THREADS:-14}"
force="${FORCE:-0}"
species=(A_digitifera A_millepora A_tenius A_cervicornis)
requested_pairs=("$@")

mkdir -p "${result_dir}"

pair_requested() {
    local pair="$1"
    local requested_pair

    if ((${#requested_pairs[@]} == 0)); then
        return 0
    fi
    for requested_pair in "${requested_pairs[@]}"; do
        if [[ "${requested_pair}" == "${pair}" ]]; then
            return 0
        fi
    done
    return 1
}

run_pair() {
    local reference_name="$1"
    local query_name="$2"
    local pair="${reference_name}_vs_${query_name}"
    local pair_dir="${result_dir}/${pair}"
    local reference="${genome_dir}/${reference_name}.chromosomes.fa"
    local query="${genome_dir}/${query_name}.chromosomes.fa"
    local bam="${pair_dir}/${pair}.sorted.bam"

    pair_requested "${pair}" || return

    mkdir -p "${pair_dir}"
    if [[ "${force}" != "1" \
        && -s "${bam}" \
        && -s "${bam}.bai" \
        && -s "${pair_dir}/${pair}.syri.out" \
        && -s "${pair_dir}/${pair}.syri.vcf" \
        && -s "${pair_dir}/${pair}.syri.summary" ]]; then
        printf 'Skipping completed comparison %s (set FORCE=1 to rerun)\n' "${pair}"
        return
    fi

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
        --lf syri.log

    test -s "${pair_dir}/${pair}.syri.out" || {
        printf 'ERROR: SyRI did not create %s\n' "${pair_dir}/${pair}.syri.out" >&2
        return 1
    }
}

for ((reference_index = 0; reference_index < ${#species[@]} - 1; reference_index++)); do
    for ((query_index = reference_index + 1; query_index < ${#species[@]}; query_index++)); do
        run_pair "${species[reference_index]}" "${species[query_index]}"
    done
done
