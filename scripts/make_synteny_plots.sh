#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
genome_dir="${project_dir}/genomes/chromosomes_only"
result_dir="${project_dir}/results/syri_results"
plot_dir="${result_dir}/plots"
species=(A_digitifera A_millepora A_tenius A_cervicornis)

mkdir -p "${plot_dir}"

write_genomes_file() {
    local output="$1"
    shift
    : > "${output}"
    local genome
    for genome in "$@"; do
        printf '%s\t%s\n' \
            "${genome_dir}/${genome}.chromosomes.fa" "${genome}" >> "${output}"
    done
}

printf '%s\n' chr{1..14} > "${plot_dir}/chromosome_order.txt"

plot_pair() {
    local reference_name="$1"
    local query_name="$2"
    local pair="${reference_name}_vs_${query_name}"
    local genomes_file="${plot_dir}/${pair}.genomes.txt"
    local syri_file="${result_dir}/${pair}/${pair}.syri.out"

    write_genomes_file "${genomes_file}" "${reference_name}" "${query_name}"
    local extension
    for extension in pdf png; do
        plotsr \
            --sr "${syri_file}" \
            --genomes "${genomes_file}" \
            --chrord "${plot_dir}/chromosome_order.txt" \
            -W 20 -H 8 -f 8 -d 600 \
            -o "${plot_dir}/${pair}.${extension}" \
            --lf "${plot_dir}/${pair}.${extension}.plotsr.log"
    done
}

for ((reference_index = 0; reference_index < ${#species[@]} - 1; reference_index++)); do
    for ((query_index = reference_index + 1; query_index < ${#species[@]}; query_index++)); do
        plot_pair "${species[reference_index]}" "${species[query_index]}"
    done
done

plot_multiple() {
    local output_stem="$1"
    local height="$2"
    shift 2
    local genomes=("$@")
    local genomes_file="${plot_dir}/${output_stem}.genomes.txt"
    local syri_arguments=()
    local index reference_name query_name pair

    write_genomes_file "${genomes_file}" "${genomes[@]}"
    for ((index = 0; index < ${#genomes[@]} - 1; index++)); do
        reference_name="${genomes[index]}"
        query_name="${genomes[index + 1]}"
        pair="${reference_name}_vs_${query_name}"
        syri_arguments+=(
            --sr "${result_dir}/${pair}/${pair}.syri.out"
        )
    done

    local extension
    for extension in pdf png; do
        plotsr \
            "${syri_arguments[@]}" \
            --genomes "${genomes_file}" \
            --chrord "${plot_dir}/chromosome_order.txt" \
            -W 20 -H "${height}" -f 8 -d 600 \
            -o "${plot_dir}/${output_stem}.${extension}" \
            --lf "${plot_dir}/${output_stem}.${extension}.plotsr.log"
    done
}

plot_multiple all_three_synteny 12 A_digitifera A_millepora A_tenius
plot_multiple all_four_synteny 16 "${species[@]}"

control_pair="A_tenius_vs_A_tenius_control"
control_syri="${result_dir}/${control_pair}/${control_pair}.syri.out"
if [[ -s "${control_syri}" ]]; then
    control_genomes="${plot_dir}/${control_pair}.genomes.txt"
    write_genomes_file "${control_genomes}" A_tenius A_tenius
    for extension in pdf png; do
        plotsr \
            --sr "${control_syri}" \
            --genomes "${control_genomes}" \
            --chrord "${plot_dir}/chromosome_order.txt" \
            -W 20 -H 8 -f 8 -d 600 \
            -o "${plot_dir}/${control_pair}.${extension}" \
            --lf "${plot_dir}/${control_pair}.${extension}.plotsr.log"
    done
fi
