#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
genome_dir="${project_dir}/genomes/chromosomes_only"
result_dir="${project_dir}/results/syri_results"
plot_dir="${result_dir}/plots"

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

plot_pair A_digitifera A_millepora
plot_pair A_digitifera A_tenius
plot_pair A_millepora A_tenius

three_genomes="${plot_dir}/all_three.genomes.txt"
write_genomes_file "${three_genomes}" A_digitifera A_millepora A_tenius
for extension in pdf png; do
    plotsr \
        --sr "${result_dir}/A_digitifera_vs_A_millepora/A_digitifera_vs_A_millepora.syri.out" \
        --sr "${result_dir}/A_millepora_vs_A_tenius/A_millepora_vs_A_tenius.syri.out" \
        --genomes "${three_genomes}" \
        --chrord "${plot_dir}/chromosome_order.txt" \
        -W 20 -H 12 -f 8 -d 600 \
        -o "${plot_dir}/all_three_synteny.${extension}" \
        --lf "${plot_dir}/all_three_synteny.${extension}.plotsr.log"
done

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
