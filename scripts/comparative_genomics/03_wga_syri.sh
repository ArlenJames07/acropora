#!/usr/bin/env bash
# Step 03 - pairwise whole-genome alignment of the 14 homologous chromosomes and
# SyRI structural annotation for all six species pairs.
#
# minimap2 asm20 is used (not asm5 as in results/syri_results): Acropora species
# differ by more than 5%, and asm5 leaves most of the genome unaligned, which
# SyRI then reports as spurious "not aligned"/translocation calls.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
out="${OUT}/03_wga"
mkdir -p "${out}/chrom"
threads="${THREADS:-10}"
species=(Adig Amil Aten Acer)
chroms=$(printf 'chr%s ' $(seq 1 14))

for sp in "${species[@]}"; do
    fa="${out}/chrom/${sp}.fa"
    [[ -s "${fa}.fai" ]] || { samtools faidx "${INPUTS}/genomes/${sp}.fa" ${chroms} > "${fa}"; samtools faidx "${fa}"; }
done

run_pair() {
    local ref="$1" qry="$2" pair="$1_vs_$2"
    local dir="${out}/${pair}"
    mkdir -p "${dir}"
    [[ -s "${dir}/${pair}.syri.out" ]] && { echo "skip ${pair}"; return; }
    if [[ ! -s "${dir}/${pair}.paf" ]]; then
        minimap2 -c -x asm20 --eqx -t "${threads}" "${out}/chrom/${ref}.fa" "${out}/chrom/${qry}.fa" \
            > "${dir}/${pair}.paf.tmp" 2> "${dir}/minimap2.log"
        mv "${dir}/${pair}.paf.tmp" "${dir}/${pair}.paf"
    fi
    # Only alignments between homologous chromosomes are given to SyRI; other
    # hits are kept in the PAF for the inter-chromosomal analysis in step 06.
    syri -c "${dir}/${pair}.paf" -F P -r "${out}/chrom/${ref}.fa" -q "${out}/chrom/${qry}.fa" \
        --nc 4 --dir "${dir}" --prefix "${pair}." --samplename "${qry}" \
        --lf syri.log > "${dir}/syri.stdout" 2>&1
}

pairs=("Adig Amil" "Adig Aten" "Adig Acer" "Amil Aten" "Amil Acer" "Aten Acer")
# Three pairs at a time keeps memory well under the machine limit.
for pair in "${pairs[@]}"; do
    run_pair ${pair} &
    while (( $(jobs -rp | wc -l) >= 3 )); do wait -n; done
done
wait
echo "03 done"
