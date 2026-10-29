#!/usr/bin/env bash
# Run the whole comparative-genomics workflow. Every step skips work whose
# outputs already exist, so the script can be re-run after an interruption.
#   THREADS=24 ./scripts/comparative_genomics/run_all.sh
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${here}/00_setup_env.sh"
source "${here}/env.sh"
mkdir -p "${OUT}/logs"
threads="${THREADS:-24}"

[[ -s "${INPUTS}/assembly_stats.tsv" ]] || python "${here}/01_prepare_inputs.py"
# 02-04 are independent; run them side by side.
THREADS=$((threads / 3)) "${here}/02_assembly_qc.sh"          > "${OUT}/logs/02_assembly_qc.log" 2>&1 &
THREADS=$((threads / 3)) "${here}/03_wga_syri.sh"             > "${OUT}/logs/03_wga_syri.log" 2>&1 &
THREADS=$((threads / 3)) "${here}/04_orthology_projection.sh" > "${OUT}/logs/04_orthology.log" 2>&1 &
wait
for step in 05_gene_loci 06_gene_synteny 08_assembly_summary 07_structural_variants \
            09_gene_conservation 10_sequence_evolution 11_technology_effects 12_report; do
    echo "== ${step}"
    THREADS="${threads}" python "${here}/${step}.py" > "${OUT}/logs/${step}.log" 2>&1
done
