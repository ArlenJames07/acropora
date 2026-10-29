#!/usr/bin/env bash
# Wait for step 02 (repeat masking + self-alignment) to finish, then rerun every
# step that uses its outputs and rebuild the report. Safe to start at any time:
#   nohup scripts/comparative_genomics/run_downstream.sh > results/comparative_genomics/logs/run_downstream.log 2>&1 &
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/env.sh"

until [[ -s "${OUT}/02_assembly_qc/self/Acer.self.paf.gz" ]] && ! pgrep -f "[0]2_assembly_qc.sh" > /dev/null; do
    sleep 60
done
for step in 07_structural_variants 08_assembly_summary 09_gene_conservation 11_technology_effects 12_report; do
    echo "== ${step} $(date +%H:%M)"
    THREADS="${THREADS:-16}" python "${here}/${step}.py" > "${OUT}/logs/${step}.log" 2>&1
done
echo "downstream done $(date +%H:%M)"
