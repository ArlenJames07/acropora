# Sourced by every shell step: isolated tool environment and shared paths.
export PYTHONNOUSERSITE=1
export PATH="/home/rare/miniconda3/envs/acropora_cg/bin:${PATH}"
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="$(cd "${SCRIPTS}/../.." && pwd)"
OUT="${PROJECT}/results/comparative_genomics"
INPUTS="${OUT}/01_inputs"
export SCRIPTS PROJECT OUT INPUTS
