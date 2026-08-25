#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "${script_dir}/prepare_chromosomes.py"
"${script_dir}/run_pairwise_syri.sh"
"${script_dir}/run_tenius_control.sh"
"${script_dir}/make_synteny_plots.sh"
