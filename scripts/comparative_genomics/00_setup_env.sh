#!/usr/bin/env bash
# Create the isolated conda environment used by every comparative_genomics script.
# PYTHONNOUSERSITE=1 is exported by the run scripts so ~/.local packages cannot leak in.
set -euo pipefail
ENV_NAME="${ENV_NAME:-acropora_cg}"
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Environment ${ENV_NAME} already exists"; exit 0
fi
mamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
    python=3.11 "busco>=6" miniprot minimap2 samtools seqkit bedtools diamond \
    orthofinder gffread mafft iqtree paml mcscanx syri plotsr \
    red pandas numpy scipy matplotlib
