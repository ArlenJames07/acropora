#!/usr/bin/env bash
# Step 04 - gene families and annotation-independent gene loci.
#
#   a) OrthoFinder on the three annotated proteomes (longest isoforms) defines
#      gene families (orthogroups).
#   b) Every annotated protein is aligned with miniprot, using identical settings,
#      to all four genomes (whole assemblies, including unplaced scaffolds).
#      Step 06 resolves these alignments into one locus set per genome. Gene
#      presence, copy number and synteny are then measured the same way in every
#      genome, removing AUGUSTUS-vs-Gnomon annotation bias, and A. cervicornis
#      (unannotated) is included.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
out="${OUT}/04_orthology"
threads="${THREADS:-10}"
mkdir -p "${out}/proteomes" "${out}/miniprot"

for sp in Adig Amil Aten; do
    ln -sf "${INPUTS}/annotations/${sp}.proteins.faa" "${out}/proteomes/${sp}.faa"
done

# DendroBLAST gene trees (-M dendroblast) skip the slow MSA step; only the
# orthogroups are used downstream. If an earlier run finished its DIAMOND
# searches, they are reused with -b.
if [[ -z "$(find "${out}/orthofinder" -name Orthogroups.tsv 2>/dev/null)" ]]; then
    previous=$(compgen -G "${out}/orthofinder/Results_*/WorkingDirectory/Blast2_2.txt.gz" | head -1 || true)
    if [[ -n "${previous}" ]]; then
        orthofinder -b "$(dirname "${previous}")" -M dendroblast -t "${threads}" -a 4 -n og \
            > "${out}/orthofinder.log" 2>&1
    else
        rm -rf "${out}/orthofinder"
        orthofinder -f "${out}/proteomes" -M dendroblast -t "${threads}" -a 4 -o "${out}/orthofinder" -n acropora \
            > "${out}/orthofinder.log" 2>&1
    fi
fi

queries="${out}/miniprot/all_annotated_proteins.faa"
cat "${INPUTS}"/annotations/{Adig,Amil,Aten}.proteins.faa > "${queries}"
for sp in Adig Amil Aten Acer; do
    gff="${out}/miniprot/${sp}.miniprot.gff.gz"
    [[ -s "${gff}" ]] && continue
    [[ -s "${out}/miniprot/${sp}.mpi" ]] || miniprot -t "${threads}" -d "${out}/miniprot/${sp}.mpi" \
        "${INPUTS}/genomes/${sp}.fa" 2> "${out}/miniprot/${sp}.index.log"
    # Secondary hits down to half the best score are kept so paralogous copies
    # are found; overlapping hits are resolved per locus in step 06.
    miniprot -t "${threads}" --gff --trans -N 50 -p 0.5 --outs=0.5 --outc=0.3 \
        "${out}/miniprot/${sp}.mpi" "${queries}" 2> "${out}/miniprot/${sp}.log" \
        | gzip > "${gff}.tmp"
    mv "${gff}.tmp" "${gff}"
    rm -f "${out}/miniprot/${sp}.mpi"
done
echo "04 done"
