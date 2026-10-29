#!/usr/bin/env bash
# Step 02 - measure methodological properties of each assembly, all with the
# same tools and parameters so the four genomes are directly comparable:
#   * BUSCO cnidaria_odb12 in genome mode (miniprot): completeness and the
#     duplicated-BUSCO rate (retained haplotypes / false duplications inflate it)
#   * BUSCO in protein mode on each annotation (annotation completeness)
#   * Red de novo repeat detection (assemblies were masked inconsistently)
#   * 20 kb fragments mapped back to their own genome: long, near-identical
#     intra-genome copies are candidate false duplications from uncollapsed haplotypes
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
out="${OUT}/02_assembly_qc"
mkdir -p "${out}"/{busco,red,self}
threads="${THREADS:-12}"
species=(Adig Amil Aten Acer)

busco_dl="${OUT}/busco_downloads"
# BUSCO's own downloader is very slow here; fetch the lineage straight from its S3 bucket.
if [[ ! -d "${busco_dl}/lineages/cnidaria_odb12" ]]; then
    mkdir -p "${busco_dl}/lineages"
    curl -sSL https://busco-data.s3.amazonaws.com/lineages/cnidaria_odb12.2026-05-22.tar.gz \
        | tar -xz -C "${busco_dl}/lineages"
fi

for sp in "${species[@]}"; do
    [[ -s "${out}/busco/genome_${sp}/short_summary.specific.cnidaria_odb12.genome_${sp}.json" ]] && continue
    busco -i "${INPUTS}/genomes/${sp}.fa" -m genome --miniprot -l cnidaria_odb12 \
        --offline --download_path "${busco_dl}" -c "${threads}" \
        --out_path "${out}/busco" -o "genome_${sp}" -f > "${out}/busco/genome_${sp}.log" 2>&1
done
for sp in Adig Amil Aten; do
    [[ -s "${out}/busco/protein_${sp}/short_summary.specific.cnidaria_odb12.protein_${sp}.json" ]] && continue
    busco -i "${INPUTS}/annotations/${sp}.proteins.faa" -m proteins -l cnidaria_odb12 \
        --offline --download_path "${busco_dl}" -c "${threads}" \
        --out_path "${out}/busco" -o "protein_${sp}" -f > "${out}/busco/protein_${sp}.log" 2>&1
done

# Red reads every *.fa in a directory; with -frm 2 it writes <name>.bed
# ("chr<TAB>start<TAB>end", zero-based, end-exclusive) into the -rpt directory.
for sp in "${species[@]}"; do
    [[ -s "${out}/red/${sp}.bed" ]] && continue
    mkdir -p "${out}/red/${sp}_in" "${out}/red/${sp}_out"
    ln -sf "${INPUTS}/genomes/${sp}.fa" "${out}/red/${sp}_in/${sp}.fa"
    [[ -s "${out}/red/${sp}_out/${sp}.bed" ]] || Red -gnm "${out}/red/${sp}_in" -rpt "${out}/red/${sp}_out" \
        -frm 2 -cor "${threads}" > "${out}/red/${sp}.log" 2>&1
    cut -f1-3 "${out}/red/${sp}_out/${sp}.bed" | sort -k1,1 -k2,2n | bedtools merge > "${out}/red/${sp}.bed"
done

# False-duplication scan: each genome is cut into 20 kb fragments that are
# mapped back to it. A highly similar secondary hit away from the fragment's own
# origin marks a duplicated segment (filtered in 08_assembly_summary.py).
# -p 0.5 keeps secondary hits down to half the best score: retained haplotypes
# differ by the ~1-3% heterozygosity and score well below a perfect self-hit.
# Much faster than a full self-alignment, whose repeat chains explode with -P.
for sp in "${species[@]}"; do
    [[ -s "${out}/self/${sp}.self.paf.gz" ]] && continue
    seqkit sliding -g -s 20000 -W 20000 "${INPUTS}/genomes/${sp}.fa" 2> /dev/null \
        | minimap2 -x asm5 -c -N 20 -p 0.5 --secondary=yes -t "${threads}" "${INPUTS}/genomes/${sp}.fa" - \
        2> "${out}/self/${sp}.log" | gzip > "${out}/self/${sp}.self.paf.gz.tmp"
    mv "${out}/self/${sp}.self.paf.gz.tmp" "${out}/self/${sp}.self.paf.gz"
done
echo "02 done"
