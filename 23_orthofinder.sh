#!/bin/bash
#SBATCH --job-name=23_orthofinder
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 23 (single job, after ALL of step 22 finished):
#   OrthoFinder on all isolate proteomes + reference proteome
#   -> core / accessory / isolate-specific gene families
#   sbatch 23_orthofinder.sh
# Key outputs in 10_orthofinder/proteomes/OrthoFinder/Results_*/:
#   Orthogroups/Orthogroups.GeneCount.tsv     presence/absence matrix
#   Orthogroups/Orthogroups_UnassignedGenes.tsv  isolate-specific genes
#   Species_Tree/SpeciesTree_rooted.txt
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_ASM}"

PROT="${ORTHO_DIR}/proteomes"
mkdir -p "${PROT}"
while read -r s; do
    f=$(ls "${ANNO_DIR}/${s}/fun/annotate_results/"*.proteins.fa 2>/dev/null | head -1 || true)
    if [[ -n "${f}" ]]; then cp "${f}" "${PROT}/${s}.fa"
    else echo "WARNING: no proteins for ${s} (step 22 not finished?)"; fi
done < "${SAMPLES}"
cp "${REF_DIR}/Ss_ref_proteins.fa" "${PROT}/Ss_1980_UF70_ref.fa"

ls "${PROT}"
orthofinder -f "${PROT}" -t "${THREADS}" -a 8

RES=$(ls -d "${PROT}"/OrthoFinder/Results_* | tail -1)
# core = in every proteome; accessory = in some
awk -F'\t' 'NR==1 {n=NF-2; next}
    { present=0; for (i=2;i<=n+1;i++) if ($i>0) present++
      if (present==n) core++; else if (present>1) acc++; else uniq++ }
    END { printf "Orthogroups: core=%d  accessory=%d  single-proteome=%d\n", core, acc, uniq }' \
    "${RES}/Orthogroups/Orthogroups.GeneCount.tsv" | tee "${ORTHO_DIR}/pangenome_summary.txt"
echo "Results: ${RES}"
