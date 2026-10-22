#!/bin/bash
#SBATCH --job-name=20_assemble
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

# ============================================================
# Step 20 (array, one task per isolate): de novo assembly.
#   trimmed reads -> (downsample to 100x if deeper) -> SPAdes
#   -> keep scaffolds >= 500 bp
#   sbatch --array=1-N 20_assemble.sh
# Output: 07_assembly/<ISOLATE>.scaffolds.fa
# Can run at the same time as steps 10-13.
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_ASM}"

S=$(get_sample)
R1="${TRIM_DIR}/${S}_R1.trim.fq.gz"
R2="${TRIM_DIR}/${S}_R2.trim.fq.gz"
WORK="${ASM_DIR}/${S}"
MAX_COV=100            # >100x adds memory/time, not quality
mkdir -p "${WORK}"
echo "Isolate: ${S}   $(date)"

# 1. downsample very deep isolates
BASES=$(seqkit stats -T -j "${THREADS}" "${R1}" "${R2}" | awk 'NR>1 {s+=$5} END {print s}')
COV=$(( BASES / GENOME_SIZE ))
echo "Coverage after trimming: ${COV}x"
if (( COV > MAX_COV )); then
    echo "Downsampling to ${MAX_COV}x"
    rasusa reads --coverage "${MAX_COV}" --genome-size "${GENOME_SIZE}" --seed 42 \
        -o "${WORK}/sub_R1.fq.gz" -o "${WORK}/sub_R2.fq.gz" "${R1}" "${R2}"
    R1="${WORK}/sub_R1.fq.gz"; R2="${WORK}/sub_R2.fq.gz"
fi

# 2. assemble
MEM_GB=$(( ${SLURM_MEM_PER_NODE:-122880} / 1024 - 8 ))
spades.py --isolate -1 "${R1}" -2 "${R2}" \
    -o "${WORK}/spades" -t "${THREADS}" -m "${MEM_GB}"

# 3. keep >= 500 bp, give contigs short clean names
seqkit seq -m 500 "${WORK}/spades/scaffolds.fasta" \
  | seqkit replace -p '.+' -r "${S}_ctg{nr}" > "${ASM_DIR}/${S}.scaffolds.fa"
seqkit stats -a -T "${ASM_DIR}/${S}.scaffolds.fa"

# 4. free space: keep only logs + final scaffolds
rm -f "${WORK}"/sub_R*.fq.gz
find "${WORK}/spades" -mindepth 1 -maxdepth 1 \
    ! -name 'spades.log' ! -name 'params.txt' ! -name 'scaffolds.fasta' ! -name 'contigs.fasta' \
    -exec rm -rf {} +

echo "Done ${S}   $(date)"
