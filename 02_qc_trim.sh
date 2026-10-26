#!/bin/bash
#SBATCH --job-name=02_qc_trim
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

# ============================================================
# Step 02 (array, one task per isolate): FastQC raw + fastp trim.
#   sbatch --array=1-N 02_qc_trim.sh
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_WGS}"

S=$(get_sample)
R1="${RAW_DIR}/${S}_R1.fq.gz"
R2="${RAW_DIR}/${S}_R2.fq.gz"
echo "Isolate: ${S}   $(date)"

# 1. raw read quality
fastqc -t "${THREADS}" --outdir "${QC_DIR}/fastqc_raw" "${R1}" "${R2}"

# 2. adapter + quality trimming
#    auto-detect adapters, trim 3' ends below Q20, drop reads < 50 bp
fastp \
    -i "${R1}" -I "${R2}" \
    -o "${TRIM_DIR}/${S}_R1.trim.fq.gz" -O "${TRIM_DIR}/${S}_R2.trim.fq.gz" \
    --detect_adapter_for_pe \
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
    --length_required 50 \
    --thread "${THREADS}" \
    --html "${QC_DIR}/fastp/${S}.fastp.html" \
    --json "${QC_DIR}/fastp/${S}.fastp.json"

echo "Done ${S}   $(date)"
