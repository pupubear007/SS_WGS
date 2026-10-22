#!/bin/bash
#SBATCH --job-name=gather_WGS_fastq
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=gather_WGS_fastq_%j.out
#SBATCH --error=gather_WGS_fastq_%j.err

# ============================================================
# Step 1 (run on MSI): collect all raw WGS FASTQ from the batch
# folders into ONE folder, then write an md5 manifest so the
# local download can be verified.
#
# Usage:
#   bash  gather_raw_fastq.sh --dry-run   # preview only, nothing moved
#   sbatch gather_raw_fastq.sh            # move + md5 (as a SLURM job)
# ============================================================

set -euo pipefail

# ---------- paths (edit if needed) ----------
BASE_DIR="/scratch.global/wan00965/Quick_test"
SRC_DIRS=(
    "${BASE_DIR}/SS_WGS_rawfastq_20260430"
    "${BASE_DIR}/SS_WGS_rawfastq_20260630"
)
OUT_DIR="${BASE_DIR}/SS_WGS_rawfastq_all"
MANIFEST="${OUT_DIR}/md5sum_rawfastq.txt"
LOG="${OUT_DIR}/move_log.tsv"

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

mkdir -p "${OUT_DIR}"

# ---------- 1. find every FASTQ in the batch folders ----------
mapfile -t FASTQS < <(
    find "${SRC_DIRS[@]}" -type f \
        \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) \
    | sort
)

echo "Found ${#FASTQS[@]} FASTQ files in ${#SRC_DIRS[@]} batch folders."
[[ ${#FASTQS[@]} -eq 0 ]] && { echo "ERROR: no FASTQ found. Check SRC_DIRS."; exit 1; }

# ---------- 2. check for duplicate file names across batches ----------
# If the same name exists in both batches, the batch date is prefixed
# so nothing gets overwritten.
DUPES=$(printf '%s\n' "${FASTQS[@]}" | xargs -n1 basename | sort | uniq -d || true)
if [[ -n "${DUPES}" ]]; then
    echo "WARNING: these names appear in more than one batch; batch date will be prefixed:"
    echo "${DUPES}"
fi

# ---------- 3. move into one folder ----------
[[ ${DRY_RUN} -eq 0 ]] && printf "original_path\tnew_path\n" > "${LOG}"

for f in "${FASTQS[@]}"; do
    name=$(basename "${f}")
    if grep -qxF "${name}" <<< "${DUPES}"; then
        batch=$(echo "${f}" | grep -oE 'rawfastq_[0-9]{8}' | cut -d_ -f2)
        name="${batch}_${name}"
    fi
    dest="${OUT_DIR}/${name}"

    if [[ -e "${dest}" ]]; then
        echo "SKIP (already exists): ${dest}"
        continue
    fi

    if [[ ${DRY_RUN} -eq 1 ]]; then
        echo "[dry-run] mv ${f} -> ${dest}"
    else
        mv -n "${f}" "${dest}"          # -n = never overwrite
        printf "%s\t%s\n" "${f}" "${dest}" >> "${LOG}"
    fi
done

[[ ${DRY_RUN} -eq 1 ]] && { echo "Dry run finished. Nothing was moved."; exit 0; }

# ---------- 4. md5 manifest (used to verify the download) ----------
echo "Computing md5 checksums..."
cd "${OUT_DIR}"
find . -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) \
    -printf '%f\n' | sort \
    | xargs -P "${SLURM_CPUS_PER_TASK:-4}" -I{} md5sum {} \
    | sort -k2 > "${MANIFEST}"

N_OUT=$(wc -l < "${MANIFEST}")
echo "Done. ${N_OUT} FASTQ files in ${OUT_DIR}"
du -sh "${OUT_DIR}"
echo "Move log : ${LOG}"
echo "Manifest : ${MANIFEST}"

# Empty batch folders are left in place on purpose; remove them yourself
# after checking:  find "${SRC_DIRS[@]}" -type f | wc -l   (should be 0 FASTQ)
