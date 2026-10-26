#!/bin/bash
# ============================================================
# Step 3 (run on your Mac, NOT on MSI): download the merged
# WGS FASTQ (one R1 + one R2 per isolate) from MSI to the
# pu_R1_4TB SSD and verify md5.
#
# Usage:
#   bash download_raw_fastq.sh
#
# Safe to re-run: rsync resumes partial files and skips finished ones.
# All output is saved to SS_WGS/log/ on the SSD.
# ============================================================

set -euo pipefail

# ---------- settings ----------
MSI_USER="wan00965"
MSI_HOST="agate.msi.umn.edu"                       # MSI login host
REMOTE_DIR="/scratch.global/wan00965/Quick_test/SS_WGS_rawfastq_merged"
LOCAL_SSD="/Volumes/pu_R1_4TB"
PROJECT_DIR="${LOCAL_SSD}/SS_WGS"
LOCAL_DIR="${PROJECT_DIR}/SS_WGS_rawfastq_merged"
LOG_DIR="${PROJECT_DIR}/log"
MANIFEST_NAME="md5sum_merged.txt"

# ---------- logging: everything printed also goes to log/ ----------
if [[ ! -d "${LOCAL_SSD}" ]]; then
    echo "ERROR: SSD not mounted at ${LOCAL_SSD}. Plug it in and check: ls /Volumes"
    exit 1
fi
mkdir -p "${LOCAL_DIR}" "${LOG_DIR}"
LOG="${LOG_DIR}/download_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1
echo "Log file: ${LOG}"
echo "Start   : $(date)"

# ---------- 1. space check ----------
REMOTE_SIZE=$(ssh "${MSI_USER}@${MSI_HOST}" "du -sk '${REMOTE_DIR}' | cut -f1")
FREE_SIZE=$(df -k "${LOCAL_SSD}" | awk 'NR==2 {print $4}')
echo "Remote data : $((REMOTE_SIZE / 1024 / 1024)) GB"
echo "SSD free    : $((FREE_SIZE  / 1024 / 1024)) GB"
if (( REMOTE_SIZE > FREE_SIZE )); then
    echo "ERROR: not enough space on SSD."
    exit 1
fi

# ---------- 2. download (resumable) ----------
# -a archive, -h human sizes, --partial keep half-done files,
# --info=progress2 one overall progress line. FASTQ.gz is already
# compressed, so no -z.
rsync -ah --partial --info=progress2 --stats \
    "${MSI_USER}@${MSI_HOST}:${REMOTE_DIR}/" \
    "${LOCAL_DIR}/" \
  || rsync -avh --partial --progress \
    "${MSI_USER}@${MSI_HOST}:${REMOTE_DIR}/" \
    "${LOCAL_DIR}/"          # fallback for the older rsync built into macOS

# ---------- 3. verify md5 against the MSI manifest ----------
MANIFEST="${LOCAL_DIR}/${MANIFEST_NAME}"
if [[ ! -f "${MANIFEST}" ]]; then
    echo "WARNING: ${MANIFEST_NAME} not found; skipping verification."
    exit 0
fi

echo "Verifying md5 checksums..."
cd "${LOCAL_DIR}"
FAIL=0
N=0
while read -r expected file; do
    N=$((N + 1))
    actual=$(md5 -q "${file}")                     # macOS md5
    if [[ "${actual}" == "${expected}" ]]; then
        echo "OK    ${file}"
    else
        echo "FAIL  ${file}  (removing so the next run re-downloads it)"
        rm -f "${file}"
        FAIL=$((FAIL + 1))
    fi
done < "${MANIFEST}"

echo "Checked ${N} files, ${FAIL} failed."
echo "Local size: $(du -sh "${LOCAL_DIR}" | cut -f1)"
echo "End     : $(date)"
if (( FAIL == 0 )); then
    echo "All files verified. Data is in ${LOCAL_DIR}"
else
    echo "Re-run this script to re-download the failed files."
    exit 1
fi
