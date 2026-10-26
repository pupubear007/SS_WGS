#!/bin/bash
#SBATCH --job-name=merge_WGS_fastq
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 2 (run on MSI): size report + merge lanes per isolate
# + rename to clean names.
#
#   BN166_CKDN260011204-1A_23K7VGLT4_L5_1.fq.gz  \
#   BN166_CKDN260011204-1A_23K7VKLT4_L6_1.fq.gz  --> BN166_R1.fq.gz
#   (same for _2 --> BN166_R2.fq.gz)
#
# Gzip files can be joined with cat; the result is a valid .fq.gz.
# R1 and R2 are joined in the SAME lane order, so read pairs stay matched.
# Originals are NOT deleted.
#
# Usage (from the repo root, so log/ is found):
#   bash   raw_file_manage/merge_rename_fastq.sh --size      # size report only
#   bash   raw_file_manage/merge_rename_fastq.sh --dry-run   # show the merge plan only
#   sbatch raw_file_manage/merge_rename_fastq.sh             # run the merge
# ============================================================

set -euo pipefail

# ---------- paths ----------
BASE_DIR="/scratch.global/wan00965/Quick_test"
IN_DIR="${BASE_DIR}/SS_WGS_rawfastq_all"
OUT_DIR="${BASE_DIR}/SS_WGS_rawfastq_merged"
LOG_DIR="${BASE_DIR}/script/log"
mkdir -p "${LOG_DIR}"

MODE="${1:-run}"
TS=$(date +%Y%m%d_%H%M%S)
PLAN="${LOG_DIR}/merge_plan_${TS}.tsv"

# sample name = everything before "_CKDN" (keeps names like BobH_NEB)
sample_of() { basename "$1" | sed 's/_CKDN.*//'; }

# ---------- 1. size report ----------
echo "===== Size report: ${IN_DIR} ====="
echo "Total size        : $(du -sh "${IN_DIR}" | cut -f1)"
echo "Number of fq.gz   : $(ls "${IN_DIR}"/*.fq.gz | wc -l)"
echo
printf "%-15s %8s %10s\n" "sample" "n_files" "size_GB"
for f in "${IN_DIR}"/*.fq.gz; do
    printf "%s\t%s\n" "$(sample_of "$f")" "$(stat -c %s "$f")"
done | awk -F'\t' '{n[$1]++; s[$1]+=$2}
    END {for (k in n) printf "%-15s %8d %10.2f\n", k, n[k], s[k]/1e9}' | sort
echo
[[ "${MODE}" == "--size" ]] && exit 0

# ---------- 2. build merge plan ----------
mapfile -t SAMPLES < <(for f in "${IN_DIR}"/*_1.fq.gz; do sample_of "$f"; done | sort -u)
echo "Isolates found: ${#SAMPLES[@]}"

printf "sample\tread\tn_lanes\tinput_files\n" > "${PLAN}"
for s in "${SAMPLES[@]}"; do
    mapfile -t R1 < <(ls "${IN_DIR}/${s}"_CKDN*_1.fq.gz | sort)
    R2=()
    for r1 in "${R1[@]}"; do
        r2="${r1%_1.fq.gz}_2.fq.gz"
        [[ -f "$r2" ]] || { echo "ERROR: missing mate for ${r1}"; exit 1; }
        R2+=("$r2")
    done
    printf "%s\tR1\t%d\t%s\n" "$s" "${#R1[@]}" "$(printf '%s,' "${R1[@]##*/}")" >> "${PLAN}"
    printf "%s\tR2\t%d\t%s\n" "$s" "${#R2[@]}" "$(printf '%s,' "${R2[@]##*/}")" >> "${PLAN}"
done
echo "Merge plan written: ${PLAN}"
column -t -s $'\t' "${PLAN}" | cut -c1-200

[[ "${MODE}" == "--dry-run" ]] && { echo "Dry run: nothing merged."; exit 0; }

# ---------- 3. merge + rename ----------
mkdir -p "${OUT_DIR}"
CHECK="${LOG_DIR}/merge_sizecheck_${TS}.tsv"
printf "file\texpected_bytes\tactual_bytes\tstatus\n" > "${CHECK}"

merge_one() {   # $1 = output, rest = inputs
    local out="$1"; shift
    if [[ -s "${out}" ]]; then echo "SKIP (exists): ${out}"; return; fi
    cat "$@" > "${out}.tmp" && mv "${out}.tmp" "${out}"
    local exp act
    exp=$(stat -c %s "$@" | awk '{s+=$1} END {print s}')
    act=$(stat -c %s "${out}")
    local st="OK"; [[ "$exp" == "$act" ]] || st="SIZE_MISMATCH"
    printf "%s\t%s\t%s\t%s\n" "${out##*/}" "$exp" "$act" "$st" >> "${CHECK}"
    echo "${st}  ${out##*/}  ($# lanes)"
}

for s in "${SAMPLES[@]}"; do
    mapfile -t R1 < <(ls "${IN_DIR}/${s}"_CKDN*_1.fq.gz | sort)
    R2=("${R1[@]/%_1.fq.gz/_2.fq.gz}")
    merge_one "${OUT_DIR}/${s}_R1.fq.gz" "${R1[@]}"
    merge_one "${OUT_DIR}/${s}_R2.fq.gz" "${R2[@]}"
done

# ---------- 4. integrity test + md5 ----------
echo "Testing gzip integrity..."
ls "${OUT_DIR}"/*.fq.gz | xargs -P "${SLURM_CPUS_PER_TASK:-4}" -I{} \
    sh -c 'gzip -t "{}" && echo "gzip OK  $(basename {})" || echo "gzip FAIL  $(basename {})"'

echo "Computing md5..."
cd "${OUT_DIR}"
ls *.fq.gz | xargs -P "${SLURM_CPUS_PER_TASK:-4}" -I{} md5sum {} | sort -k2 > md5sum_merged.txt

echo
echo "===== Done ====="
echo "Merged files : $(ls "${OUT_DIR}"/*.fq.gz | wc -l)  (expected $(( ${#SAMPLES[@]} * 2 )))"
echo "Merged size  : $(du -sh "${OUT_DIR}" | cut -f1)"
echo "Size check   : ${CHECK}  ($(grep -c SIZE_MISMATCH "${CHECK}" || true) mismatches)"
echo "md5          : ${OUT_DIR}/md5sum_merged.txt"
echo "Originals in ${IN_DIR} were kept. Delete them only after you have checked the merge."
