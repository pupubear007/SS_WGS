#!/bin/bash
#SBATCH --job-name=11_align_call
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

# ============================================================
# Step 11 (array, one task per isolate):
#   bwa-mem2 align -> sort -> mark duplicates -> mapping stats
#   -> GATK HaplotypeCaller (haploid, GVCF mode)
#   sbatch --array=1-N 11_align_call.sh
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_WGS}"

S=$(get_sample)
R1="${TRIM_DIR}/${S}_R1.trim.fq.gz"
R2="${TRIM_DIR}/${S}_R2.trim.fq.gz"
BAM="${ALN_DIR}/${S}.markdup.bam"
TMP="${TMPDIR:-/tmp}/${S}_$$"; mkdir -p "${TMP}"
echo "Isolate: ${S}   $(date)"

# 1. align + duplicate marking in one stream
if [[ ! -s "${BAM}.bai" ]]; then
    bwa-mem2 mem -t "${THREADS}" \
        -R "@RG\tID:${S}\tSM:${S}\tLB:${S}\tPL:ILLUMINA" \
        "${REF}" "${R1}" "${R2}" \
    | samtools collate -@ 4 -O -u -T "${TMP}/collate" - \
    | samtools fixmate -@ 4 -m -u - - \
    | samtools sort    -@ 4 -m 2G -u -T "${TMP}/sort" - \
    | samtools markdup -@ 4 -T "${TMP}/md" -f "${ALN_DIR}/${S}.markdup_stats.txt" - "${BAM}"
    samtools index -@ "${THREADS}" "${BAM}"
fi

# 2. mapping QC
samtools flagstat -@ "${THREADS}" "${BAM}" > "${ALN_DIR}/${S}.flagstat.txt"
samtools coverage "${BAM}"                 > "${ALN_DIR}/${S}.coverage_per_contig.tsv"
awk -v s="${S}" 'NR>1 {L+=$3; cov+=$6*$3; dp+=$7*$3}
    END {printf "%s\tbreadth_pct=%.2f\tmean_depth=%.1f\n", s, cov/L, dp/L}' \
    "${ALN_DIR}/${S}.coverage_per_contig.tsv" | tee "${ALN_DIR}/${S}.coverage_summary.txt"

# 3. per-isolate variant calling (S. sclerotiorum is haploid -> ploidy 1)
gatk --java-options "-Xmx32g -Djava.io.tmpdir=${TMP}" HaplotypeCaller \
    -R "${REF}" -I "${BAM}" \
    -O "${GVCF_DIR}/${S}.g.vcf.gz" \
    -ERC GVCF -ploidy 1 \
    --native-pair-hmm-threads "${THREADS}"

rm -rf "${TMP}"
echo "Done ${S}   $(date)"
