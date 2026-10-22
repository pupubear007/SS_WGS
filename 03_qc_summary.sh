#!/bin/bash
#SBATCH --job-name=03_qc_summary
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 03 (single job, after ALL of step 02 finished):
# read counts, bases, coverage per isolate + MultiQC report.
#   sbatch 03_qc_summary.sh
# Key outputs:
#   01_qc/seqkit/coverage_per_isolate.tsv
#   01_qc/multiqc_report.html
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_WGS}"

seqkit stats -a -T -j "${THREADS}" "${RAW_DIR}"/*.fq.gz   > "${QC_DIR}/seqkit/raw_stats.tsv"
seqkit stats -a -T -j "${THREADS}" "${TRIM_DIR}"/*.fq.gz  > "${QC_DIR}/seqkit/trimmed_stats.tsv"

# coverage = (R1 + R2 bases) / genome size
summarize() {   # $1 = stats file, $2 = suffix to strip
    awk -F'\t' -v G="${GENOME_SIZE}" -v suf="$2" 'NR>1 {
        n=split($1,p,"/"); f=p[n]; sub(suf,"",f)
        reads[f]+=$4; bases[f]+=$5 }
        END { for (s in reads) printf "%s\t%d\t%.0f\t%.1f\n", s, reads[s], bases[s], bases[s]/G }' "$1" | sort
}
{
    printf "isolate\traw_reads\traw_bases\traw_coverage_x\ttrim_reads\ttrim_bases\ttrim_coverage_x\tpct_bases_kept\n"
    join -t $'\t' \
        <(summarize "${QC_DIR}/seqkit/raw_stats.tsv"     '_R[12]\\.fq\\.gz$') \
        <(summarize "${QC_DIR}/seqkit/trimmed_stats.tsv" '_R[12]\\.trim\\.fq\\.gz$') \
    | awk -F'\t' -v OFS='\t' '{print $0, sprintf("%.1f", 100*$6/$3)}'
} > "${QC_DIR}/seqkit/coverage_per_isolate.tsv"

column -t "${QC_DIR}/seqkit/coverage_per_isolate.tsv"

multiqc -f -o "${QC_DIR}" -n multiqc_report.html "${QC_DIR}/fastqc_raw" "${QC_DIR}/fastp"
echo "MultiQC: ${QC_DIR}/multiqc_report.html"
