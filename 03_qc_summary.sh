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

# seqkit tables (reused if a previous run already made them)
RAW_TSV="${QC_DIR}/seqkit/raw_stats.tsv"
TRIM_TSV="${QC_DIR}/seqkit/trimmed_stats.tsv"
[[ -s "${RAW_TSV}" ]]  || seqkit stats -a -T -j "${THREADS}" "${RAW_DIR}"/*.fq.gz  > "${RAW_TSV}"
[[ -s "${TRIM_TSV}" ]] || seqkit stats -a -T -j "${THREADS}" "${TRIM_DIR}"/*.fq.gz > "${TRIM_TSV}"

# coverage = (R1 + R2 bases) / genome size, raw vs trimmed, one row per isolate
python3 - "${RAW_TSV}" "${TRIM_TSV}" "${GENOME_SIZE}" "${QC_DIR}/seqkit/coverage_per_isolate.tsv" <<'PYEOF'
import csv, os, re, sys
raw_f, trim_f, G, out_f = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
def per_isolate(path, pattern):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            name = re.sub(pattern, "", os.path.basename(r["file"]))
            reads, bases = d.get(name, (0, 0))
            d[name] = (reads + int(r["num_seqs"]), bases + int(r["sum_len"]))
    return d
raw  = per_isolate(raw_f,  r"_R[12]\.fq\.gz$")
trim = per_isolate(trim_f, r"_R[12]\.trim\.fq\.gz$")
with open(out_f, "w") as out:
    out.write("isolate\traw_reads\traw_bases\traw_coverage_x\ttrim_reads\ttrim_bases\ttrim_coverage_x\tpct_bases_kept\n")
    for s in sorted(raw):
        rr, rb = raw[s]; tr, tb = trim.get(s, (0, 0))
        out.write(f"{s}\t{rr}\t{rb}\t{rb/G:.1f}\t{tr}\t{tb}\t{tb/G:.1f}\t{100*tb/rb if rb else 0:.1f}\n")
missing = sorted(set(raw) - set(trim))
if missing:
    print("WARNING: no trimmed reads for:", ", ".join(missing))
print(f"{len(raw)} isolates written to {out_f}")
PYEOF

column -t "${QC_DIR}/seqkit/coverage_per_isolate.tsv"
echo "Lowest trimmed coverage:"
sort -t$'\t' -k7,7g "${QC_DIR}/seqkit/coverage_per_isolate.tsv" | grep -v '^isolate' | head -5 | cut -f1,7

multiqc -f -o "${QC_DIR}" -n multiqc_report.html "${QC_DIR}/fastqc_raw" "${QC_DIR}/fastp"
echo "MultiQC: ${QC_DIR}/multiqc_report.html"
