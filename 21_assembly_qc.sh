#!/bin/bash
#SBATCH --job-name=21_assembly_qc
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 21 (single job, after ALL of step 20 finished):
#   QUAST (size, N50, contigs, vs reference) + BUSCO completeness
#   sbatch 21_assembly_qc.sh
# Key outputs in 08_assembly_qc/:
#   quast/report.tsv                    one column per isolate
#   busco_summary.tsv                   complete / fragmented / missing %
# Good assembly: ~38-40 Mb total, BUSCO complete > 95%
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_ASM}"
cd "${ASMQC_DIR}"

ASMS=( "${ASM_DIR}"/*.scaffolds.fa )
echo "Assemblies: ${#ASMS[@]}"

# 1. QUAST against the reference
#    QUAST needs Python <=3.11 (it imports distutils), so it has its own env
QUAST="$(conda info --base)/envs/pu_quast/bin/quast.py"
[[ -x "${QUAST}" ]] || QUAST="quast.py"
"${QUAST}" "${ASMS[@]}" -r "${REF}" --fungus --large -t "${THREADS}" -o quast \
    -l "$(printf '%s\n' "${ASMS[@]}" | xargs -n1 basename | sed 's/\.scaffolds\.fa//' | paste -sd,)"

# 2. BUSCO (reference included as a benchmark)
mkdir -p busco busco_downloads busco_short
for fa in "${ASMS[@]}" "${REF}"; do
    name=$(basename "${fa}" | sed 's/\.scaffolds\.fa$//; s/\.fasta$//')
    [[ -d "busco/${name}" ]] && { echo "BUSCO done: ${name}"; continue; }
    busco -i "${fa}" -m genome -l "${BUSCO_LINEAGE}" -c "${THREADS}" \
        -o "${name}" --out_path busco --download_path busco_downloads
    cp busco/"${name}"/short_summary*.txt busco_short/ 2>/dev/null || true
done

# 3. one-table BUSCO summary
{
  printf "assembly\tcomplete_pct\tsingle_pct\tduplicated_pct\tfragmented_pct\tmissing_pct\n"
  for d in busco/*/; do
      name=$(basename "$d")
      line=$(grep -h "C:" "$d"/short_summary*.txt | head -1)
      echo "${line}" | sed -E "s/.*C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%.*/${name}\t\1\t\2\t\3\t\4\t\5/"
  done
} > busco_summary.tsv
column -t busco_summary.tsv

busco --plot busco_short 2>/dev/null || generate_plot.py -wd busco_short 2>/dev/null || true
echo "QUAST: ${ASMQC_DIR}/quast/report.html"
