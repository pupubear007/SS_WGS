#!/bin/bash
#SBATCH --job-name=00_setup_env
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 00 (run ONCE): build conda environments + databases.
#   sbatch 00_setup_env.sh
# Three separate envs avoid dependency conflicts.
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
init_conda
echo "Using conda: $(command -v conda)"

CH="-c conda-forge -c bioconda"

# mamba solves large bioconda envs in minutes; the old conda solver can take hours
# (MSI module's mamba 0.5.1 is too old to download today's channels, so only
#  use mamba from your own Miniforge; otherwise fall back to conda)
if [[ -x "${MINIFORGE}/bin/mamba" ]]; then CREATE="mamba create"; else CREATE="conda create"; fi
echo "Solver: ${CREATE}"

# 1. QC + reference-based
conda env list | grep -q "^${ENV_WGS} " || \
${CREATE} -y -n "${ENV_WGS}" ${CH} \
    fastqc multiqc fastp seqkit pigz \
    bwa-mem2 samtools bcftools gatk4 snpeff \
    iqtree python=3.11 numpy

# 2. de novo assembly + assembly QC + orthology
conda env list | grep -q "^${ENV_ASM} " || \
${CREATE} -y -n "${ENV_ASM}" ${CH} \
    spades rasusa seqkit quast busco orthofinder

# 3. annotation
conda env list | grep -q "^${ENV_ANNO} " || \
${CREATE} -y -n "${ENV_ANNO}" ${CH} funannotate

# funannotate databases (UniProt, Pfam, BUSCO dikarya, ...)
mkdir -p "${FUNANNOTATE_DB}"
set +u; conda activate "${ENV_ANNO}"; set -u
funannotate setup -i all -d "${FUNANNOTATE_DB}"
set +u; conda deactivate; set -u

echo "Environments ready:"
conda env list | grep -E "${ENV_WGS}|${ENV_ASM}|${ENV_ANNO}"
