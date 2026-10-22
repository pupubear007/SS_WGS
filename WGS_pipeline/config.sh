#!/bin/bash
# ============================================================
# Shared settings for every step of the SS_WGS pipeline.
# Every script does:  source ./config.sh
# Edit the paths here ONCE; nothing else needs path changes.
# ============================================================

# ---------- SLURM / account ----------
# (the #SBATCH lines in each script use: account mmccaghe, partition msismall)

# ---------- input ----------
RAW_DIR="/scratch.global/wan00965/Quick_test/SS_WGS_rawfastq_merged"   # <ISOLATE>_R1.fq.gz / _R2.fq.gz

# ---------- project output ----------
PROJECT="/scratch.global/wan00965/Quick_test/SS_WGS_analysis"
SAMPLES="${PROJECT}/samples.txt"

# ---------- reference (S. sclerotiorum 1980 UF-70) ----------
# Reference folder, made by 00b_get_reference.sh
REF_STORE="/scratch.global/wan00965/reference/Ss_1980_UF70"
REF_FASTA_SRC="${REF_STORE}/current/Ss_genome.fa"
REF_GFF_SRC="${REF_STORE}/current/Ss_genes.gff"
REF_PROT_SRC="${REF_STORE}/current/Ss_proteins.fa"   # optional; built from genome + GFF in step 10 if missing
REF_DIR="${PROJECT}/ref"
REF="${REF_DIR}/Ss_ref.fasta"
GENOME_SIZE=38900000          # ~38.9 Mb, used for coverage estimates

# ---------- conda environments (made by 00_setup_env.sh) ----------
ENV_WGS="pu_WGS"              # QC, alignment, variant calling, tree
ENV_ASM="pu_WGS_asm"          # SPAdes, QUAST, BUSCO, OrthoFinder
ENV_ANNO="pu_funannotate"     # funannotate
FUNANNOTATE_DB="${PROJECT}/db/funannotate_db"
BUSCO_LINEAGE="leotiomycetes_odb10"

# ---------- output folders ----------
QC_DIR="${PROJECT}/01_qc"
TRIM_DIR="${PROJECT}/02_trimmed"
ALN_DIR="${PROJECT}/03_align"
GVCF_DIR="${PROJECT}/04_gvcf"
VCF_DIR="${PROJECT}/05_vcf"
POP_DIR="${PROJECT}/06_popgen"
ASM_DIR="${PROJECT}/07_assembly"
ASMQC_DIR="${PROJECT}/08_assembly_qc"
ANNO_DIR="${PROJECT}/09_annotation"
ORTHO_DIR="${PROJECT}/10_orthofinder"

# ---------- helpers ----------
# Conda setup:
#   1) your own Miniforge (~/miniforge3, from 00a_install_miniforge.sh): current conda + mamba
#   2) otherwise the UG100_RNAseq way: source ~/.bashrc ; module load conda
MINIFORGE="$HOME/miniforge3"
init_conda() {
    set +u
    if [[ -f "${MINIFORGE}/etc/profile.d/conda.sh" ]]; then
        source "${MINIFORGE}/etc/profile.d/conda.sh"
    else
        [[ -f "$HOME/.bashrc" ]] && source "$HOME/.bashrc"
        module load conda
        source "$(conda info --base)/etc/profile.d/conda.sh"
    fi
    set -u
    command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found (run 00a_install_miniforge.sh)"; exit 1; }
}

activate_env() {
    init_conda
    set +u
    conda activate "$1"
    set -u
    echo "Active env: ${CONDA_DEFAULT_ENV:-none}"
}

# isolate name for this array task (1-based line of samples.txt)
get_sample() {
    sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES}"
}

THREADS="${SLURM_CPUS_PER_TASK:-4}"
