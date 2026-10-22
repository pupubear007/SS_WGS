#!/bin/bash
#SBATCH --job-name=00b_get_reference
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 00b (once): put the S. sclerotiorum 1980 UF-70 reference in
# /scratch.global/wan00965/reference/Ss_1980_UF70 (REF_STORE in config.sh).
#
# Option A  (same reference as UG100_RNAseq, recommended)
#   JGI MycoCosm needs a login, so download in your browser:
#     https://mycocosm.jgi.doe.gov/Sclsc1  ->  Download
#       Assembly/Genome Assembly (unmasked)/Sclsc1_AssemblyScaffolds.fasta.gz
#       Annotation/Filtered Models ("best")/Genes/Sclsc1_GeneCatalog_genes_20110903.gff.gz
#       Annotation/Filtered Models ("best")/Proteins/Sclsc1_GeneCatalog_proteins_*.aa.fasta.gz  (optional)
#   Upload them (OnDemand "Upload" button) into ${REF_STORE}/jgi_upload/
#   then:   bash 00b_get_reference.sh jgi
#
# Option B  (no login, automatic): NCBI RefSeq GCF_000146945.2, same strain
#   sbatch 00b_get_reference.sh ncbi        (or: bash 00b_get_reference.sh ncbi)
#
# Afterwards config.sh already points at ${REF_STORE}/current/
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh

SOURCE="${1:-jgi}"
mkdir -p "${REF_STORE}/jgi_upload" "${REF_STORE}/current"
cd "${REF_STORE}"
chmod u+w current/* 2>/dev/null || true     # allow re-runs
echo "Reference store: ${REF_STORE}   source: ${SOURCE}   $(date)"

if [[ "${SOURCE}" == "jgi" ]]; then
    G=$(ls jgi_upload/Sclsc1_AssemblyScaffolds.fasta* 2>/dev/null | head -1 || true)
    A=$(ls jgi_upload/Sclsc1_GeneCatalog_genes_*.gff* 2>/dev/null | head -1 || true)
    P=$(ls jgi_upload/Sclsc1_GeneCatalog_proteins_*.fasta* 2>/dev/null | head -1 || true)
    if [[ -z "$G" || -z "$A" ]]; then
        echo "ERROR: put the JGI files in ${REF_STORE}/jgi_upload/ first:"
        echo "  Sclsc1_AssemblyScaffolds.fasta.gz"
        echo "  Sclsc1_GeneCatalog_genes_20110903.gff.gz"
        echo "Found: $(ls jgi_upload/ 2>/dev/null | tr '\n' ' ')"
        exit 1
    fi
    unz() { if [[ "$1" == *.gz ]]; then gunzip -c "$1"; else cat "$1"; fi; }
    unz "$G" > current/Ss_genome.fa
    unz "$A" > current/Ss_genes.gff
    [[ -n "$P" ]] && unz "$P" > current/Ss_proteins.fa
    echo "JGI Sclsc1 (same as UG100_RNAseq)" > current/SOURCE.txt

elif [[ "${SOURCE}" == "ncbi" ]]; then
    ACC="GCF_000146945.2"
    init_conda
    # NCBI datasets CLI in its own small env
    conda env list | grep -q "^ncbi_datasets " || \
        { command -v mamba >/dev/null && mamba create -y -n ncbi_datasets -c conda-forge ncbi-datasets-cli \
          || conda create -y -n ncbi_datasets -c conda-forge ncbi-datasets-cli; }
    set +u; conda activate ncbi_datasets; set -u
    rm -rf ncbi_tmp && mkdir ncbi_tmp && cd ncbi_tmp
    datasets download genome accession "${ACC}" --include genome,gff3,protein --filename ref.zip
    unzip -q -o ref.zip
    D="ncbi_dataset/data/${ACC}"
    cat "${D}"/*_genomic.fna  > ../current/Ss_genome.fa
    cp  "${D}"/genomic.gff       ../current/Ss_genes.gff
    cp  "${D}"/protein.faa       ../current/Ss_proteins.fa
    cd .. && rm -rf ncbi_tmp
    echo "NCBI RefSeq ${ACC} (ASM14694v2)" > current/SOURCE.txt
else
    echo "Usage: bash 00b_get_reference.sh [jgi|ncbi]"; exit 1
fi

# read-only so nothing overwrites it by accident
chmod 444 current/*
echo "===== reference ready ====="
cat current/SOURCE.txt
echo "Genome   : $(grep -c '>' current/Ss_genome.fa) sequences, $(awk '!/^>/{n+=length($0)} END{printf "%.1f Mb", n/1e6}' current/Ss_genome.fa)"
echo "GFF lines: $(grep -vc '^#' current/Ss_genes.gff)"
[[ -f current/Ss_proteins.fa ]] && echo "Proteins : $(grep -c '>' current/Ss_proteins.fa)"
ls -la current/
echo "Next: bash 01_make_sample_list.sh"
