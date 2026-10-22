#!/bin/bash
#SBATCH --job-name=10_ref_prep
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 10 (single job, once): index the reference for
# bwa-mem2 / samtools / GATK and build a SnpEff database.
#   sbatch 10_ref_prep.sh
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_WGS}"

cd "${REF_DIR}"
bwa-mem2 index Ss_ref.fasta
samtools faidx Ss_ref.fasta
rm -f Ss_ref.dict
gatk CreateSequenceDictionary -R Ss_ref.fasta -O Ss_ref.dict

# ---------- gene models: JGI 2011 GFF -> GTF + protein FASTA ----------
command -v gffread >/dev/null 2>&1 || mamba install -y -n "${ENV_WGS}" -c conda-forge -c bioconda gffread \
  || conda install -y -n "${ENV_WGS}" -c conda-forge -c bioconda gffread
gffread Ss_ref.gff -T -o Ss_ref.gtf
echo "GTF transcripts: $(cut -f9 Ss_ref.gtf | grep -o 'transcript_id "[^"]*"' | sort -u | wc -l)"
if [[ ! -s Ss_ref_proteins.fa ]]; then
    gffread Ss_ref.gff -g Ss_ref.fasta -y Ss_ref_proteins.fa
    sed -i 's/\.$//; /^[^>]/ s/\.//g' Ss_ref_proteins.fa      # drop stop-codon dots
fi
echo "Reference proteins: $(grep -c '>' Ss_ref_proteins.fa)"

# ---------- SnpEff database "Sclsc" ----------
SNPEFF_DIR="${REF_DIR}/snpeff"
mkdir -p "${SNPEFF_DIR}/data/Sclsc"
cp Ss_ref.fasta "${SNPEFF_DIR}/data/Sclsc/sequences.fa"
cp Ss_ref.gtf   "${SNPEFF_DIR}/data/Sclsc/genes.gtf"
cat > "${SNPEFF_DIR}/snpEff.config" <<EOF
data.dir = ${SNPEFF_DIR}/data/
Sclsc.genome : Sclerotinia_sclerotiorum_1980_UF-70
EOF
snpEff build -c "${SNPEFF_DIR}/snpEff.config" -gtf22 -noCheckCds -noCheckProtein -v Sclsc \
  || echo "WARNING: SnpEff build failed; variant calling still works, only gene-effect annotation is skipped."

grep -c ">" Ss_ref.fasta | xargs echo "Reference contigs:"
echo "Reference ready in ${REF_DIR}"
