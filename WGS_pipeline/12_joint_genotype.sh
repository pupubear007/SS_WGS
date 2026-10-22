#!/bin/bash
#SBATCH --job-name=12_joint_genotype
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 12 (single job, after ALL of step 11 finished):
#   combine GVCFs -> joint genotyping -> GATK hard filters
#   -> clean SNP set for tree/PCA -> SnpEff annotation
#   sbatch 12_joint_genotype.sh
# Key outputs in 05_vcf/:
#   SS_WGS.filtered.vcf.gz         all variants, FILTER column set
#   SS_WGS.PASS.snpeff.vcf.gz      PASS SNPs + indels with gene effects
#   SS_WGS.snps.clean.vcf.gz       biallelic SNPs for phylogeny / PCA
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_WGS}"
cd "${VCF_DIR}"
JAVA="-Xmx40g -Djava.io.tmpdir=${TMPDIR:-/tmp}"

# 1. combine + joint genotype
V_ARGS=()
while read -r s; do V_ARGS+=(-V "${GVCF_DIR}/${s}.g.vcf.gz"); done < "${SAMPLES}"
gatk --java-options "${JAVA}" CombineGVCFs  -R "${REF}" "${V_ARGS[@]}" -O SS_WGS.combined.g.vcf.gz
gatk --java-options "${JAVA}" GenotypeGVCFs -R "${REF}" -V SS_WGS.combined.g.vcf.gz -O SS_WGS.raw.vcf.gz

# 2. GATK recommended hard filters (SNPs and indels separately)
gatk SelectVariants -V SS_WGS.raw.vcf.gz -select-type SNP   -O raw.snps.vcf.gz
gatk SelectVariants -V SS_WGS.raw.vcf.gz -select-type INDEL -O raw.indels.vcf.gz

gatk VariantFiltration -V raw.snps.vcf.gz -O filt.snps.vcf.gz \
    -filter "QD < 2.0"              --filter-name "QD2" \
    -filter "FS > 60.0"             --filter-name "FS60" \
    -filter "SOR > 3.0"             --filter-name "SOR3" \
    -filter "MQ < 40.0"             --filter-name "MQ40" \
    -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

gatk VariantFiltration -V raw.indels.vcf.gz -O filt.indels.vcf.gz \
    -filter "QD < 2.0"               --filter-name "QD2" \
    -filter "FS > 200.0"             --filter-name "FS200" \
    -filter "SOR > 10.0"             --filter-name "SOR10" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

gatk MergeVcfs -I filt.snps.vcf.gz -I filt.indels.vcf.gz -O SS_WGS.filtered.vcf.gz

# 3. PASS variants + SnpEff gene effects
bcftools view -f PASS -Oz -o SS_WGS.PASS.vcf.gz SS_WGS.filtered.vcf.gz
bcftools index -t SS_WGS.PASS.vcf.gz
SNPEFF_CFG="${REF_DIR}/snpeff/snpEff.config"
if [[ -d "${REF_DIR}/snpeff/data/Sclsc" ]] && ls "${REF_DIR}/snpeff/data/Sclsc"/snpEffectPredictor.bin >/dev/null 2>&1; then
    snpEff ann -c "${SNPEFF_CFG}" -stats SS_WGS.snpeff_summary.html Sclsc SS_WGS.PASS.vcf.gz \
        | bgzip > SS_WGS.PASS.snpeff.vcf.gz
    bcftools index -t SS_WGS.PASS.snpeff.vcf.gz
    # table of moderate/high-impact changes, per gene
    bcftools +split-vep SS_WGS.PASS.snpeff.vcf.gz -a ANN \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%Annotation\t%Annotation_Impact\t%Gene_Name\t%HGVS.p\t[%GT,]\n' \
        -i 'Annotation_Impact="HIGH" || Annotation_Impact="MODERATE"' -d \
        > SS_WGS.high_moderate_effects.tsv 2>/dev/null \
      || echo "NOTE: split-vep table skipped; use SS_WGS.PASS.snpeff.vcf.gz directly"
else
    echo "WARNING: SnpEff database not built (see 10_ref_prep log); skipping annotation."
fi

# 4. clean SNP set for phylogeny / PCA
#    - PASS biallelic SNPs
#    - genotypes with depth < 5 set to missing
#    - sites missing in > 10% of isolates removed
#    - both alleles must still be present
bcftools view -f PASS -m2 -M2 -v snps SS_WGS.filtered.vcf.gz \
  | bcftools +setGT -- -t q -n . -i 'FMT/DP<5' \
  | bcftools view -i 'F_MISSING<0.1' \
  | bcftools view --min-ac 1:minor -Oz -o SS_WGS.snps.clean.vcf.gz
bcftools index -t SS_WGS.snps.clean.vcf.gz

rm -f raw.snps.vcf.gz* raw.indels.vcf.gz* filt.snps.vcf.gz* filt.indels.vcf.gz*

echo "===== variant counts ====="
echo "raw           : $(bcftools view -H SS_WGS.raw.vcf.gz      | wc -l)"
echo "PASS          : $(bcftools view -H SS_WGS.PASS.vcf.gz     | wc -l)"
echo "clean SNPs    : $(bcftools view -H SS_WGS.snps.clean.vcf.gz | wc -l)"
bcftools stats -s - SS_WGS.PASS.vcf.gz > SS_WGS.PASS.bcftools_stats.txt
echo "Per-isolate stats: ${VCF_DIR}/SS_WGS.PASS.bcftools_stats.txt (PSC lines)"
