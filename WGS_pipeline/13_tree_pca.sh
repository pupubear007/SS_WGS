#!/bin/bash
#SBATCH --job-name=13_tree_pca
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 13 (single job, after step 12):
#   SNP alignment -> IQ-TREE ML phylogeny (1000 ultrafast bootstraps)
#   SNP PCA + pairwise SNP distance matrix
#   sbatch 13_tree_pca.sh
# Key outputs in 06_popgen/:
#   SS_WGS_snps.treefile        tree (open in FigTree / iTOL / ggtree)
#   pca.tsv, pca_variance.tsv   plot PC1 vs PC2 in R
#   snp_distance.tsv            pairwise SNP differences
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_WGS}"
cd "${POP_DIR}"

VCF="${VCF_DIR}/SS_WGS.snps.clean.vcf.gz"
bcftools query -l "${VCF}" > isolates.txt
bcftools query -f '%REF\t%ALT[\t%GT]\n' "${VCF}" > genotypes.tsv

# ---------- VCF -> FASTA alignment, PCA, distance (python) ----------
python3 - <<'PY'
import numpy as np
samples = [l.strip() for l in open("isolates.txt") if l.strip()]
seqs = [[] for _ in samples]
mat = []
for line in open("genotypes.tsv"):
    f = line.rstrip("\n").split("\t")
    ref, alt, gts = f[0], f[1], f[2:]
    row = []
    for i, g in enumerate(gts):
        a = g.replace("|", "/").split("/")[0]       # haploid: "0", "1" or "."
        if a == "0":   seqs[i].append(ref); row.append(0.0)
        elif a == "1": seqs[i].append(alt); row.append(1.0)
        else:          seqs[i].append("N"); row.append(np.nan)
    mat.append(row)

with open("SS_WGS_snps.fasta", "w") as out:
    for s, sq in zip(samples, seqs):
        out.write(f">{s}\n{''.join(sq)}\n")

X = np.array(mat).T                                  # isolates x SNPs
# pairwise SNP distance (sites where both called)
n = len(samples)
D = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        ok = ~np.isnan(X[i]) & ~np.isnan(X[j])
        D[i, j] = int(np.sum(X[i, ok] != X[j, ok]))
with open("snp_distance.tsv", "w") as out:
    out.write("isolate\t" + "\t".join(samples) + "\n")
    for i, s in enumerate(samples):
        out.write(s + "\t" + "\t".join(map(str, D[i])) + "\n")

# PCA: mean-impute missing, center, SVD
col_mean = np.nanmean(X, axis=0)
Xi = np.where(np.isnan(X), col_mean, X) - col_mean
U, S, Vt = np.linalg.svd(Xi, full_matrices=False)
k = min(10, n)
pcs = U[:, :k] * S[:k]
var = S**2 / np.sum(S**2)
with open("pca.tsv", "w") as out:
    out.write("isolate\t" + "\t".join(f"PC{i+1}" for i in range(k)) + "\n")
    for i, s in enumerate(samples):
        out.write(s + "\t" + "\t".join(f"{v:.5f}" for v in pcs[i]) + "\n")
with open("pca_variance.tsv", "w") as out:
    out.write("PC\tpct_variance\n")
    for i in range(k):
        out.write(f"PC{i+1}\t{100*var[i]:.2f}\n")
print(f"{n} isolates, {X.shape[1]} SNPs; PC1 {100*var[0]:.1f}%, PC2 {100*var[1]:.1f}%")
PY

# ---------- ML tree ----------
# +ASC corrects for using variable sites only
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
"${IQTREE}" -s SS_WGS_snps.fasta -m MFP+ASC -B 1000 -T "${THREADS}" \
    --prefix SS_WGS_snps -redo

rm -f genotypes.tsv
echo "Tree : ${POP_DIR}/SS_WGS_snps.treefile"
echo "PCA  : ${POP_DIR}/pca.tsv"
