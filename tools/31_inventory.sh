#!/bin/bash
#SBATCH --job-name=31_inventory
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 31 (single job): list every file the pipeline produced.
# Slow on the login node because funannotate / OrthoFinder / SPAdes
# work folders hold hundreds of thousands of small files.
#
#   sbatch 31_inventory.sh
#
# Outputs (SS_WGS_analysis/inventory/, summaries also copied into
# SS_WGS_results/ so they download with the results):
#   inventory.tsv          every file: size_bytes, date, path
#   by_folder.tsv          per top-level folder: files, total GB
#   by_folder_type.tsv     per folder x file type: files, total MB
#   largest_files.tsv      the 50 largest files
#   busiest_dirs.tsv       folders with the most files (work folders)
# The summary tables are also printed in log/31_inventory_<jobid>.out
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh

INV="${PROJECT}/inventory"
mkdir -p "${INV}"
cd "${PROJECT}"
echo "Inventory of ${PROJECT}   $(date)"

# 1. full listing (skip the inventory folder itself)
find . -path ./inventory -prune -o -type f -printf '%s\t%TY-%Tm-%Td\t%p\n' \
    | sort -t$'\t' -k3,3 > "${INV}/inventory.tsv"
echo "Files listed: $(wc -l < "${INV}/inventory.tsv")"

# 2-5. summaries
python3 - "${INV}" <<'PY'
import os, sys, collections
inv = sys.argv[1]
by_dir = collections.defaultdict(lambda: [0, 0])
by_type = collections.defaultdict(lambda: [0, 0])
by_leaf = collections.Counter()
rows = []
for line in open(os.path.join(inv, "inventory.tsv")):
    size, date, path = line.rstrip("\n").split("\t", 2)
    size = int(size)
    parts = path.lstrip("./").split("/")
    top = parts[0] if len(parts) > 1 else "(top level)"
    name = parts[-1]
    # file type = last extension, keeping .gz pairs (e.g. fq.gz, vcf.gz)
    bits = name.split(".")
    ext = "none" if len(bits) == 1 else (bits[-2] + ".gz" if bits[-1] == "gz" and len(bits) > 2 else bits[-1])
    by_dir[top][0] += 1;            by_dir[top][1] += size
    by_type[(top, ext)][0] += 1;    by_type[(top, ext)][1] += size
    by_leaf["/".join(parts[:-1]) or "."] += 1
    rows.append((size, path))

with open(os.path.join(inv, "by_folder.tsv"), "w") as o:
    o.write("folder\tfiles\tsize_GB\n")
    for k, (n, s) in sorted(by_dir.items()):
        o.write(f"{k}\t{n}\t{s/1e9:.2f}\n")
with open(os.path.join(inv, "by_folder_type.tsv"), "w") as o:
    o.write("folder\tfile_type\tfiles\tsize_MB\n")
    for (d, e), (n, s) in sorted(by_type.items()):
        o.write(f"{d}\t{e}\t{n}\t{s/1e6:.1f}\n")
with open(os.path.join(inv, "largest_files.tsv"), "w") as o:
    o.write("size_GB\tpath\n")
    for s, p in sorted(rows, reverse=True)[:50]:
        o.write(f"{s/1e9:.2f}\t{p}\n")
with open(os.path.join(inv, "busiest_dirs.tsv"), "w") as o:
    o.write("files\tfolder\n")
    for d, n in by_leaf.most_common(40):
        o.write(f"{n}\t{d}\n")
total = sum(s for s, _ in rows)
print(f"Total: {len(rows)} files, {total/1e9:.1f} GB")
PY

show() {
    echo; echo "===== $1 ====="
    if command -v column >/dev/null 2>&1; then column -t -s $'\t' "${INV}/$2" | head -"${3:-200}" || true
    else head -"${3:-200}" "${INV}/$2" || true; fi
}
show "Per folder"                  by_folder.tsv
show "Per folder and file type"    by_folder_type.tsv 300
show "50 largest files"            largest_files.tsv 52
show "Folders with the most files" busiest_dirs.tsv 42

# copy the small summaries next to the results so they download together
mkdir -p "${PROJECT}/SS_WGS_results/inventory"
cp "${INV}"/by_folder.tsv "${INV}"/by_folder_type.tsv "${INV}"/largest_files.tsv "${INV}"/busiest_dirs.tsv \
   "${PROJECT}/SS_WGS_results/inventory/"
gzip -c "${INV}/inventory.tsv" > "${PROJECT}/SS_WGS_results/inventory/inventory.tsv.gz"
echo
echo "Done $(date). Tables in ${INV}/ and SS_WGS_results/inventory/"
