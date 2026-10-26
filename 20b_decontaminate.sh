#!/bin/bash
#SBATCH --job-name=20b_decontaminate
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 20b (single job, after step 20): find and remove non-Sclerotinia
# contigs from every assembly, and report what the contaminant is.
#
# Per contig it measures:
#   aln_frac  fraction of the contig that aligns to the 1980 UF-70 reference (minimap2)
#   gc        GC content
#   cov       SPAdes k-mer coverage, relative to that isolate's median
# Rule (conservative, keeps isolate-specific Sclerotinia DNA):
#   KEEP   aln_frac >= 0.30
#   KEEP   unaligned but GC within +/-4% of the isolate's Ss GC AND cov >= 0.25x median
#   REMOVE everything else
# Unaligned contigs >= 1 kb are also searched (DIAMOND blastx) against the
# UniProt/SwissProt DB that funannotate installed, to name the contaminant.
#
#   sbatch 20b_decontaminate.sh
# Outputs (07_assembly/decontam/):
#   summary.tsv                one row per isolate: size before/after, Mb removed, top contaminant
#   <ISO>.contigs.tsv          every contig with aln_frac, gc, cov, decision
#   <ISO>.unaligned_hits.tsv   DIAMOND top hit per unaligned contig (organism in "OS=")
#   ../clean/<ISO>.scaffolds.fa  cleaned assemblies  -> used by step 22
# Originals in 07_assembly/ are not changed.
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_ASM}"

DEC="${ASM_DIR}/decontam"
CLEAN="${ASM_DIR}/clean"
mkdir -p "${DEC}" "${CLEAN}"

# tools: minimap2 into the assembly env if missing; DIAMOND comes from pu_funannotate
command -v minimap2 >/dev/null 2>&1 || mamba install -y -n "${ENV_ASM}" -c conda-forge -c bioconda minimap2
DIAMOND="$(conda info --base)/envs/${ENV_ANNO}/bin/diamond"
UNIPROT=$(ls "${FUNANNOTATE_DB}"/uniprot*.dmnd 2>/dev/null | head -1 || true)
[[ -x "${DIAMOND}" && -n "${UNIPROT}" ]] || { echo "NOTE: DIAMOND/UniProt not found, contaminant naming skipped"; DIAMOND=""; }

printf "isolate\tcontigs\tsize_Mb\tremoved_contigs\tremoved_Mb\tclean_Mb\tSs_GC\ttop_contaminant_hits\n" > "${DEC}/summary.tsv"

while read -r S; do
    ASM="${ASM_DIR}/${S}.scaffolds.fa"
    SPADES="${ASM_DIR}/${S}/spades/scaffolds.fasta"   # original headers carry coverage
    [[ -s "${ASM}" && -s "${SPADES}" ]] || { echo "SKIP ${S}: assembly files missing"; continue; }
    echo "== ${S}  $(date +%T)"

    # 1. align contigs to the reference
    minimap2 -x asm10 -c --secondary=no -t "${THREADS}" "${REF}" "${ASM}" > "${DEC}/${S}.paf" 2>/dev/null

    # 2. classify every contig
    python3 - "${S}" "${ASM}" "${SPADES}" "${DEC}/${S}.paf" "${DEC}" "${CLEAN}" <<'PY'
import sys, statistics
S, asm, spades, paf, dec, clean = sys.argv[1:7]

def fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name: yield name, "".join(seq)

contigs = list(fasta(asm))
# SPAdes coverage: same order as the >=500 bp filter used in step 20
cov_list = []
for h, s in fasta(spades):
    if len(s) >= 500:
        try: cov_list.append(float(h.split("_cov_")[1].split("_")[0]))
        except (IndexError, ValueError): cov_list.append(float("nan"))
if len(cov_list) != len(contigs):
    print(f"WARNING {S}: coverage list ({len(cov_list)}) != contigs ({len(contigs)}); coverage ignored")
    cov_list = [float("nan")] * len(contigs)

# aligned bases per contig (merge overlapping blocks)
blocks = {}
for line in open(paf):
    f = line.split("\t")
    blocks.setdefault(f[0], []).append((int(f[2]), int(f[3])))
def covered(iv):
    tot, end = 0, -1
    for a, b in sorted(iv):
        if b <= end: continue
        tot += b - max(a, end); end = b
    return tot

rows = []
for (name, seq), cov in zip(contigs, cov_list):
    L = len(seq); s = seq.upper()
    gc = (s.count("G") + s.count("C")) / max(1, L - s.count("N"))
    aln = covered(blocks.get(name, [])) / L
    rows.append([name, L, aln, gc, cov])

ss = [r for r in rows if r[2] >= 0.30]
ss_gc  = statistics.median([r[3] for r in ss]) if ss else 0.42
_c = [r[4] for r in ss if r[4] == r[4]]
ss_cov = statistics.median(_c) if _c else float("nan")

keep = set()
with open(f"{dec}/{S}.contigs.tsv", "w") as out:
    out.write("contig\tlength\taln_frac\tgc\tcov_rel\tdecision\n")
    for name, L, aln, gc, cov in rows:
        rel = cov / ss_cov if ss_cov == ss_cov and ss_cov > 0 and cov == cov else float("nan")
        if aln >= 0.30:
            d = "keep_aligned"
        elif abs(gc - ss_gc) <= 0.04 and (rel != rel or rel >= 0.25):
            d = "keep_Ss_like"
        else:
            d = "remove"
        if d != "remove": keep.add(name)
        out.write(f"{name}\t{L}\t{aln:.2f}\t{gc:.3f}\t{rel:.2f}\t{d}\n")

with open(f"{clean}/{S}.scaffolds.fa", "w") as out:
    for name, seq in contigs:
        if name in keep:
            out.write(f">{name}\n" + "\n".join(seq[i:i+80] for i in range(0, len(seq), 80)) + "\n")
# unaligned contigs >= 1 kb for DIAMOND
with open(f"{dec}/{S}.unaligned.fa", "w") as out:
    for (name, seq), r in zip(contigs, rows):
        if r[2] < 0.30 and len(seq) >= 1000:
            out.write(f">{name}\n{seq}\n")

tot = sum(r[1] for r in rows); rem = [r for r in rows if r[0] not in keep]
rm_bp = sum(r[1] for r in rem)
with open(f"{dec}/{S}.stats", "w") as out:
    out.write(f"{len(rows)}\t{tot/1e6:.2f}\t{len(rem)}\t{rm_bp/1e6:.2f}\t{(tot-rm_bp)/1e6:.2f}\t{100*ss_gc:.1f}")
print(f"{S}: {len(rows)} contigs, {tot/1e6:.2f} Mb; removed {len(rem)} contigs = {rm_bp/1e6:.2f} Mb; Ss GC {100*ss_gc:.1f}%")
PY

    # 3. name the contaminant (top organisms among unaligned-contig hits)
    TOP="NA"
    if [[ -n "${DIAMOND}" && -s "${DEC}/${S}.unaligned.fa" ]]; then
        "${DIAMOND}" blastx -d "${UNIPROT}" -q "${DEC}/${S}.unaligned.fa" \
            --fast --max-target-seqs 1 --evalue 1e-10 -p "${THREADS}" --quiet \
            -f 6 qseqid pident length evalue bitscore stitle \
            -o "${DEC}/${S}.unaligned_hits.tsv"
        TOP=$(grep -oE 'OS=[^=]+ OX=' "${DEC}/${S}.unaligned_hits.tsv" 2>/dev/null \
              | sed -E 's/OS=//; s/ OX=//' | awk '{print $1" "$2}' \
              | sort | uniq -c | sort -rn | head -3 \
              | awk '{printf "%s%s %s (%d)", (NR>1?"; ":""), $2, $3, $1}')
        [[ -n "${TOP}" ]] || TOP="no hits"
    fi
    printf "%s\t%s\t%s\n" "${S}" "$(cat "${DEC}/${S}.stats")" "${TOP}" >> "${DEC}/summary.tsv"
    rm -f "${DEC}/${S}.stats" "${DEC}/${S}.unaligned.fa"
done < "${SAMPLES}"

echo
echo "===== Decontamination summary ====="
column -t -s $'\t' "${DEC}/summary.tsv"
echo
echo "Clean assemblies : ${CLEAN}/"
echo "Per-contig tables: ${DEC}/<ISOLATE>.contigs.tsv"
