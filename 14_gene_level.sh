#!/bin/bash
#SBATCH --job-name=14_gene_level
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

# ============================================================
# Step 14 (single job, after steps 11-12): gene-level genomic variation,
# keyed by the SAME JGI gene IDs used in the RNA-seq (SS1G_xxxxx / SS1G_xxxxxT0).
#
#   sbatch 14_gene_level.sh
#
# 1. converts the JGI 2011 GFF (GTF-like, name "SS1G_..T0") to a proper GTF
#    (gffread could not, so the step-10 SnpEff database was never built)
# 2. builds SnpEff database + annotates the PASS VCF  -> protein-level effects
# 3. per gene x isolate: HIGH / MODERATE variants, promoter (<=1 kb upstream) variants
# 4. per gene x isolate: read depth relative to the isolate's genome median
#    (~0 = gene deleted, ~1 = single copy, >=1.8 = duplicated)
# 5. gene_info: protein length, SignalP / small secreted / Cys-rich, KOG, InterPro, KEGG
#    from the JGI MycoCosm annotation files (no extra tools needed)
#
# Outputs: SS_WGS_analysis/14_gene_level/ and copied to SS_WGS_results/gene_level/
#   variant_effects.tsv.gz        one row per variant x gene, genotype string for all isolates
#   gene_isolate_HIGH.tsv         count of HIGH-impact ALT alleles each isolate carries
#   gene_isolate_MODERATE.tsv
#   gene_isolate_promoter1kb.tsv
#   gene_isolate_relcov.tsv       relative read depth (copy number / deletion)
#   gene_info.tsv                 gene coordinates + function
#   SS_WGS.PASS.snpeff.vcf.gz     annotated VCF (also into 05_vcf/)
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_WGS}"

# JGI MycoCosm download folder with *_SigP.tab.gz, *_IPR.tab.gz, *_KOG.tab.gz, *_KEGG.tab.gz
JGI_DIR="${JGI_DIR:-/scratch.global/wan00965/JGI_ref_genome_n_transcripts/Ss_genome}"

GL="${PROJECT}/14_gene_level"
mkdir -p "${GL}"
cd "${GL}"
echo "Start $(date)"

# ---------- 1. proper GTF + gene BED from the JGI GFF ----------
awk -F'\t' 'BEGIN{OFS="\t"}
    !/^#/ && NF>=9 && ($3=="exon"||$3=="CDS"||$3=="start_codon"||$3=="stop_codon") {
        if (!match($9, /name "[^"]+"/)) next
        t = substr($9, RSTART+6, RLENGTH-7); g = t; sub(/T[0-9]+$/, "", g)
        print $1,$2,$3,$4,$5,$6,$7,$8, "gene_id \"" g "\"; transcript_id \"" t "\"; gene_name \"" g "\";"
    }' "${REF_DIR}/Ss_ref.gff" > Ss_ref.fixed.gtf

awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" {
        match($9, /transcript_id "[^"]+"/); t = substr($9, RSTART+15, RLENGTH-16)
        g = t; sub(/T[0-9]+$/, "", g); k = $1 SUBSEP g SUBSEP t SUBSEP $7
        if (!(k in s) || $4 < s[k]) s[k] = $4
        if (!(k in e) || $5 > e[k]) e[k] = $5
    } END { for (k in s) { split(k, a, SUBSEP); print a[1], s[k]-1, e[k], a[2], a[3], a[4] } }' \
    Ss_ref.fixed.gtf | sort -k1,1 -k2,2n > genes.bed
echo "GTF lines: $(wc -l < Ss_ref.fixed.gtf); genes: $(wc -l < genes.bed)"
[[ $(wc -l < genes.bed) -gt 10000 ]] || { echo "ERROR: too few genes parsed from ${REF_DIR}/Ss_ref.gff"; exit 1; }
cp Ss_ref.fixed.gtf "${REF_DIR}/Ss_ref.gtf"      # replace the empty step-10 GTF

# ---------- 2. SnpEff database + annotation ----------
SE="${GL}/snpeff"
mkdir -p "${SE}/data/Sclsc"
ln -sf "${REF}" "${SE}/data/Sclsc/sequences.fa"
cp Ss_ref.fixed.gtf "${SE}/data/Sclsc/genes.gtf"
printf "data.dir = %s/data/\nSclsc.genome : Sclerotinia_sclerotiorum_1980_UF-70_Sclsc1\n" "${SE}" > "${SE}/snpEff.config"
snpEff -Xmx16g build -c "${SE}/snpEff.config" -gtf22 -noCheckCds -noCheckProtein -v Sclsc > snpeff_build.log 2>&1 \
    || { tail -30 snpeff_build.log; echo "ERROR: SnpEff build failed"; exit 1; }
echo "SnpEff database built"

snpEff -Xmx16g ann -c "${SE}/snpEff.config" -noLog -stats SS_WGS.snpeff_summary.html \
    Sclsc "${VCF_DIR}/SS_WGS.PASS.vcf.gz" | bgzip -@ 4 > SS_WGS.PASS.snpeff.vcf.gz
tabix -f -p vcf SS_WGS.PASS.snpeff.vcf.gz
cp SS_WGS.PASS.snpeff.vcf.gz SS_WGS.PASS.snpeff.vcf.gz.tbi SS_WGS.snpeff_summary* "${VCF_DIR}/" 2>/dev/null || true
echo "Annotated variants: $(bcftools view -H SS_WGS.PASS.snpeff.vcf.gz | wc -l)"

# ---------- 3. variant x gene table + gene x isolate matrices ----------
python3 - SS_WGS.PASS.snpeff.vcf.gz genes.bed <<'PY'
import gzip, sys, collections
vcf, bed = sys.argv[1], sys.argv[2]
rank = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
genes = sorted({l.split("\t")[3] for l in open(bed)})
samples = None
mat = {k: collections.defaultdict(lambda: None) for k in ("HIGH", "MODERATE", "promoter1kb")}
out = gzip.open("variant_effects.tsv.gz", "wt")
out.write("chrom\tpos\tref\talt\tgene\tcategory\teffect\timpact\thgvs_c\thgvs_p\tdistance\tn_alt\tn_called\tgenotypes\n")
n_rows = 0
for line in gzip.open(vcf, "rt"):
    if line.startswith("##"): continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#"):
        samples = f[9:]
        for k in mat:
            mat[k] = {g: [0] * len(samples) for g in genes}
        continue
    chrom, pos, ref, alts = f[0], f[1], f[3], f[4].split(",")
    info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in f[7].split(";"))
    if "ANN" not in info: continue
    gts = [s.split(":", 1)[0].replace("|", "/").split("/")[0] for s in f[9:]]
    best = {}   # (allele, gene) -> annotation fields of most severe effect
    for a in info["ANN"].split(","):
        p = a.split("|")
        if len(p) < 15: continue
        allele, eff, imp, gene, ftype = p[0], p[1], p[2], p[3], p[5]
        if not gene.startswith("SS1G"): continue
        dist = p[14]
        if "upstream_gene_variant" in eff:
            if not dist or int(dist) > 1000: continue
            cat = "promoter1kb"
        elif "downstream_gene_variant" in eff or "intergenic" in eff:
            continue
        elif imp in ("HIGH", "MODERATE"):
            cat = imp
        elif imp == "LOW":
            cat = "LOW"
        else:
            cat = "genic_noncoding"          # intron / UTR / splice-region modifiers
        order = {"HIGH": 5, "MODERATE": 4, "LOW": 3, "genic_noncoding": 2, "promoter1kb": 1}[cat]
        k = (allele, gene)
        if k not in best or order > best[k][0]:
            best[k] = (order, cat, eff, imp, p[9], p[10], dist)
    for (allele, gene), (_, cat, eff, imp, hc, hp, dist) in best.items():
        if allele not in alts: continue
        idx = str(alts.index(allele) + 1)
        g = "".join("1" if x == idx else ("." if x in (".", "") else ("0" if x == "0" else "o")) for x in gts)
        n_alt = g.count("1"); n_called = len(g) - g.count(".")
        out.write(f"{chrom}\t{pos}\t{ref}\t{allele}\t{gene}\t{cat}\t{eff}\t{imp}\t{hc}\t{hp}\t{dist}\t{n_alt}\t{n_called}\t{g}\n")
        n_rows += 1
        if cat in mat and gene in mat[cat]:
            row = mat[cat][gene]
            for i, c in enumerate(g):
                if c == "1": row[i] += 1
out.close()
for k, m in mat.items():
    with open(f"gene_isolate_{k}.tsv", "w") as o:
        o.write("gene\t" + "\t".join(samples) + "\n")
        for g in genes:
            o.write(g + "\t" + "\t".join(map(str, m[g])) + "\n")
print(f"variant x gene rows: {n_rows}; genes: {len(genes)}; isolates: {len(samples)}")
PY

# ---------- 4. relative read depth per gene (deletions / duplications) ----------
ls "${ALN_DIR}"/*.markdup.bam | sort > bams.txt
mkdir -p cov
export GL
cat bams.txt | xargs -P "${THREADS}" -I{} sh -c \
    's=$(basename {} .markdup.bam); samtools bedcov -Q 10 "${GL}/genes.bed" {} | cut -f7 > "${GL}/cov/${s}.sum"'
python3 - genes.bed bams.txt <<'PY'
import os, statistics, sys
bed = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
iso = [os.path.basename(l.strip()).replace(".markdup.bam", "") for l in open(sys.argv[2])]
lens = [int(b[2]) - int(b[1]) for b in bed]
cols = {}
for s in iso:
    sums = [float(x) for x in open(f"cov/{s}.sum")]
    depth = [v / L if L > 0 else 0 for v, L in zip(sums, lens)]
    med = statistics.median([d for d in depth if d > 0]) or 1
    cols[s] = [d / med for d in depth]
with open("gene_isolate_relcov.tsv", "w") as o:
    o.write("gene\t" + "\t".join(iso) + "\n")
    for i, b in enumerate(bed):
        o.write(b[3] + "\t" + "\t".join(f"{cols[s][i]:.3f}" for s in iso) + "\n")
absent = {s: sum(v < 0.1 for v in cols[s]) for s in iso}
print("genes with relative depth < 0.1 (likely deleted), per isolate:")
for s in iso: print(f"  {s}\t{absent[s]}")
PY

# ---------- 5. gene function from JGI annotation files ----------
python3 - genes.bed "${REF_DIR}/Ss_ref_proteins.fa" "${JGI_DIR}" <<'PY'
import glob, gzip, os, sys, collections
bed, prot, jgi = sys.argv[1:4]
op = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)
# protein id -> transcript, sequence
pid2t, seq, cur = {}, {}, None
for l in open(prot):
    if l.startswith(">"):
        h = l[1:].strip().split("|"); cur = h[-1]; pid2t[h[2] if len(h) > 2 else cur] = cur; seq[cur] = []
    elif cur: seq[cur].append(l.strip())
seq = {t: "".join(s).rstrip("*") for t, s in seq.items()}
def find(pat):
    hits = [f for f in glob.glob(os.path.join(jgi, "**", pat), recursive=True)]
    return hits[0] if hits else None
def table(f):
    rows = []
    with op(f) as fh:
        hl = None
        for l in fh:
            if hl is None:                       # header = first line naming the protein id column
                if "proteinid" in l.lower().replace(" ", ""):
                    hl = [h.strip().lower() for h in l.lstrip("#").rstrip("\n").split("\t")]
                continue
            if l.startswith("#") or not l.strip():
                continue
            rows.append(dict(zip(hl, l.rstrip("\n").split("\t"))))
    return rows
sigp, ipr, kog, kegg = set(), collections.defaultdict(set), {}, collections.defaultdict(set)
f = find("*SigP.tab*")
if f: sigp = {r.get("proteinid") for r in table(f)}
f = find("*IPR.tab*")
if f:
    for r in table(f):
        if r.get("iprdesc"): ipr[r.get("proteinid")].add(r["iprdesc"])
f = find("*KOG.tab*")
if f:
    for r in table(f): kog[r.get("proteinid")] = (r.get("kogclass", ""), r.get("kogdefline", ""))
f = find("*KEGG.tab*")
if f:
    for r in table(f):
        if r.get("definition"): kegg[r.get("proteinid")].add(r["definition"])
print(f"JGI annotation: SignalP {len(sigp)}, InterPro {len(ipr)}, KOG {len(kog)}, KEGG {len(kegg)} proteins")
t2pid = {t: p for p, t in pid2t.items()}
with open("gene_info.tsv", "w") as o:
    o.write("gene\ttranscript\tprotein_id\tcontig\tstart\tend\tstrand\taa_length\tcysteines\t"
            "signalP\tsmall_secreted\tcys_rich_SSP\tkog_class\tkog_def\tinterpro\tkegg\n")
    for l in open(bed):
        c, s, e, g, t, strand = l.rstrip("\n").split("\t")
        p = t2pid.get(t, ""); aa = seq.get(t, "")
        sp = p in sigp; L = len(aa); cys = aa.count("C")
        ssp = sp and 0 < L <= 300; crich = ssp and cys >= 4
        kc, kd = kog.get(p, ("", ""))
        o.write("\t".join(map(str, [g, t, p, c, int(s) + 1, e, strand, L, cys, int(sp), int(ssp), int(crich),
                                    kc, kd, "; ".join(sorted(ipr.get(p, []))), "; ".join(sorted(kegg.get(p, [])))])) + "\n")
print("gene_info.tsv written")
PY

# ---------- copy to the results folder ----------
OUTGL="${PROJECT}/SS_WGS_results/gene_level"
mkdir -p "${OUTGL}"
cp variant_effects.tsv.gz gene_isolate_*.tsv gene_info.tsv genes.bed Ss_ref.fixed.gtf SS_WGS.snpeff_summary* "${OUTGL}/" 2>/dev/null || true
echo
echo "Files in ${OUTGL}:"; ls -la "${OUTGL}"
echo "Done $(date)"
