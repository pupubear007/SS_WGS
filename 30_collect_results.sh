#!/bin/bash
# ============================================================
# Step 30 (bash on the login node, a few minutes; re-run any time):
# gather ALL useful results into one folder for download to the Mac.
#
#   bash 30_collect_results.sh                    # tables, reports, VCFs, gene models (~0.5-1 GB)
#   bash 30_collect_results.sh --with-assemblies  # + 44 cleaned assemblies (~1.6 GB more)
#
# Left on MSI on purpose (too large): FASTQ, BAMs, gVCFs, SPAdes work dirs,
# funannotate work dirs, OrthoFinder sequence/tree folders.
#
# Then on your Mac:
#   rsync -av wan00965@agate.msi.umn.edu:/scratch.global/wan00965/Quick_test/SS_WGS_analysis/SS_WGS_results/ ~/SS_WGS/results/
#
# Folder names qc/ assembly/ variants/ popgen/ annotation/ pangenome/ are the
# ones R/SS_WGS_results.Rmd reads, so the report keeps working.
# ============================================================
set -euo pipefail
cd "$(dirname "$0")"; SCRIPT_DIR=$(pwd); source ./config.sh
mkdir -p log
exec > >(tee -a "log/30_collect_results_$(date +%Y%m%d_%H%M%S).log") 2>&1

WITH_ASM=0; [[ "${1:-}" == "--with-assemblies" ]] && WITH_ASM=1
OUT="${PROJECT}/SS_WGS_results"
mkdir -p "${OUT}"/{reference,qc/fastp,qc/multiqc_data,alignment,variants,popgen,assembly/decontam,assembly/busco,assembly/quast,annotation/per_isolate,pangenome/orthologues_vs_reference,logs}
echo "Collecting into ${OUT}   $(date)"

N_OK=0; N_MISS=0
# cp_f SRC DEST_DIR [NEW_NAME] : copy one file if it exists
cp_f() {
    local src="$1" dst="$2" name="${3:-$(basename "$1")}"
    if [[ -s "${src}" ]]; then cp -pf "${src}" "${dst}/${name}"; N_OK=$((N_OK+1))
    else echo "  - missing: ${src#${PROJECT}/}"; N_MISS=$((N_MISS+1)); fi
}
# cp_g DEST_DIR GLOB... : copy every file matching the globs (quiet if none)
cp_g() {
    local dst="$1"; shift; local n=0 f
    for f in "$@"; do [[ -f "$f" ]] && cp -pf "$f" "${dst}/" && n=$((n+1)); done
    echo "  + ${dst#${OUT}/}: ${n} files"; N_OK=$((N_OK+n))
}

# ---------- reference ----------
echo "== reference"
cp_g "${OUT}/reference" "${REF_DIR}"/Ss_ref.gff "${REF_DIR}"/Ss_ref.gtf "${REF_DIR}"/Ss_ref_proteins.fa "${REF_DIR}"/Ss_ref.fasta.fai
[[ ${WITH_ASM} -eq 1 ]] && cp_f "${REF}" "${OUT}/reference"

# ---------- QC ----------
echo "== QC (steps 02-03)"
cp_f "${QC_DIR}/seqkit/coverage_per_isolate.tsv" "${OUT}/qc"
cp_f "${QC_DIR}/seqkit/raw_stats.tsv"            "${OUT}/qc"
cp_f "${QC_DIR}/seqkit/trimmed_stats.tsv"        "${OUT}/qc"
cp_f "${QC_DIR}/multiqc_report.html"             "${OUT}/qc"
cp_g "${OUT}/qc/multiqc_data" "${QC_DIR}"/multiqc_report_data/*.txt "${QC_DIR}"/multiqc_data/*.txt
cp_g "${OUT}/qc/fastp" "${QC_DIR}"/fastp/*.json "${QC_DIR}"/fastp/*.html

# ---------- alignment ----------
echo "== alignment (step 11)"
cp_g "${OUT}/alignment" "${ALN_DIR}"/*.flagstat.txt "${ALN_DIR}"/*.coverage_summary.txt \
                        "${ALN_DIR}"/*.coverage_per_contig.tsv "${ALN_DIR}"/*.markdup_stats.txt
if ls "${ALN_DIR}"/*.flagstat.txt >/dev/null 2>&1; then
    { printf "isolate\tmapped_pct\tmean_depth\tbreadth_pct\tduplicate_pct\n"
      for f in "${ALN_DIR}"/*.flagstat.txt; do
          s=$(basename "$f" .flagstat.txt)
          m=$(grep -m1 'primary mapped' "$f" | grep -oE '[0-9.]+%' | tr -d %)
          tot=$(grep -m1 'primary$' "$f" | awk '{print $1}'); dup=$(grep -m1 'primary duplicates' "$f" | awk '{print $1}')
          dp=$(awk -v t="${tot:-0}" -v d="${dup:-0}" 'BEGIN{if(t>0) printf "%.2f",100*d/t; else print "NA"}')
          c="${ALN_DIR}/${s}.coverage_summary.txt"
          d=$(grep -oE 'mean_depth=[0-9.]+' "$c" 2>/dev/null | cut -d= -f2); b=$(grep -oE 'breadth_pct=[0-9.]+' "$c" 2>/dev/null | cut -d= -f2)
          printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$m" "${d:-NA}" "${b:-NA}" "$dp"
      done; } > "${OUT}/qc/mapping_rates.tsv"
    echo "  + qc/mapping_rates.tsv"
fi

# ---------- variants ----------
echo "== variants (step 12)"
for f in SS_WGS.filtered.vcf.gz SS_WGS.PASS.vcf.gz SS_WGS.snps.clean.vcf.gz SS_WGS.PASS.snpeff.vcf.gz; do
    cp_f "${VCF_DIR}/${f}" "${OUT}/variants"
    [[ -s "${VCF_DIR}/${f}.tbi" ]] && cp -p "${VCF_DIR}/${f}.tbi" "${OUT}/variants/"
done
cp_f "${VCF_DIR}/SS_WGS.PASS.bcftools_stats.txt"   "${OUT}/variants" bcftools_stats.txt
cp_g "${OUT}/variants" "${VCF_DIR}"/SS_WGS.snpeff_summary* "${VCF_DIR}"/*snpEff_genes.txt "${VCF_DIR}"/*snpeff*genes.txt

# ---------- population ----------
echo "== population (step 13)"
cp_g "${OUT}/popgen" "${POP_DIR}"/*.tsv "${POP_DIR}"/isolates.txt "${POP_DIR}"/SS_WGS_snps.fasta \
     "${POP_DIR}"/SS_WGS_snps.{iqtree,treefile,contree,log,mldist,bionj,model.gz}

# ---------- assemblies ----------
echo "== assembly (steps 20, 20b, 21)"
cp_f "${ASM_DIR}/decontam/summary.tsv" "${OUT}/assembly" decontam_summary.tsv
cp_g "${OUT}/assembly/decontam" "${ASM_DIR}"/decontam/*.contigs.tsv "${ASM_DIR}"/decontam/*.unaligned_hits.tsv
cp_f "${ASMQC_DIR}/busco_summary.tsv" "${OUT}/assembly"
cp_g "${OUT}/assembly/busco" "${ASMQC_DIR}"/busco_short/*.txt "${ASMQC_DIR}"/busco_short/*.png
cp_f "${ASMQC_DIR}/quast/transposed_report.tsv" "${OUT}/assembly" quast_report.tsv
cp_g "${OUT}/assembly/quast" "${ASMQC_DIR}"/quast/report.{html,tsv,txt,pdf} "${ASMQC_DIR}"/quast/transposed_report.tsv \
     "${ASMQC_DIR}"/quast/genome_stats/*.txt "${ASMQC_DIR}"/quast/aligned_stats/*.pdf
# assembly stats for the cleaned assemblies (sizes after decontamination)
if ls "${ASM_DIR}"/clean/*.scaffolds.fa >/dev/null 2>&1; then
    activate_env "${ENV_ASM}" >/dev/null
    seqkit stats -a -T -j 8 "${ASM_DIR}"/clean/*.scaffolds.fa > "${OUT}/assembly/clean_assembly_stats.tsv" \
        && echo "  + assembly/clean_assembly_stats.tsv"
    set +u; conda deactivate 2>/dev/null || true; set -u
fi
if [[ ${WITH_ASM} -eq 1 ]]; then
    mkdir -p "${OUT}/assembly/clean_fasta"
    cp_g "${OUT}/assembly/clean_fasta" "${ASM_DIR}"/clean/*.scaffolds.fa
fi

# ---------- annotation ----------
echo "== annotation (step 22)"
python3 - "${ANNO_DIR}" "${OUT}/annotation" "${SAMPLES}" <<'PY'
import csv, glob, json, os, shutil, sys
anno_dir, out, samples = sys.argv[1:4]
csv.field_size_limit(10**9)
seq_cols = {"gDNA", "mRNA", "CDS-transcript", "Translation"}
summary = []
for s in [l.strip() for l in open(samples) if l.strip()]:
    res = os.path.join(anno_dir, s, "fun", "annotate_results")
    ann = glob.glob(os.path.join(res, "*.annotations.txt"))
    if not ann:
        print(f"  - missing annotation: {s}"); continue
    d = os.path.join(out, "per_isolate", s); os.makedirs(d, exist_ok=True)
    # gene models, proteins, CDS and stats for local use (IGV, Geneious, effector tools)
    for pat in ("*.gff3", "*.proteins.fa", "*.cds-transcripts.fa", "*.stats.json"):
        for f in glob.glob(os.path.join(res, pat)):
            shutil.copy2(f, os.path.join(d, f"{s}{os.path.splitext(f)[1] if pat=='*.gff3' else '.' + pat.split('.',1)[1]}"))
    # annotation table without sequence columns
    n = dict(secreted=0, cazyme=0, protease=0, pfam=0, interpro=0, eggnog=0, membrane=0, go=0, named=0)
    genes = set()
    with open(ann[0], newline="") as fh, open(os.path.join(d, f"{s}.annotations.tsv"), "w", newline="") as o:
        r = csv.DictReader(fh, delimiter="\t"); keep = [c for c in r.fieldnames if c not in seq_cols]
        w = csv.writer(o, delimiter="\t"); w.writerow(keep)
        for row in r:
            w.writerow([row.get(c, "") for c in keep])
            if row.get("Feature", "mRNA") != "mRNA": continue
            genes.add(row.get("GeneID"))
            for k, col in [("secreted","Secreted"),("cazyme","CAZyme"),("protease","Protease"),("pfam","PFAM"),
                           ("interpro","InterPro"),("eggnog","EggNog"),("membrane","Membrane"),("go","GO Terms")]:
                if row.get(col, "").strip(): n[k] += 1
            if row.get("Product", "").strip() and row["Product"] != "hypothetical protein": n["named"] += 1
    # tRNA / gene counts from funannotate stats if present
    trna = ""
    for f in glob.glob(os.path.join(res, "*.stats.json")):
        try: trna = json.load(open(f)).get("annotation", {}).get("tRNA", "")
        except Exception: pass
    summary.append([s, len(genes), n["named"], n["pfam"], n["interpro"], n["eggnog"], n["go"],
                    n["secreted"], n["membrane"], n["cazyme"], n["protease"], trna])
if summary:
    with open(os.path.join(out, "annotation_summary.tsv"), "w") as o:
        o.write("isolate\tgenes\twith_product_name\twith_pfam\twith_interpro\twith_eggnog\twith_GO\t"
                "secreted\tmembrane\tcazyme\tprotease\ttRNA\n")
        for r in summary: o.write("\t".join(map(str, r)) + "\n")
    print(f"  + annotation/annotation_summary.tsv and per_isolate/ ({len(summary)} isolates)")
    empty = [c for i, c in enumerate(["with_interpro","with_eggnog","with_GO","secreted","membrane"], start=4)
             if all(r[i] == 0 for r in summary)]
    if empty:
        print("  (funannotate did not add", ", ".join(empty) + "; secreted / InterPro / KOG labels come from gene_level/gene_info.tsv, step 14)")
PY
# keep the old flat layout too, so the current report still finds them
cp -p "${OUT}"/annotation/per_isolate/*/*.annotations.tsv "${OUT}/annotation/" 2>/dev/null || true

# ---------- pangenome ----------
echo "== pangenome (step 23)"
RES=$(ls -d "${ORTHO_DIR}"/proteomes/OrthoFinder/Results_* 2>/dev/null | tail -1 || true)
if [[ -n "${RES}" ]]; then
    cp_g "${OUT}/pangenome" "${RES}"/Orthogroups/Orthogroups.tsv "${RES}"/Orthogroups/Orthogroups.GeneCount.tsv \
         "${RES}"/Orthogroups/Orthogroups_UnassignedGenes.tsv "${RES}"/Orthogroups/Orthogroups_SingleCopyOrthologues.txt \
         "${RES}"/Comparative_Genomics_Statistics/*.tsv "${RES}"/Species_Tree/SpeciesTree_rooted*.txt \
         "${RES}"/Log.txt "${ORTHO_DIR}/pangenome_summary.txt"
    # one-to-one ortholog tables: reference (SS1G_ IDs) vs each isolate (FUN_ IDs) -> links to RNA-seq gene IDs
    REFDIR=$(ls -d "${RES}"/Orthologues/Orthologues_*ref* 2>/dev/null | head -1 || true)
    [[ -n "${REFDIR}" ]] && cp_g "${OUT}/pangenome/orthologues_vs_reference" "${REFDIR}"/*.tsv
else
    echo "  - no OrthoFinder results yet (step 23)"
fi

# ---------- logs ----------
echo "== logs"
cp_g "${OUT}/logs" "${SCRIPT_DIR}"/log/*.out "${SCRIPT_DIR}"/log/*.err "${SCRIPT_DIR}"/log/*.log
for s in $(cat "${SAMPLES}"); do
    f="${ANNO_DIR}/${s}/fun/logfiles/funannotate-predict.log"
    [[ -s "$f" ]] && cp -p "$f" "${OUT}/logs/${s}.funannotate-predict.log"
done

# ---------- manifest + readme ----------
( cd "${OUT}" && find . -type f ! -name MANIFEST.tsv -printf '%s\t%p\n' | sort -k2 \
    | awk -F'\t' 'BEGIN{print "size_MB\tpath"} {printf "%.2f\t%s\n",$1/1e6,$2}' > MANIFEST.tsv )
cat > "${OUT}/README.txt" <<'EOF'
SS_WGS results (made by 30_collect_results.sh)
reference/   Ss_ref gene models (GFF/GTF) and proteins (SS1G_ IDs, same as RNA-seq)
qc/          coverage per isolate, mapping rates, MultiQC report, fastp reports
alignment/   per-isolate flagstat, depth, per-contig coverage, duplicate stats
variants/    VCFs (filtered, PASS, clean SNPs, SnpEff), bcftools stats, gene-effect table
popgen/      PCA, SNP distance matrix, SNP alignment, IQ-TREE tree + report
assembly/    decontamination summary + per-contig calls, BUSCO, QUAST, cleaned-assembly stats
             clean_fasta/ only with --with-assemblies
annotation/  annotation_summary.tsv; per_isolate/<ISO>/ gff3, proteins, CDS, annotation table
pangenome/   OrthoFinder orthogroups (Orthogroups.tsv = gene IDs per family), statistics,
             species tree; orthologues_vs_reference/ = FUN_ <-> SS1G_ gene ID links
gene_level/  (from 14_gene_level.sh) variant effects per gene, gene x isolate HIGH/MODERATE/promoter
             variant counts, relative gene depth (deletions), gene function (SignalP, KOG, InterPro)
logs/        SLURM logs and funannotate predict logs
MANIFEST.tsv every file with its size
EOF

echo
echo "Copied ${N_OK} files; ${N_MISS} expected files missing (listed above)."
du -sh "${OUT}"
du -sh "${OUT}"/*/ | sort -k2
echo
echo "Download on your Mac with:"
echo "  rsync -av wan00965@agate.msi.umn.edu:${OUT}/ ~/SS_WGS/results/"
