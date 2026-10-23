#!/bin/bash
# ============================================================
# Step 30 (bash on the login node, ~1 min; re-run any time):
# gather all SMALL result files into one folder for download.
# Big files (FASTQ, BAM, assemblies, full VCFs) stay on MSI.
#
#   bash 30_collect_results.sh
#
# Then on your Mac:
#   rsync -av wan00965@agate.msi.umn.edu:/scratch.global/wan00965/Quick_test/SS_WGS_analysis/SS_WGS_results/ ~/SS_WGS/results/
# Missing steps are skipped, so run again after annotation / OrthoFinder finish.
# ============================================================
set -euo pipefail
cd "$(dirname "$0")"; source ./config.sh
mkdir -p log
OUT="${PROJECT}/SS_WGS_results"
exec > >(tee -a "log/30_collect_results_$(date +%Y%m%d_%H%M%S).log") 2>&1
mkdir -p "${OUT}"/{qc,assembly,variants,popgen,annotation,pangenome}

cp_if() { [[ -s "$1" ]] && cp "$1" "$2" && echo "  + ${2#${OUT}/}" || echo "  - missing: $1"; }

echo "== QC"
cp_if "${QC_DIR}/seqkit/coverage_per_isolate.tsv" "${OUT}/qc/coverage_per_isolate.tsv"
cp_if "${QC_DIR}/multiqc_report.html"             "${OUT}/qc/multiqc_report.html"
if ls "${ALN_DIR}"/*.flagstat.txt >/dev/null 2>&1; then
    { printf "isolate\tmapped_pct\tmean_depth\tbreadth_pct\n"
      for f in "${ALN_DIR}"/*.flagstat.txt; do
          s=$(basename "$f" .flagstat.txt)
          m=$(grep -m1 'primary mapped' "$f" | grep -oE '[0-9.]+%' | tr -d %)
          c="${ALN_DIR}/${s}.coverage_summary.txt"
          d=$(grep -oE 'mean_depth=[0-9.]+' "$c" 2>/dev/null | cut -d= -f2); b=$(grep -oE 'breadth_pct=[0-9.]+' "$c" 2>/dev/null | cut -d= -f2)
          printf "%s\t%s\t%s\t%s\n" "$s" "$m" "${d:-NA}" "${b:-NA}"
      done; } > "${OUT}/qc/mapping_rates.tsv" && echo "  + qc/mapping_rates.tsv"
fi

echo "== Assembly"
cp_if "${ASM_DIR}/decontam/summary.tsv"       "${OUT}/assembly/decontam_summary.tsv"
cp_if "${ASMQC_DIR}/busco_summary.tsv"        "${OUT}/assembly/busco_summary.tsv"
cp_if "${ASMQC_DIR}/quast/transposed_report.tsv" "${OUT}/assembly/quast_report.tsv"
cp_if "${ASMQC_DIR}/quast/report.html"        "${OUT}/assembly/quast_report.html"

echo "== Variants"
cp_if "${VCF_DIR}/SS_WGS.high_moderate_effects.tsv" "${OUT}/variants/high_moderate_effects.tsv"
cp_if "${VCF_DIR}/SS_WGS.snpeff_summary.html"       "${OUT}/variants/snpeff_summary.html"
cp_if "${VCF_DIR}/SS_WGS.PASS.bcftools_stats.txt"   "${OUT}/variants/bcftools_stats.txt"

echo "== Population"
for f in SS_WGS_snps.treefile SS_WGS_snps.iqtree pca.tsv pca_variance.tsv snp_distance.tsv isolates.txt; do
    cp_if "${POP_DIR}/${f}" "${OUT}/popgen/${f}"
done

echo "== Annotation"
# per-isolate gene table without the sequence columns + one summary row per isolate
python3 - "${ANNO_DIR}" "${OUT}/annotation" "${SAMPLES}" <<'PY'
import csv, glob, os, sys
anno_dir, out, samples = sys.argv[1:4]
csv.field_size_limit(10**9)
seq_cols = {"gDNA", "mRNA", "CDS-transcript", "Translation"}
rows_out = []
for s in [l.strip() for l in open(samples) if l.strip()]:
    hits = glob.glob(os.path.join(anno_dir, s, "fun", "annotate_results", "*.annotations.txt"))
    if not hits:
        continue
    with open(hits[0], newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        keep = [c for c in r.fieldnames if c not in seq_cols]
        genes = set(); n = {"secreted": 0, "cazyme": 0, "protease": 0, "pfam": 0, "membrane": 0, "go": 0}
        with open(os.path.join(out, f"{s}.annotations.tsv"), "w", newline="") as o:
            w = csv.writer(o, delimiter="\t"); w.writerow(keep)
            for row in r:
                w.writerow([row.get(c, "") for c in keep])
                if row.get("Feature", "mRNA") != "mRNA":
                    continue
                genes.add(row.get("GeneID"))
                for key, col in [("secreted", "Secreted"), ("cazyme", "CAZyme"), ("protease", "Protease"),
                                 ("pfam", "PFAM"), ("membrane", "Membrane"), ("go", "GO Terms")]:
                    if row.get(col, "").strip():
                        n[key] += 1
    rows_out.append([s, len(genes)] + [n[k] for k in ["secreted", "cazyme", "protease", "pfam", "membrane", "go"]])
    print(f"  + annotation/{s}.annotations.tsv  ({len(genes)} genes)")
if rows_out:
    with open(os.path.join(out, "annotation_summary.tsv"), "w") as o:
        o.write("isolate\tgenes\tsecreted\tcazyme\tprotease\twith_pfam\tmembrane\twith_GO\n")
        for r in rows_out:
            o.write("\t".join(map(str, r)) + "\n")
    print(f"  + annotation/annotation_summary.tsv ({len(rows_out)} isolates)")
else:
    print("  - no annotation results yet (step 22)")
PY

echo "== Pangenome"
RES=$(ls -d "${ORTHO_DIR}"/proteomes/OrthoFinder/Results_* 2>/dev/null | tail -1 || true)
if [[ -n "${RES}" ]]; then
    cp_if "${RES}/Orthogroups/Orthogroups.GeneCount.tsv"         "${OUT}/pangenome/Orthogroups.GeneCount.tsv"
    cp_if "${RES}/Orthogroups/Orthogroups_UnassignedGenes.tsv"   "${OUT}/pangenome/Orthogroups_UnassignedGenes.tsv"
    cp_if "${RES}/Comparative_Genomics_Statistics/Statistics_Overall.tsv" "${OUT}/pangenome/Statistics_Overall.tsv"
    cp_if "${RES}/Species_Tree/SpeciesTree_rooted.txt"           "${OUT}/pangenome/SpeciesTree_rooted.txt"
    cp_if "${ORTHO_DIR}/pangenome_summary.txt"                   "${OUT}/pangenome/pangenome_summary.txt"
else
    echo "  - no OrthoFinder results yet (step 23)"
fi

echo
du -sh "${OUT}"
echo "Download on your Mac with:"
echo "  rsync -av wan00965@agate.msi.umn.edu:${OUT}/ ~/SS_WGS/results/"
