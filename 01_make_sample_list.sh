#!/bin/bash
# ============================================================
# Step 01 (bash, seconds): make samples.txt + output folders
# + copy/rename the reference.
#   bash 01_make_sample_list.sh
# ============================================================
set -euo pipefail
cd "$(dirname "$0")"; source ./config.sh
mkdir -p log

mkdir -p "${PROJECT}" "${QC_DIR}"/{fastqc_raw,fastp,seqkit} "${TRIM_DIR}" \
         "${ALN_DIR}" "${GVCF_DIR}" "${VCF_DIR}" "${POP_DIR}" \
         "${ASM_DIR}" "${ASMQC_DIR}" "${ANNO_DIR}" "${ORTHO_DIR}" "${REF_DIR}"

# isolate names from <ISOLATE>_R1.fq.gz
ls "${RAW_DIR}"/*_R1.fq.gz | xargs -n1 basename | sed 's/_R1\.fq\.gz$//' | sort > "${SAMPLES}"

# check every R1 has an R2
while read -r s; do
    [[ -f "${RAW_DIR}/${s}_R2.fq.gz" ]] || { echo "ERROR: ${s}_R2.fq.gz missing"; exit 1; }
done < "${SAMPLES}"

# reference copies with fixed names
for f in "${REF_FASTA_SRC}" "${REF_GFF_SRC}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: reference file not found: $f"
        echo "Run 00b_get_reference.sh first (see its header for the JGI or NCBI option)."
        exit 1
    fi
done
copy_maybe_gz() { if [[ "$1" == *.gz ]]; then gunzip -c "$1" > "$2"; else cp "$1" "$2"; fi; }
copy_maybe_gz "${REF_FASTA_SRC}" "${REF}"
copy_maybe_gz "${REF_GFF_SRC}"   "${REF_DIR}/Ss_ref.gff"
if [[ -n "${REF_PROT_SRC}" && -f "${REF_PROT_SRC}" ]]; then
    copy_maybe_gz "${REF_PROT_SRC}" "${REF_DIR}/Ss_ref_proteins.fa"
fi
echo "Reference : $(grep -c '>' "${REF}") contigs, $(du -h "${REF}" | cut -f1)"
echo "Gene GFF  : $(grep -vc '^#' "${REF_DIR}/Ss_ref.gff") lines"

N=$(wc -l < "${SAMPLES}")
echo "Isolates: ${N}"
cat "${SAMPLES}"
echo
echo "For array jobs use:  --array=1-${N}"
