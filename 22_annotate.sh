#!/bin/bash
#SBATCH --job-name=22_annotate
#SBATCH --account=mmccaghe
#SBATCH --partition=msismall
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

# ============================================================
# Step 22 (array, one task per isolate, after step 21 looks good):
#   funannotate clean -> sort -> repeat mask -> gene prediction
#   -> functional annotation (Pfam, UniProt, MEROPS, CAZy, BUSCO)
#   sbatch --array=1-N 22_annotate.sh
# Protein evidence: reference 1980 UF-70 proteins + UniProt/SwissProt.
# Output: 09_annotation/<ISOLATE>/annotate_results/*.gff3, *.proteins.fa
# ============================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"; source ./config.sh
activate_env "${ENV_ANNO}"
export FUNANNOTATE_DB

S=$(get_sample)
ASM="${ASM_DIR}/clean/${S}.scaffolds.fa"        # decontaminated (step 20b)
[[ -s "${ASM}" ]] || { echo "WARNING: no clean assembly for ${S}; using raw (run 20b first)"; ASM="${ASM_DIR}/${S}.scaffolds.fa"; }
OUT="${ANNO_DIR}/${S}"
mkdir -p "${OUT}"
cd "${OUT}"
echo "Isolate: ${S}   $(date)"

# contig names must be <= 16 characters for GenBank-style output
SHORT=$(echo "${S}" | tr -cd '[:alnum:]' | cut -c1-10)

# 1. remove duplicate small contigs, sort/rename, soft-mask repeats
funannotate clean -i "${ASM}" -o "${S}.clean.fa" --minlen 1000 --cpus "${THREADS}"
funannotate sort  -i "${S}.clean.fa" -o "${S}.sorted.fa" -b "${SHORT}" --minlen 1000
funannotate mask  -i "${S}.sorted.fa" -o "${S}.masked.fa" --cpus "${THREADS}"

# 2. gene prediction
#    Botrytis cinerea (closest Augustus model) seeds training.
funannotate predict \
    -i "${S}.masked.fa" -o "${OUT}/fun" \
    -s "Sclerotinia sclerotiorum" --strain "${S}" \
    --busco_seed_species botrytis_cinerea \
    --protein_evidence "${REF_DIR}/Ss_ref_proteins.fa" "${FUNANNOTATE_DB}/uniprot_sprot.fasta" \
    --cpus "${THREADS}"

# 3. functional annotation (Pfam, UniProt, MEROPS, CAZymes, GO via eggNOG if installed)
funannotate annotate -i "${OUT}/fun" --cpus "${THREADS}"

ls "${OUT}/fun/annotate_results/"
echo "Done ${S}   $(date)"
