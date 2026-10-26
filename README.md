# SS_WGS — *Sclerotinia sclerotiorum* whole-genome analysis (44 isolates, MSI)

Whole-genome sequencing of 44 *S. sclerotiorum* isolates: read QC, reference-based variant calling
and population structure, de novo assembly and annotation, pangenome, and gene-level variation keyed
to the JGI `SS1G_` gene IDs used in the soybean / sunflower dual RNA-seq.

All steps are SLURM jobs run from this folder; logs go to `log/`. Paths are set only in `config.sh`.
Reference: *S. sclerotiorum* 1980 UF-70, JGI MycoCosm Sclsc1 (same as the RNA-seq).

## Run order

```bash
mkdir -p log
# --- one-time setup ---
bash   00a_install_miniforge.sh        # own conda + mamba in ~/miniforge3
sbatch 00_setup_env.sh                 # conda envs + funannotate databases
bash   00b_get_reference.sh jgi        # reference -> /scratch.global/wan00965/reference/Ss_1980_UF70
bash   01_make_sample_list.sh          # samples.txt; prints N (number of isolates)

# --- QC ---
sbatch --array=1-N 02_qc_trim.sh       # FastQC + fastp
sbatch 03_qc_summary.sh                # coverage table + MultiQC

# --- reference-based ---
sbatch 10_ref_prep.sh                  # index reference
sbatch --array=1-N 11_align_call.sh    # bwa-mem2 + markdup + HaplotypeCaller (ploidy 1)
sbatch 12_joint_genotype.sh            # joint genotyping + hard filters
sbatch 13_tree_pca.sh                  # IQ-TREE phylogeny, PCA, SNP distances
sbatch 14_gene_level.sh                # SnpEff effects, gene deletions, gene function (SS1G_ IDs)

# --- de novo --- (can run in parallel with the reference-based steps)
sbatch --array=1-N 20_assemble.sh      # SPAdes
sbatch 20b_decontaminate.sh            # remove non-Sclerotinia contigs, name the contaminant
sbatch --time=48:00:00 21_assembly_qc.sh   # QUAST + BUSCO
sbatch --array=1-N 22_annotate.sh      # funannotate (uses the cleaned assemblies)
sbatch 23_orthofinder.sh               # pangenome: core / accessory gene families

# --- results ---
bash   30_collect_results.sh           # gather results -> SS_WGS_analysis/SS_WGS_results/
```

Then on the Mac:

```bash
rsync -av wan00965@agate.msi.umn.edu:/scratch.global/wan00965/Quick_test/SS_WGS_analysis/SS_WGS_results/ ~/SS_WGS/results/
```

and knit `Rmd/SS_WGS_results.Rmd` in RStudio.

## Folder layout

| Path | Contents |
|---|---|
| `00_…`–`30_…`, `config.sh` | the pipeline (MSI) |
| `raw_file_manage/` | gather / merge / download the raw FASTQ (run before the pipeline) |
| `Rmd/SS_WGS_results.Rmd` | report: clonal lineages, RNA-seq isolates in genomic context, Xtra7 vs SSPotter, gene deletions, protein-changing variation, link to RNA-seq DE genes, phenotype vs genotype |
| `rnaseq_DE/` | pathogen DE lists from the dual RNA-seq, written by `~/RNAseq_paper/CH5_DE_tables/gather_DE_tables.py` (not tracked in git) |
| `results/` | downloaded results (not tracked in git) |
| `tools/` | optional utilities (`31_inventory.sh`: list every file the pipeline produced, as a SLURM job) |

## Notes
- *S. sclerotiorum* is haploid → GATK `-ploidy 1`.
- 10 isolates (library series CKDN2600112xx) carried bacterial DNA (*Acinetobacter*, *E. coli*, *Pseudomonas*);
  `20b` removes it before annotation. Variant calling is unaffected (bacterial reads do not map).
- Genes with low read depth in Model_1980 (the reference strain) are treated as unmappable, not deleted.
- funannotate runs without GeneMark; secreted-protein / KOG / InterPro labels come from the JGI annotation files (step 14).

## License

The code in this repository is released under the [MIT License](LICENSE).
Sequencing reads and analysis results are not included in the repository and are not covered by this license.
To cite this pipeline, see [`CITATION.cff`](CITATION.cff).
