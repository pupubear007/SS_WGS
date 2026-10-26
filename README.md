# SS_WGS — *Sclerotinia sclerotiorum* whole-genome analysis (44 isolates, MSI)

A SLURM pipeline for whole-genome sequencing of 44 *S. sclerotiorum* isolates. It covers read QC,
reference-based variant calling and population structure, de novo assembly and annotation, the
pangenome, and gene-level variation keyed to the JGI `SS1G_` gene IDs used in the soybean /
sunflower dual RNA-seq ([UG100_RNAseq](https://github.com/pupubear007/UG100_RNAseq)).

## Contents

- [Study design](#study-design)
- [Pipeline](#pipeline)
- [Requirements](#requirements)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)
- [Repository layout](#repository-layout)
- [Notes](#notes)
- [License](#license)
- [Citation](#citation)

## Study design

| Item | Value |
|---|---|
| Organism | *Sclerotinia sclerotiorum* (haploid) |
| Isolates | 44, paired-end Illumina WGS |
| Reference | *S. sclerotiorum* 1980 UF-70, JGI MycoCosm Sclsc1 (same as the RNA-seq) |
| Gene IDs | JGI `SS1G_` |
| Compute | MSI (University of Minnesota), SLURM |

## Pipeline

| Step | Script | Job type | What it does |
|---|---|---|---|
| 00a | `00a_install_miniforge.sh` | bash, login node, once | own conda + mamba in `~/miniforge3` |
| 00 | `00_setup_env.sh` | sbatch, once | conda envs + funannotate databases |
| 00b | `00b_get_reference.sh` | bash, once | download the reference |
| 01 | `01_make_sample_list.sh` | bash | `samples.txt`; prints N (number of isolates) |
| 02 | `02_qc_trim.sh` | sbatch array | FastQC + fastp |
| 03 | `03_qc_summary.sh` | sbatch | coverage table + MultiQC |
| 10 | `10_ref_prep.sh` | sbatch | index reference |
| 11 | `11_align_call.sh` | sbatch array | bwa-mem2 + markdup + HaplotypeCaller (ploidy 1) |
| 12 | `12_joint_genotype.sh` | sbatch | joint genotyping + hard filters |
| 13 | `13_tree_pca.sh` | sbatch | IQ-TREE phylogeny, PCA, SNP distances |
| 14 | `14_gene_level.sh` | sbatch | SnpEff effects, gene deletions, gene function (`SS1G_` IDs) |
| 20 | `20_assemble.sh` | sbatch array | SPAdes |
| 20b | `20b_decontaminate.sh` | sbatch | remove non-*Sclerotinia* contigs, name the contaminant |
| 21 | `21_assembly_qc.sh` | sbatch | QUAST + BUSCO |
| 22 | `22_annotate.sh` | sbatch array | funannotate (uses the cleaned assemblies) |
| 23 | `23_orthofinder.sh` | sbatch | pangenome: core / accessory gene families |
| 30 | `30_collect_results.sh` | bash, login node | gather results into one folder |

Steps 0x are setup, 1x reference-based, 2x de novo, 30 results. The 1x and 2x branches can run in parallel.

## Requirements

- An HPC cluster with SLURM (written for MSI, `agate`)
- Miniforge / conda; `00_setup_env.sh` creates the environments:

| Env (in `config.sh`) | Tools |
|---|---|
| `pu_WGS` | fastqc, multiqc, fastp, seqkit, pigz, bwa-mem2, samtools, bcftools, gatk4, snpeff, iqtree, python 3.11, numpy |
| `pu_WGS_asm` | spades, rasusa, seqkit, quast, busco, orthofinder |
| `pu_quast` | quast on python 3.11 |
| `pu_funannotate` | funannotate + its databases |

- For the report: R and RStudio with ggplot2, ape, phangorn, pheatmap, ggrepel, dplyr, tidyr, knitr, scales, vegan (the Rmd installs any that are missing)

## Configuration

All paths are set in `config.sh`; no script needs path edits.

| Setting | Current value |
|---|---|
| `RAW_DIR` | `/scratch.global/wan00965/Quick_test/SS_WGS_rawfastq_merged` (`<ISOLATE>_R1.fq.gz` / `_R2.fq.gz`) |
| `PROJECT` | `/scratch.global/wan00965/Quick_test/SS_WGS_analysis` |
| `REF_STORE` | `/scratch.global/wan00965/reference/Ss_1980_UF70` |
| `BUSCO_LINEAGE` | `leotiomycetes_odb10` |

The `#SBATCH` lines in each script use account `mmccaghe`, partition `msismall`, and mail to
`wan00965@umn.edu`. Change these in the script headers if you run under another account.

The raw-file scripts in `raw_file_manage/` set their own paths at the top of each script.

## Usage

Run everything from the repo root, so every job finds `config.sh` and writes its logs to `log/`.

```bash
mkdir -p log

# --- raw FASTQ (once, before the pipeline) ---
sbatch raw_file_manage/gather_raw_fastq.sh       # MSI: collect batch folders + md5 manifest
sbatch raw_file_manage/merge_rename_fastq.sh     # MSI: merge lanes, rename to <ISOLATE>_R1/_R2
bash   raw_file_manage/download_raw_fastq.sh     # Mac: optional local copy to the SSD

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

## Output

Under `PROJECT` on MSI:

| Folder | Contents |
|---|---|
| `01_qc/`, `02_trimmed/` | FastQC, fastp reports, trimmed reads |
| `03_align/`, `04_gvcf/`, `05_vcf/` | BAMs, per-isolate gVCFs, joint filtered VCF |
| `06_popgen/` | phylogeny, PCA, SNP distance matrix |
| `07_assembly/`, `08_assembly_qc/` | SPAdes assemblies (raw + cleaned), QUAST / BUSCO |
| `09_annotation/`, `10_orthofinder/` | funannotate gene models, pangenome gene families |
| `SS_WGS_results/` | everything the report needs, from step 30 (with `MANIFEST.tsv`) |

The report (`Rmd/SS_WGS_results.Rmd`) covers clonal lineages, the RNA-seq isolates in genomic
context, Xtra7 vs SSPotter, gene deletions, protein-changing variation, the link to RNA-seq DE
genes, and phenotype vs genotype. Tables go to `results/report_tables/`, figures to `results/figures/`.

## Repository layout

| Path | Contents |
|---|---|
| `00_…`–`30_…`, `config.sh` | the pipeline (MSI) |
| `raw_file_manage/` | gather / merge / download the raw FASTQ (run before the pipeline) |
| `Rmd/SS_WGS_results.Rmd` | the results report |
| `tools/` | optional utilities (`31_inventory.sh`: list every file the pipeline produced, as a SLURM job) |
| `log/` | SLURM `.out` / `.err` files (contents not tracked) |
| `rnaseq_DE/` | pathogen DE lists from the dual RNA-seq, written by `~/RNAseq_paper/CH5_DE_tables/gather_DE_tables.py` (not tracked in git) |
| `results/` | downloaded results (not tracked in git) |

## Notes

- *S. sclerotiorum* is haploid → GATK `-ploidy 1`.
- 10 isolates (library series CKDN2600112xx) carried bacterial DNA (*Acinetobacter*, *E. coli*, *Pseudomonas*);
  `20b` removes it before annotation. Variant calling is unaffected (bacterial reads do not map).
- Genes with low read depth in Model_1980 (the reference strain) are treated as unmappable, not deleted.
- funannotate runs without GeneMark; secreted-protein / KOG / InterPro labels come from the JGI annotation files (step 14).

## License

The code in this repository is released under the [MIT License](LICENSE).
Sequencing reads and analysis results are not included in the repository and are not covered by this license.

## Citation

See [`CITATION.cff`](CITATION.cff), or use **Cite this repository** on the GitHub page.
