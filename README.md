# SS_WGS pipeline — *Sclerotinia sclerotiorum* whole-genome analysis (MSI)

All steps are SLURM jobs. Run every command from this folder; logs go to `log/`.
Edit paths only in `config.sh`. The reference is kept in `/scratch.global/wan00965/reference/Ss_1980_UF70`.

## Run order

```bash
mkdir -p log
bash   00a_install_miniforge.sh        # once: your own current conda + mamba
sbatch 00_setup_env.sh                 # once: conda envs + funannotate DB
bash   00b_get_reference.sh jgi        # once: reference -> /scratch.global/wan00965/reference (or: ... ncbi)
bash   01_make_sample_list.sh          # prints N = number of isolates

# --- QC ---
sbatch --array=1-N 02_qc_trim.sh       # FastQC + fastp trimming
sbatch 03_qc_summary.sh                # coverage table + MultiQC

# --- Reference-based ---
sbatch 10_ref_prep.sh                  # index reference + SnpEff DB
sbatch --array=1-N 11_align_call.sh    # bwa-mem2 + markdup + HaplotypeCaller (ploidy 1)
sbatch 12_joint_genotype.sh            # joint genotyping, filters, SnpEff
sbatch 13_tree_pca.sh                  # IQ-TREE phylogeny, PCA, SNP distances

# --- De novo ---  (can run in parallel with the reference-based steps)
sbatch --array=1-N 20_assemble.sh      # SPAdes per isolate
sbatch 21_assembly_qc.sh               # QUAST + BUSCO
sbatch --array=1-N 22_annotate.sh      # funannotate gene models + function
sbatch 23_orthofinder.sh               # core / accessory genes across isolates
```

Wait for each step to finish before the next one in its branch. To chain
automatically: `sbatch --dependency=afterok:<JOBID> next_step.sh`.

## Where to look

| Question | File |
|---|---|
| Read quality, coverage per isolate | `01_qc/multiqc_report.html`, `01_qc/seqkit/coverage_per_isolate.tsv` |
| Mapping rate / depth per isolate | `03_align/*.flagstat.txt`, `*.coverage_summary.txt` |
| Variants with gene effects | `05_vcf/SS_WGS.PASS.snpeff.vcf.gz`, `SS_WGS.high_moderate_effects.tsv` |
| Relationships among isolates | `06_popgen/SS_WGS_snps.treefile`, `pca.tsv`, `snp_distance.tsv` |
| Assembly quality | `08_assembly_qc/quast/report.html`, `busco_summary.tsv` |
| Gene content differences | `10_orthofinder/.../Orthogroups.GeneCount.tsv` |

## Notes
- *S. sclerotiorum* is haploid → GATK `-ploidy 1`.
- Assemblies are downsampled to 100x before SPAdes.
- funannotate runs without GeneMark (licensed); add it later if wanted.
- Effector prediction (SignalP + EffectorP) and antiSMASH can be run on the
  `22_annotate` protein files as a follow-up.
