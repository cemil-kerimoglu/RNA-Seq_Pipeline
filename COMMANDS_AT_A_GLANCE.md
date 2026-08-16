# Commands at a glance

```bash
# 1. Create workflow environment
mamba env create -f bootstrap/environment.yaml
conda activate setd1b-rnaseq-workflow

# 2. Inspect the plan
snakemake --snakefile workflow/Snakefile --dry-run --printshellcmds

# 3. Run on Slurm (after adapting the profile)
snakemake --snakefile workflow/Snakefile --profile profiles/slurm

# 4. Inspect final QC first
# results/multiqc/multiqc_report.html

# 5. Main DE table
# results/deseq2/deseq2_all.tsv

# 6. Paper-like downregulated set
# results/deseq2/deseq2_paper_down.tsv

# 7. GO enrichment emphasized by the paper
# results/topgo/paper_down_BP.tsv
```
