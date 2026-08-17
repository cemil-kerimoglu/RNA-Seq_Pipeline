# SETD1B neuronal nuclear RNA-seq — modern reproducible reanalysis scaffold

This repository is an **inspectable, end-to-end Snakemake workflow** for the six main neuronal nuclear RNA-seq libraries from the SETD1B study (`GSE180326`). It starts from SRA accessions, creates FASTQ files, performs QC and trimming, maps to the mouse genome, generates two alternative count matrices, performs DESeq2 differential expression, and runs topGO enrichment.

It is designed as a **modernized reanalysis**, not as a claim of byte-for-byte reproduction of the 2022 analysis. The paper reports STAR + featureCounts + DESeq2 and says that spliced and unspliced transcripts were counted, but it does not report the exact GTF release or the exact featureCounts command-line switches. Those missing details are therefore made explicit and configurable here.

## 1. Samples

The workflow uses only the six neuronal RNA-seq samples underlying the main control-vs-Setd1b-cKO comparison, not the whole-CA libraries and not the single-nucleus RNA-seq experiment.

| Sample               | GEO        | SRA experiment | Run         | Condition |
| -------------------- | ---------- | -------------- | ----------- | --------- |
| Neuron_RNA_cKO_1     | GSM5460849 | SRX11490981    | SRR15183990 | cKO       |
| Neuron_RNA_cKO_2     | GSM5460850 | SRX11490982    | SRR15183991 | cKO       |
| Neuron_RNA_cKO_3     | GSM5460851 | SRX11490983    | SRR15183992 | cKO       |
| Neuron_RNA_Control_1 | GSM5460852 | SRX11490984    | SRR15183993 | control   |
| Neuron_RNA_Control_2 | GSM5460853 | SRX11490985    | SRR15183994 | control   |
| Neuron_RNA_Control_3 | GSM5460854 | SRX11490986    | SRR15183995 | control   |

All six public runs are single-end Illumina RNA-seq libraries. The paper describes the principal neuronal libraries as 75-bp single-end NextSeq 550 data.

## 2. Pipeline overview

```text
NCBI SRA
  |
  | prefetch -> vdb-validate -> fasterq-dump -> pigz
  v
FASTQ.gz
  |
  |-- FastQC (raw)
  |
  |-- fastp (adapter/quality filtering + explicit poly-G trimming)
  |       |
  |       `-- FastQC (post-trim)
  v
trimmed FASTQ.gz
  |
  | STAR 2-pass -> coordinate-sorted BAM -> samtools index/QC
  v
BAMs
  |
  |-- featureCounts on exons --------------------> exon count matrix
  |
  `-- featureCounts on whole gene bodies -------> nuclear-RNA count matrix [default]
                                                      |
                                                      v
                                                   DESeq2
                                                      |
                       +------------------------------+------------------+
                       |                              |                  |
                       v                              v                  v
                      PCA                       DE tables          MA/volcano
                                                      |
                                                      v
                                                    topGO

FastQC + fastp + STAR + samtools + featureCounts metrics -> MultiQC
Raw FASTQs + reference FASTA/GTF -> SHA256 provenance manifest
```

## 3. Why two featureCounts modes?

The paper states that the bulk **nuclear** RNA-seq reads were counted while considering both spliced and unspliced transcripts. Nuclear RNA contains substantial intronic signal from pre-mRNA.

Because the exact original featureCounts flags/annotation were not published, the workflow generates both:

1. **`exon_counts.tsv`** — conventional gene-level exon counting (`-t exon -g gene_id`).
2. **`gene_body_counts.tsv`** — a deliberately transparent nuclear-RNA approximation in which each annotated gene's complete genomic span is converted to SAF and counted. This allows intronic and exonic reads to contribute.

`config/config.yaml` defaults to `count_mode: gene_body` to stay conceptually close to the paper. Change it to `exon` to perform a conventional exonic bulk RNA-seq analysis. The two matrices are always generated so they can be compared.

**Important caveat:** whole-gene-body counting is not asserted to be the undocumented command used in the original paper. It is a reproducible reconstruction of the stated goal of including unspliced nuclear transcripts.

## 4. Reference choice

The workflow deliberately stays on **GRCm38/mm10** for comparability with the paper and pins an archival **GENCODE M25 / GRCm38 primary assembly** FASTA and GTF. This avoids silently switching the biological coordinate system to GRCm39 in a study-reproduction exercise.

The exact original annotation release was not specified in the paper, so M25 is a documented project choice, not a reconstructed historical fact.

For 75-bp reads, the STAR index uses `sjdbOverhang = 74` (`read length - 1`).

## 5. Modernizations relative to the original analysis

The biological/statistical backbone remains intentionally recognizable: **STAR -> featureCounts -> DESeq2 -> topGO**. Modern engineering and QC are added around it:

- Snakemake DAG instead of hand-run commands.
- Per-rule Conda environments and pinned major tool versions.
- A single-VM Google Compute Engine profile sized for 8 vCPUs, 64 GB RAM,
  and a 500 GB boot disk.
- A Slurm executor profile for HPC execution.
- `prefetch`, `vdb-validate`, and `fasterq-dump` for traceable SRA retrieval.
- FastQC both before and after filtering.
- fastp for adapter/quality handling; poly-G trimming is explicitly enabled because these public runs were produced on a NextSeq instrument.
- STAR two-pass alignment.
- samtools `flagstat` and `stats` for alignment QC.
- MultiQC as a single QC dashboard.
- No PCR-duplicate removal from ordinary RNA-seq BAMs: without UMIs, identical RNA-seq fragments can represent genuine high expression and should not automatically be treated as PCR artifacts.
- DESeq2 is run on **raw integer counts**, never RPKM/TPM.
- `apeglm` shrinkage is added for more stable effect-size ranking and plotting. P-values/padj still come from the original DESeq2 Wald test.
- Both the paper-like threshold (`padj < 0.10` and `|FC| > 1.2`) and a conventional `FDR < 0.05` reporting set are written.
- topGO results include `classic`, `weight`, and `weight01` Fisher tests. The paper reports a weighted topGO analysis with weighted P < 0.005 but does not unambiguously name the exact topGO algorithm. The output therefore provides both `paper_like_weight_pass` and `paper_like_weight01_pass` at that numeric cutoff; current topGO documents `weight01` as its default algorithm.
- The default GO universe is the genes that actually entered DESeq2 after requiring at least 10 counts in at least the smallest group size (`tested`), avoiding enrichment bias from genes that could never have been selected. Set `universe: all_annotated` to make the universe closer to the paper's stated all-genome reference.
- SHA256 checksums and R `sessionInfo()` are preserved for provenance.

## 6. Folder structure

```text
RNA-Seq_Pipeline/
├── README.md
├── Makefile
├── bootstrap/
│   ├── environment.yaml
│   └── setup_gcp_vm.sh
├── config/
│   ├── config.yaml
│   └── samples.tsv
├── profiles/
│   ├── gcp-vm/
│   │   └── profile.v9+.yaml
│   └── slurm/
│       └── profile.v9+.yaml
├── resources/
│   └── multiqc_config.yaml
├── workflow/
│   ├── Snakefile
│   ├── envs/
│   │   ├── alignment.yaml
│   │   ├── counts.yaml
│   │   ├── download.yaml
│   │   ├── qc.yaml
│   │   ├── reference.yaml
│   │   └── rnaseq_r.yaml
│   └── scripts/
│       ├── build_gene_body_saf.py
│       ├── clean_featurecounts.py
│       ├── deseq2.R
│       └── topgo.R
├── data/raw/       # generated; gitignored
├── work/           # generated; gitignored
├── results/        # generated; gitignored
└── logs/           # generated; gitignored
```

## 7. Bootstrap the workflow environment

The workflow requires Linux. On a system that already has Conda or Mamba:

```bash
cd RNA-Seq_Pipeline
mamba env create -f bootstrap/environment.yaml
conda activate setd1b-rnaseq-workflow
```

Snakemake will create the smaller per-rule environments automatically when `--software-deployment-method conda` is enabled.

The environment files use `conda-forge` and `bioconda`, with `nodefaults` to
prevent an existing user-level `defaults` channel from being mixed into the
solve. The Makefile applies strict channel priority only to each workflow
command. This does **not** remove Conda, rewrite the user's global Conda
configuration, or alter other environments.

The downstream R environment is pinned as one internally consistent
Bioconductor 3.22 family: R 4.5, DESeq2 1.50.2, apeglm 1.32.0, topGO 2.62.0,
org.Mm.eg.db 3.22.0, and AnnotationDbi 1.72.0. Bioconductor releases are tied
to specific R minor versions, so these pins must be advanced together rather
than individually.

## 8. Run on the Google Compute Engine VM

The included `gcp-vm` profile is tuned to this project's current VM: 8 vCPUs,
about 62 GiB usable RAM, and about 484 GiB usable disk. It runs Snakemake's
local executor on the remote Linux machine; Google Compute Engine provides the
machine, while Snakemake schedules the individual rules on it.

### First-time VM preparation

Open an SSH terminal to the VM in Google Cloud Console and run:

```bash
sudo apt-get update
sudo apt-get install -y git curl bzip2 tmux make
git clone https://github.com/cemil-kerimoglu/RNA-Seq_Pipeline.git
cd RNA-Seq_Pipeline
bash bootstrap/setup_gcp_vm.sh
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate setd1b-rnaseq-workflow
```

If the repository was already cloned, use `cd` to enter it and run `git pull`
instead of cloning it again. The setup script installs Miniforge only inside
the current Linux user's home directory, applies strict channel priority only
while solving this environment, and creates or updates the workflow's bootstrap
environment. It is safe to rerun.

### Validate before spending compute time

```bash
make gcp-dryrun
```

This constructs and validates the complete DAG without executing it.

### Run persistently inside tmux

Browser SSH sessions can disconnect. Start the workflow inside `tmux`:

```bash
tmux new -s setd1b
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate setd1b-rnaseq-workflow
cd "$HOME/RNA-Seq_Pipeline"
make gcp-envs
make gcp-run
```

`make gcp-envs` pre-creates all rule-specific Conda environments, so dependency
problems are discovered before large downloads and alignment start. It is kept
inside `tmux` because the initial dependency solve can take some time.

Detach without stopping the workflow by pressing `Ctrl-b`, releasing the keys,
and then pressing `d`. Later, reconnect over SSH and run:

```bash
tmux attach -t setd1b
```

Snakemake uses `rerun-incomplete: true`, so after an interrupted VM or process
it is normally sufficient to activate the environment, return to the repository,
and run `make gcp-run` again. Completed outputs are retained.

### Monitor, package, and download the results

In a second SSH terminal, useful checks are:

```bash
free -h
df -h /
du -sh data work results 2>/dev/null
```

After successful completion:

```bash
make archive-results
ls -lh exports/setd1b-results.tar.gz
```

Download `exports/setd1b-results.tar.gz` through the Compute Engine SSH file
transfer menu. Alternatively, from a laptop terminal with Google Cloud CLI
installed, run:

```bash
gcloud compute scp \
  rnaseq-vm:~/RNA-Seq_Pipeline/exports/setd1b-results.tar.gz . \
  --zone=europe-west10-a
```

The raw reads and large intermediate BAM/index files are deliberately not
included in the archive.

Finally, stop the VM from the Google Cloud Console when it is not in use.
Stopping ends vCPU/RAM charges, although the persistent boot disk continues to
incur storage charges until deleted. Do not delete the VM or disk until the
results have been downloaded and checked.

## 9. Inspect before running

Always start with a dry run:

```bash
snakemake --snakefile workflow/Snakefile \
  --dry-run \
  --printshellcmds
```

or:

```bash
make dryrun
```

You can also render the DAG on a machine with Graphviz:

```bash
snakemake --snakefile workflow/Snakefile --dag | dot -Tpdf > dag.pdf
```

## 10. Run on Slurm

Edit `profiles/slurm/profile.v9+.yaml` first. In particular, replace the generic partition name `compute` and adapt memory/runtime defaults to the local cluster.

Then:

```bash
snakemake --snakefile workflow/Snakefile \
  --profile profiles/slurm
```

The profile uses the Snakemake Slurm executor plugin. Cluster credentials, modules, account/QoS rules and internet access are site-specific and intentionally not hard-coded.

For a small local test instead:

```bash
snakemake --snakefile workflow/Snakefile \
  --cores 8 \
  --software-deployment-method conda \
  --printshellcmds
```

## 11. Key outputs

### QC

```text
results/multiqc/multiqc_report.html
```

Inspect this **before** interpreting differential expression. At minimum, check raw/post-trim sequence quality, read loss during fastp, mapping rate, multimapping, featureCounts assignment rate, and whether any replicate is an obvious outlier.

### Count matrices

```text
results/counts/exon_counts.tsv
results/counts/gene_body_counts.tsv
```

The matrix chosen by `quantification.count_mode` is copied into `work/selected_counts.tsv` and used downstream.

### Differential expression

```text
results/deseq2/deseq2_all.tsv
results/deseq2/deseq2_paper_significant.tsv
results/deseq2/deseq2_paper_up.tsv
results/deseq2/deseq2_paper_down.tsv
results/deseq2/deseq2_fdr_significant.tsv
results/deseq2/deseq2_fdr_up.tsv
results/deseq2/deseq2_fdr_down.tsv
results/deseq2/normalized_counts.tsv
results/deseq2/vst_counts.tsv
results/deseq2/PCA.pdf
results/deseq2/sample_distance_heatmap.pdf
results/deseq2/MA_apeglm.pdf
results/deseq2/volcano.pdf
results/deseq2/dds.rds
results/deseq2/vsd.rds
results/deseq2/sessionInfo.txt
```

`deseq2_all.tsv` contains both the unshrunk DESeq2 log2 fold change used for the paper-like threshold and an `apeglm`-shrunken log2 fold change for ranking/visualization.

### GO enrichment

```text
results/topgo/paper_down_BP.tsv
results/topgo/paper_up_BP.tsv
results/topgo/fdr_down_BP.tsv
results/topgo/fdr_up_BP.tsv
```

For reproducing the biological question emphasized in the paper, start with `paper_down_BP.tsv`.

### Provenance

```text
results/provenance/sha256sums.txt
results/deseq2/sessionInfo.txt
```

## 12. Important analysis decisions to review before real use

### Strandedness

The public metadata and paper do not give a sufficiently explicit library-strand setting to infer a featureCounts `-s` value with complete confidence. The project therefore defaults conservatively to **unstranded (`-s 0`)**. On a real reanalysis, verify strandedness empirically (for example with RSeQC `infer_experiment.py`) and change `quantification.strandedness` if needed.

### Biological covariates

The public sample sheet used here contains only the condition and replicate labels needed for the published comparison. Do **not** invent batch, sex, age, lane, or other covariates. If reliable sample-level metadata are recovered, add them to `samples.tsv` and update the DESeq2 design formula accordingly (for example `~ batch + condition`).

### Read trimming

Modern RNA-seq pipelines increasingly avoid aggressive trimming. This workflow performs conservative fastp filtering and retains post-trim FastQC so any effect is visible. If raw data are already clean, the filtering settings can be relaxed.

### Duplicate reads

The workflow intentionally **does not deduplicate RNA-seq BAM files**. Duplicate removal was appropriate in the paper's ChIP-seq workflow, but it is generally not a default step for ordinary non-UMI RNA-seq.

### RPKM/TPM

The original paper calculated RPKM for descriptive expression plots. This workflow does not use RPKM for differential testing. DESeq2 receives raw counts. If TPM/RPKM values are wanted for visualization or cross-gene descriptive purposes, add them as a separate derived output rather than feeding them to DESeq2.

## 13. Paper-like versus modern GO universes

`config/config.yaml` defaults to:

```yaml
gene_ontology:
  universe: tested
```

This asks: among genes that had enough information to enter this DESeq2 analysis and that have GO mappings, are selected genes unusually enriched for a GO term?

To more closely mimic the Methods statement that all genomic genes were the reference set:

```yaml
gene_ontology:
  universe: all_annotated
```

Both choices are explicit because gene-universe definition can materially affect enrichment results.

## 14. What this scaffold deliberately does not pretend to know

A high-quality reanalysis should distinguish **known provenance** from **reasonable reconstruction**. The following are not claimed to be known from the paper:

- the exact historical GTF release used in 2020/2021;
- the exact featureCounts options used to include spliced and unspliced nuclear transcripts;
- an experimentally confirmed strandedness value for these public FASTQ files;
- hidden batch variables not represented in the public metadata.

Those choices are surfaced in configuration rather than silently guessed.

## 15. If I were deploying this for production

For a real institute/core-facility analysis I would additionally consider:

- institutional container images (Apptainer/Singularity) instead of solving Conda environments on every cluster;
- automated strandedness inference;
- reference checksums from the reference provider;
- a CI test with tiny synthetic FASTQs;
- workflow linting and a frozen lock/container digest;
- a project-specific MultiQC acceptance checklist;
- export of a complete software/environment manifest with the final report.

The present repository is intentionally readable enough to study line-by-line while still following an HPC workflow-engine pattern.

## 16. Why a custom Snakemake workflow instead of simply calling nf-core/rnaseq?

As of August 2026, nf-core/rnaseq is an excellent turnkey choice for routine production RNA-seq and its current stable branch provides several alignment/quantification modes plus extensive QC. For this study, however, a small custom workflow is useful pedagogically and scientifically because you asked for code that can be inspected line-by-line and because the original experiment is **nuclear RNA-seq with an explicit requirement to include unspliced signal**. The custom `gene_body` featureCounts branch makes that project-specific decision visible rather than hiding it behind a generic transcript-quantification default.

For a new conventional bulk RNA-seq project with no historical-analysis constraint, I would seriously consider starting from nf-core/rnaseq and adding only the project-specific downstream statistics. For this SETD1B reconstruction, keeping STAR + featureCounts + DESeq2 also makes the relationship to the published analysis much easier to explain in an interview.
