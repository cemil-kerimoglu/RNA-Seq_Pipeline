# Commands at a glance

```bash
# First-time setup on the Google Compute Engine VM
sudo apt-get update
sudo apt-get install -y git curl bzip2 tmux make
git clone https://github.com/cemil-kerimoglu/RNA-Seq_Pipeline.git
cd RNA-Seq_Pipeline
bash bootstrap/setup_gcp_vm.sh
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate setd1b-rnaseq-workflow

# Validate the DAG
make gcp-dryrun

# Build rule environments and run in a persistent terminal
tmux new -s setd1b
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate setd1b-rnaseq-workflow
cd "$HOME/RNA-Seq_Pipeline"
make gcp-envs
make gcp-run
# Detach: Ctrl-b, then d
# Reattach later: tmux attach -t setd1b

# Monitor from a second SSH terminal
free -h
df -h /
du -sh data work results 2>/dev/null

# Package final results after successful completion
make archive-results
# Download: exports/setd1b-results.tar.gz
# With Google Cloud CLI on the laptop:
gcloud compute scp rnaseq-vm:~/RNA-Seq_Pipeline/exports/setd1b-results.tar.gz . --zone=europe-west10-a

# Resume after interruption
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate setd1b-rnaseq-workflow
cd "$HOME/RNA-Seq_Pipeline"
make gcp-run

# Slurm alternative (after adapting that profile to the cluster)
snakemake --snakefile workflow/Snakefile --profile profiles/slurm

# Inspect final QC first
# results/multiqc/multiqc_report.html

# Main DE table
# results/deseq2/deseq2_all.tsv

# Paper-like downregulated set
# results/deseq2/deseq2_paper_down.tsv

# GO enrichment emphasized by the paper
# results/topgo/paper_down_BP.tsv
```
