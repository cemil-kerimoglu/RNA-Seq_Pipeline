SHELL := /bin/bash

.PHONY: dryrun local gcp-dryrun gcp-envs gcp-run archive-results slurm clean

dryrun:
	snakemake --snakefile workflow/Snakefile --dry-run --printshellcmds

local:
	CONDA_CHANNEL_PRIORITY=strict snakemake --snakefile workflow/Snakefile --cores 8 --software-deployment-method conda --printshellcmds

gcp-dryrun:
	CONDA_CHANNEL_PRIORITY=strict snakemake --snakefile workflow/Snakefile --profile profiles/gcp-vm --dry-run

gcp-envs:
	CONDA_CHANNEL_PRIORITY=strict snakemake --snakefile workflow/Snakefile --profile profiles/gcp-vm --conda-create-envs-only

gcp-run:
	CONDA_CHANNEL_PRIORITY=strict snakemake --snakefile workflow/Snakefile --profile profiles/gcp-vm

archive-results:
	test -d results
	mkdir -p exports
	tar -czf exports/setd1b-results.tar.gz results logs config workflow/Snakefile profiles/gcp-vm

slurm:
	CONDA_CHANNEL_PRIORITY=strict snakemake --snakefile workflow/Snakefile --profile profiles/slurm

clean:
	rm -rf work results logs .snakemake resources/reference data/raw
