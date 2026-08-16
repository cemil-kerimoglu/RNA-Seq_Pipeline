SHELL := /bin/bash

.PHONY: dryrun local slurm clean

dryrun:
	snakemake --snakefile workflow/Snakefile --dry-run --printshellcmds

local:
	snakemake --snakefile workflow/Snakefile --cores 8 --software-deployment-method conda --printshellcmds

slurm:
	snakemake --snakefile workflow/Snakefile --profile profiles/slurm

clean:
	rm -rf work results logs .snakemake resources/reference data/raw
