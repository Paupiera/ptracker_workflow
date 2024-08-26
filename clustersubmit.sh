#!/bin/bash

# The number of CPUs (cores) used by your task.
#SBATCH --cpus-per-task=1
# The amount of RAM used by your task.
#SBATCH --mem=15G
# Set a maximum runtime in hours:minutes:seconds. Not required by esrum
##SBATCH --time=50:00:00:00
#SBATCH --job-name snakemake_workflow
#SBATCH --output=snakemake_output/%x.o
#SBATCH --error=snakemake_output/%x.e

## TOADD
# --keep-going

# REMOVE
# --conda-base-path ~/bxc755/miniconda3 --use-conda
# module load snakemake/7.30.1

# Just added use-conda..  and base path remove if not works and run scapp alone

snakemake --rerun-triggers mtime --nolock --conda-base-path ~/bxc755/miniconda3 --use-conda --rerun-incomplete -p --snakefile snakefile.py --keep-incomplete --use-conda --rerun-triggers mtime \
  --jobs 2 --max-jobs-per-second 5 --max-status-checks-per-second 5 --latency-wait 60 --keep-incomplete \
  --cluster "sbatch  --output=snakemake_output/{rule}.%j.o --error=snakemake_output/{rule}.%j.e --time={resources.walltime} --job-name {rule}  --cpus-per-task {threads} --mem {resources.mem_gb}G "
&>snakemake_output/snakemake.oe
