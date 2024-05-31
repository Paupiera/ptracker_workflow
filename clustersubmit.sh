#!/bin/bash

# The number of CPUs (cores) used by your task. 
#SBATCH --cpus-per-task=1
# The amount of RAM used by your task. 
#SBATCH --mem=15G
# Set a maximum runtime in hours:minutes:seconds. Not required by esrum 
#SBATCH --time=05:00:00
#SBATCH --job-name test_snakemake_workflow
#SBATCH --output=snakemake_output/%x_%j.o
#SBATCH --error=snakemake_output/%x_%j.e

# source activate ~/bxc755/miniconda3/envs/ptracker_pipeline4; 
rm snakemake_output/snakemake.oe
snakemake -p --snakefile snakefile.smk \
--jobs 2 --latency-wait 60 \
--cluster "sbatch  --output={log.o}.cluster --error={log.e}.cluster --time=00:00:10 --job-name {rule}  --cpus-per-task {threads} --mem {resources.mem_gb}G "
&> snakemake_output/snakemake.oe

