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

snakemake -p --snakefile snakefile.smk \
--jobs 8 --max-jobs-per-second 5 --max-status-checks-per-second 5 --latency-wait 60 --keep-going --keep-incomplete \
--cluster "sbatch  --output=snakemake_output/{rule}.%j.o --error=snakemake_output/{rule}.%j.e --time={resources.walltime} --job-name {rule}  --cpus-per-task {threads} --mem {resources.mem_gb}G "
&> snakemake_output/snakemake.oe
