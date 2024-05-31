#!/bin/bash


# The number of CPUs (cores) used by your task. 
#SBATCH --cpus-per-task=1
# The amount of RAM used by your task. 
#SBATCH --mem=15G
# Set a maximum runtime in hours:minutes:seconds. Not required by esrum 
#SBATCH --time=05:00:00

# conda activate ptracker_pipeline4
snakemake -p --snakefile snakefile.smk \
--use-conda  \
--jobs 2 --latency-wait 60 \
--cluster "sbatch --cpus-per-task {threads} --mem {params.mem_gb}G" &> snakemake_output/snakemake.oe

