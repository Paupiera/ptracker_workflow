rm snakemake_output/snakemake.oe
dir="/maps/projects/rasmussen/scratch/ptracker/ptracker/snakemake_output"

snakemake -n --rerun-incomplete  --snakefile snakefile.smk \
--jobs 2 --latency-wait 60 \
--cluster "sbatch  --output={log.o}.cluster --error={log.e}.cluster --time=00:00:10 --job-name {rule}  --cpus-per-task {threads} --mem {resources.mem_gb}G "

# --cluster "sbatch --time=00:00:10 --output $dir/snakemake_output/%j.o --error $dir/snakemake_output/%j.o --job-name {rule} --cpus-per-task {threads} --mem {resources.mem_gb}G" 
#--error=$dir/%x.out 
