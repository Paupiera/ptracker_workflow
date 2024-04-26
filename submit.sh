snakemake -np --jobs 40 -k --snakefile snakefile.py -p --use-conda --latency-wait 60 --rerun-incomplete \
  --cluster "qsub -l walltime={resources.walltime},nodes=1:thinnode:ppn={threads},mem={resources.mem_gb}gb"\
" -W group_list=ku_00197 -A ku_00197 " \
# TODO add qsub-status.py script thing
