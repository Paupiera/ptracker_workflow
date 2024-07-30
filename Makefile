all:
	snakemake -np --snakefile snakefile.py
	
runsnakemake:
	snakemake -p -c4 --snakefile snakefile.py

unlockdir:
	snakemake -p --unlock --snakefile snakefile.py

clustersubmit:
	mkdir -p snakemake_output
	rm snakemake_output/snakemake_workflow.e
	touch snakemake_output/snakemake_workflow.e
	sbatch clustersubmit.sh
	tail -f snakemake_output/snakemake_workflow.e

