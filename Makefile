all:
	snakemake -np --snakefile snakefile.py
	
runsnakemake:
	snakemake -p -c4 --snakefile snakefile.py
