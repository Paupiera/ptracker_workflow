all:
	snakemake -np --snakefile snakefile.py
	
runsnakemake:
	snakemake -p -c4 --snakefile snakefile.py

unlockdir:
	snakemake -p --unlock --snakefile snakefile.py

clustersubmit:
	rm -rf results/Airways/vamb_from_strobealign_default_params_1
	mkdir -p snakemake_output
	rm snakemake_output/snakemake_workflow.e
	touch snakemake_output/snakemake_workflow.e
	sbatch clustersubmit.sh
	tail -f snakemake_output/snakemake_workflow.e

parse_snake: 
	python parse_snakemake_output.py

view:
	python parse_snakemake_output.py tst.e | fzf --height=20% --ansi --preview "cat {}.e" | xargs -I {} less -S {}.e {}.o


