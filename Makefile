all:
	snakemake -np --snakefile snakefile.py --rerun-triggers mtime
	
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
	python parse_snakemake_output.py snakemake_output/snakemake_workflow.e | fzf --height=19% --ansi --preview "cat {}.e" | xargs -I {} less -S {}.e {}.o

scancel:
	squeue  --me | grep -v NODELIST | fzf -m --bind ctrl-a:select-all,ctrl-d:deselect-all,ctrl-t:toggle-all | awk '{print $$1}' | xargs scancel



