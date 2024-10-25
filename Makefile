SFILE = snakefile.py

all:
	snakemake -np  --snakefile  $(SFILE)  --rerun-incomplete  --rerun-triggers mtime

touch:
	snakemake -p --snakefile  $(SFILE) --touch  --rerun-incomplete -c4

split_fst:
	python split_fasta.py --help

runsnakemake:
	snakemake -p -c4 --snakefile $(SFILE) --rerun-triggers mtime

unlockdir:
	snakemake -p --unlock --snakefile $(SFILE)

clustersubmit:
	mkdir -p snakemake_output
	rm snakemake_output/snakemake_workflow.e
	touch snakemake_output/snakemake_workflow.e
	sbatch clustersubmit.sh
	tail -f snakemake_output/snakemake_workflow.e

parse_snake: 
	python parse_snakemake_output.py

view:
	python parse_snakemake_output.py snakemake_output/snakemake_workflow.e | fzf --height=100%  --ansi --preview "cat {}.e <(echo '------------ std.out ------------') {}.o" | xargs -I {} less -S {}.e {}.o

scancel:
	squeue  --me | grep -v NODELIST | fzf -m --bind ctrl-a:select-all,ctrl-d:deselect-all,ctrl-t:toggle-all | awk '{print $$1}' | xargs scancel




