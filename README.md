# Pipeline for running PLAMB WIP
Pipeline for running Plamb: https://github.com/RasmussenLab/vamb/tree/vamb_n2v_asy


## Running the pipeline
```
# clone this repository
git clone https://github.com/Las02/ptracker_workflow -b clean_up_the_code
cd ptracker_workflow

# Clone the Plamb directory and the Plamb helper script directory
mkdir bin
cd bin
git clone https://github.com/RasmussenLab/vamb -b vamb_n2v_asy
git clone https://github.com/Paupiera/ptracker

```



The pipeline can be configurated in: ``` config/config.yaml ```
Here the resources for each rule can be configurated as follows
```
mpSpades:
  walltime: "15-00:00:00"
  threads: 60
  mem_gb: 950
```
if no resourcess are configurated for a rule the defaults will be used which are also defined in: ``` config/config.yaml ```

The input files for the pipeline can be configurated in ``` config/accesions.txt ``` 
As an example this could look like:
```
SAMPLE ID READ1 READ2
Airways 4 reads/errorfree/Airways/reads/4/fw.fq.gz reads/errorfree/Airways/reads/4/rv.fq.gz
Airways 5 reads/errorfree/Airways/reads/5/fw.fq.gz reads/errorfree/Airways/reads/5/rv.fq.gz
```

with an installation of snakemake run the following to dry-run the pipeline
```
	snakemake -np --snakefile snakefile.smk
```
and running the pipeline with 4 threads
```
	snakemake -p -c4 --snakefile snakefile.smk --use-conda
```

```
## Important Files
- snakefile.smk: The snakemake pipeline
- utils.py: utils used by the pipeline
- config: directory with the configuration files
  - accesions.txt: Sample information
  - config.yaml: configuration for the pipeline eg. resourcess
- envs: directory with the conda environment descriptions

## Misc files
Makefile - various small scripts for running the pipeline
clustersubmit.sh - script for submitting the snakefile to SLURM
parse_snakemake_output.py - small script for viewing snakefile logs
