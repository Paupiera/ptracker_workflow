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
if no resourcess are configurated for a rule the defaults will be used

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

## Creating the env for plamb
```
# create the conda env from the yaml file: pipeline_conda.yaml
  conda env create --file=pipeline_conda.yaml

# activate the conda env
  conda activate ptracker_pipeline4

# manually install packages for vamb through pip,
# making sure to use the pip-version associated with the conda env
  pip install numpy==1.24.2 torch==1.13.1 pycoverm==0.6.0 networkx==3.1 scikit-learn==1.2.2 pandas==2.0.0 dadaptation==3.0 loguru==0.7.2 fastnode2vec scipy==1.10.1 fastnode2vec pysam

# clone avamb and and pip install it
  git clone https://github.com/RasmussenLab/vamb -b vamb_n2v_asy
  pip install -e path/to/vamb_n2v_asy
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
