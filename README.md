# The pipeline for running PLAMB

## Important Files
- snakefile.smk: The snakemake pipeline
- utils.py: utils used by the pipeline
- config: directory with the configuration files
  - accesions.txt: Sample information
  - config.yaml: configuration for the pipeline eg. resourcess
- envs: directory with the conda environment descriptions

## Creating the env for plamb
```
# create the conda env from the yaml file: pipeline_conda.yaml
  conda env create --file=pipeline_conda.yaml

# activate the conda env
  conda activate ptracker_pipeline4

# manually install packages for vamb through pip,
# making sure to use the pip-version associated with the conda env
  pip install numpy==1.24.2 torch==1.13.1 pycoverm==0.6.0 networkx==3.1 scikit-learn==1.2.2 pandas==2.0.0 dadaptation==3.0 loguru==0.7.2 fastnode2vec scipy==1.10.1

# clone avamb and and pip install it
  git clone https://github.com/RasmussenLab/vamb -b vamb_n2v_asy
  pip install -e path/to/vamb_n2v_asy
```
## Misc files
Makefile - various small scripts for running the pipeline
clustersubmit.sh - script for submitting the snakefile to SLURM
parse_snakemake_output.py - small script for viewing snakefile logs
