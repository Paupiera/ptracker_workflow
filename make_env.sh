
# create the conda env from the yaml file
conda env create --file=pipeline_conda.yaml
# activate the conda env
conda activate ptracker_pipeline4

# manually install packages for vamb through pip as these are wierd in the pau version
pip install numpy==1.24.2 torch==1.13.1 pycoverm==0.6.0 networkx==3.1 scikit-learn==1.2.2 pandas==2.0.0 dadaptation==3.0 loguru==0.7.2 

# cd to the vamb dir and install it or just give it the path eg. 
pip install -e bin/vamb 

