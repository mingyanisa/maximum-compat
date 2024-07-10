Maximum Compatability
===
## Input data 
- Salmonella data in the folder `clusters`
## Installing packages
```
source /hps/software/opt/linux-rocky8-cascadelake/anaconda3/etc/profile.d/conda.sh

# conda create -n ${ENV_NAME}
conda activate ${ENV_NAME}
conda install conda-forge/label/h6f9ffa1_0::gcc
cd compat-0.8.6
pip3 install -e .
```
## Running the script on codon-slurm
```
./start.slurm.sh
```

References 
[Link to Paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1520-4#Sec15)