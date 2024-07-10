Maximum Compatability
===
## Input data 
- Salmonella data in the folder `clusters`
## Output data 
- Generated maximum compatability tree for Salmonella data in the folder `results`
## Installing packages
```
source /hps/software/opt/linux-rocky8-cascadelake/anaconda3/etc/profile.d/conda.sh
module load python 

# conda create -n ${ENV_NAME}
conda activate ${ENV_NAME}
conda install conda-forge/label/h6f9ffa1_0::gcc setuptools
cd compat-0.8.6
python3 ./setup.py build
python3 ./setup.py install
```
## Running the script on codon-slurm
```
./start.slurm.sh
```

References 
[Link to Paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1520-4#Sec15)