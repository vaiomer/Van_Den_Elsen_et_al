# Van_Den_Elsen_et_al
Python scripts for 16S diversity analyses from article Van Den Elsen et al. (PC0121)


## Virtual environment setup (1st time run only)
Prior to run the analysis scripts, the Python virtual environment needs to be created using the following commands: 
- python2.7 -m "virtualenv" venv
- source venv/bin/activate
- pip install -r requirements.txt 
- deactivate


## 0. Virtual environment activation
source venv/bin/activate


## 1. Alpha diversities
The command './alphaDiversity.py' will create 2 files:
- alpha_diversity_Genus.tsv: the calculated alpha diversity values
- alpha_diversity_kruskal.xlsx: the Kruskal-Wallis statistical test across analysis groups


## 2. Beta diversities
The command './betaDiversity.py' will create 2 directories:
- betaDiversityMatrices: the calculated beta diversity values
- betaDiversityPermanova: the PERMANOVA statistical test across analysis groups


## 3. Ordinations
The command './ordination.py' will create 1 directory:
- betaDiversityOrdination: the MDS ordination images
