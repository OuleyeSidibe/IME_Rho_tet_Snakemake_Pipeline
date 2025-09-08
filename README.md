# IME_Rho_tet Snakemake Pipeline

**Bioinformatic pipeline of IME_Rho_tet diversity and distribution in human and animal bacterial gut**

![Snakemake](https://img.shields.io/badge/Snakemake-8.15.1-yellow)
![Python](https://img.shields.io/badge/Python-3.9%2B-informational?logo=python)
![Conda](https://img.shields.io/badge/Conda-22.11.1-green)
![Bash](https://img.shields.io/badge/Bash-5.2-violet)

---

## Overview
This repository contains the **Snakemake Pipeline** developed for the analysis of IME_Rho_tet element carrying _tet(W)_ and _tet(32)_ genes, focusing on their diversity and distribution in bacterial gut genomes from human and animal hosts.

The pipeline automates all steps up to protein clustering and includes to ensuring full reproducibility and scalability

---

## Features

1. **Import RefSeq Genomic Data** : Migale platform, NCBI FTP taxonomy

2. **In silico search of IME_Rho_tet** : BLASTp, BLASTn, GenoFiga

3. **Protein clustering** : ResFinder, CD-HIT, BLASTp, SankeyMATIC

4. **Boundaries and integration site characterization**  
   - TIR detection : Clustal W, FIMO MEME Suite, CD-HIT, BLASTn  
   - Integration site : BLASTn, nr

5. **Host–microbiota association profiling** : CIRCOS

---

## Directory Structure


IME_Rho_tet_pipeline/

├── config/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Configuration YAMLs files


├── data/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Example input datasets (e.g., rankedlineage.dmp/)


├── Picture/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Pipeline picture


├── scripts/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Python and bash scripts use in Snakemake


├── Snakefiles/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Snakefiles for pipeline steps (import_data, IN_silico IME_Rho_tet, Proteins Clustering)


└── README.md



## Contact
Created and maintained by Ouléye Sidibé
Questions or Feature requests ? Open an issue or contact ouleyehelena@gmail.com


## Citation

## License
This project is licensed under the GNU General Pu

## 
