# IME_Rho_tet Snakemake Pipeline


**Snakemake pipeline for IME_Rho_tet analysis to run on the high-performance computing (HPC) cluster of the Migale platform.**


[![Migale](https://img.shields.io/badge/Migale-Cluster-red)](https://migale.inrae.fr/cluster)
![Snakemake](https://img.shields.io/badge/Snakemake-8.15.1-yellow)
![Python](https://img.shields.io/badge/Python-3.9%2B-informational?logo=python)
![Conda](https://img.shields.io/badge/Conda-22.11.1-green)
![Bash](https://img.shields.io/badge/Bash-5.2-violet)


---


## Overview
This repository contains the **Snakemake Pipeline** developed for the analysis of IME_Rho_tet element carrying _tet(W)_ and _tet(32)_ genes, focusing on their diversity and distribution in bacterial gut genomes from human and animal hosts.

The pipeline automates all steps up to protein clustering ensuring full reproducibility and scalability


---

## Features

1. **Import RefSeq Genomic Data** : Migale platform, NCBI FTP taxonomy

2. **In silico search of IME_Rho_tet** : BLASTp, BLASTn, GenoFig

3. **Protein clustering** : ResFinder, CD-HIT, BLASTp, SankeyMATIC

4. **Boundaries and integration site characterization**  
   - TIR detection : Clustal W, FIMO MEME Suite, CD-HIT, BLASTn  
   - Integration site : BLASTn, nr

5. **Sources analysis** : Excel, CIRCOS
   
7. **Statistical tests** : Chi-Square, R


---


## Directory Structure


IME_Rho_tet_pipeline/

├── configs/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Configuration YAMLs files

├── data/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Example input datasets (e.g., rankedlineage.dmp/)

├── scripts/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Python-bash use in Snk and FIMO ; bash-R use in source/hosts isolation and statistical tests

├── Snakefiles/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Snakefiles for pipeline steps (import_data, IN_silico IME_Rho_tet, Proteins Clustering)

└── README.md

---


## 🧾 Instructions

Before executing the workflow, please follow these guidelines:

1. ⚙️ Use **Conda** for environment management to ensure proper installation of Python dependencies and the **Snakemake** package.
   
2. 📂 Execute the **Snakefiles** sequentially, respecting the pipeline order (**step 1 → step 2 → step 3**).
   
3. 📌 Always navigate to the designated **working directory** prior to workflow execution.
   
4. 🛠️ Adjust the **file and script directory paths** in the configuration file to reflect your working environment.  

---


## Contact

Created and maintained by Ouléye Sidibé

Questions or Feature requests ? Open an ![issue](https://github.com/OuleyeSidibe/IME_Rho_tet_Snakemake_Pipeline/issues/new) or contact ouleyehelena@gmail.com


## Citation
If you use  IME_Rho_tet Snakemake Pipeline in your work, please cite: .....



## License
This project is licensed under the GNU General Public License v3.0 (GPL-3.0).


# Project Title

[![Overview](https://img.shields.io/badge/Overview-Project%20summary-blue?style=flat-square)](#overview)
[![Features](https://img.shields.io/badge/Features-Highlights-%232ecc71?style=flat-square)](#features)
[![Directory](https://img.shields.io/badge/Directory-structure-lightgrey?style=flat-square)](#directory-structure)
[![Contact](https://img.shields.io/badge/Contact-Get%20in%20touch-orange?style=flat-square)](mailto:ton.email@exemple.com)
[![Citation](https://img.shields.io/badge/Citation-How%20to%20cite-blueviolet?style=flat-square)](<URL-to-citation-or-DOI>)
[![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](LICENSE)


