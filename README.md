# IME_Rho_tet Snakemake Pipeline


**Snakemake pipeline for IME_Rho_tet analysis to run on the high-performance computing (HPC) cluster of the Migale platform.**


[![Migale](https://img.shields.io/badge/Migale-Cluster-red)](https://migale.inrae.fr/cluster)
![Snakemake](https://img.shields.io/badge/Snakemake-8.15.1-yellow)
![Python](https://img.shields.io/badge/Python-3.9%2B-informational?logo=python)
![Conda](https://img.shields.io/badge/Conda-22.11.1-green)
![Bash](https://img.shields.io/badge/Bash-5.2-violet)


---


## 🧭 Overview
This repository contains the **Snakemake Pipeline** developed for the analysis of IME_Rho_tet element carrying _tet(W)_ and _tet(32)_ genes, focusing on their diversity and distribution in bacterial gut genomes from human and animal hosts.

The pipeline automates all steps up to protein clustering ensuring full reproducibility and scalability


---

## 🧩 Features

1. **Download RefSeq Genomic Data** : Migale platform, NCBI FTP taxonomy

2. **In silico search of IME_Rho_tet** : BLASTp, BLASTn, GenoFig

3. **Protein clustering** : ResFinder, CD-HIT, BLASTp, SankeyMATIC

4. **Boundaries and integration site characterization**  
   - TIR detection : Clustal W, FIMO MEME Suite, CD-HIT, BLASTn  
   - Integration site : BLASTn, nr

5. **Sources analysis** : Excel, CIRCOS
   
7. **Statistical tests** : Chi-Square, R


---


## 📁 Directory Structure


IME_Rho_tet_pipeline/

├── configs/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Configuration YAMLs files

├── data/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Example input datasets (e.g., rankedlineage.dmp/)

├── scripts/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# Python-bash for Snk and FIMO ; bash-R for source/hosts isolation and statistical tests

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


## 🚀 Pipeline Launch

Follow these steps to execute the **IME_Rho_tet Snakemake pipeline**:

Exemple : Refseq accession (GCF*) annotations of tree strains
   - Roseburia hominis GCF_000225345.1_ASM22534v1
   - Bifidobacterium longum GCF_033344655.1_ASM3334465v1
   - Streptococcus parasanguinis GCF_024295855.1_ASM2429585v1


### 1️⃣ Download RefSeq Genomic Data


- Retrieve translated amino acid and genbank annotation files of _Bacillota_ and _Actinomycetota_ from RefSeq.  
- Refseq annotation files exemple:  
 
   - path to download refseq annotation of translated CDS file (_translated_cds.faa.gz) and gengank annotation file (_genomic_gbff.gz) :
     
     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/225/345/GCF_000225345.1_ASM22534v1/
     
     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/033/344/655/GCF_033344655.1_ASM3334465v1/
     
     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/699/395/GCF_036699395.1_ASM3669939v1/
     
   - Use "wget" command to download files  : https://doc.ubuntu-fr.org/wget



- Snakemake documentation : snakemake.readthedocs.io/en/stable/
     

     
     
### 2️⃣ In Silico Search of IME_Rho_tet


**Snakefile:** `1_InSilico_IME_Rho_tet.smk`
   - config : `configs/1_config.yaml`
   - scripts :  available in scripts repository


**Inputs:**  
*Example files are available in the `data/inputs` repository.*

- Query file: `data/query_p.fa`   
- RefSeq downloaded files : `data/inputs/GCF_000225345.1_ASM22534v1_translated_cds.fa` and `data/inputs/GCF_000225345.1_ASM22534v1_genomic_gbff`
- Query files for truncated sequences :  `data/query_t.fa` 
- Reference files for PPR (r1), relaxase(r2), and recombinase(r3) :  `data/inputs/r1` `data/inputs/r2` `data/inputs/r3`
- NCBI taxonomie file : `data/rankedlineage_head.dmp`
- NCBI taxonomie reworked for _Bacillota_ and _Actinomycetota_ : `data/inputs/taxonomie_bacillota_actinomycetota`
  

**Outputs:**   
*Example files are available in the `data/outputs` repository.*

- BLASTp result tables after analysis and filtering : `data/outputs/GCF_000225345.1_ASM22534v1` `data/outputs/GCF_024295855.1_ASM2429585v1` `data/outputs/GCF_033344655.1_ASM3334465v1`
- Summary table with metadata  (separator: `,`) : `data/outputs/Summary_table_1`



  
### 3️⃣ Protein Clustering


**Snakefile:** `2_clustering.smk`  
   - config file : `configs/2_config.yaml`


**Inputs:**  
*Example files are available in the `data/inputs` repository.*

- Created database fasta files for PPR, REL and REC : `data/inputs/PPR_DB.fa` `data/inputs/REL_DB.fa` `data/inputs/REC_DB.fa`
  

**Outputs:**   
*Example files are available in the `data/outputs` repository.*  

- Summary table with metadata  (separator: `,`) : `data/outputs/Summary_table_2`


   
### 4️⃣ Boundaries and Integration Site Characterization

#### a) TIR Motif Detection
- Scripts:  
1. `1_run_fimo.sh` – Search TIR motifs using `/scripts/FIMO/TIR_model.meme`  
2. `2_parse_TIRs.py` – Extract TIRs for IME_Rho_tet and RPP(tet) groups  
3. `3_add_TIRs.py` – Add TIRs to summary table  

- Output: Summary table with TIRsequence and coordinate ex : /data/4_summary_table_with_TIR

#### b) Integration Site Analysis
- Scripts:  
1. `1_integrationSite_IME_Rho_tet.py` – Extract 150 bp flanking TIRs  
2. `2_clustering_blastn.sh` + `blastn_analyse.py` – Cluster sequences & search homologs in nr  
3. `3_Refseq_GB_remote.py` – Annotate homologous sequences using RefSeq/GenBank  

### 5️⃣ Sources Analysis
- Metadata about sequence sources is extracted and categorized using bash scripts.


### 6️⃣ Statistical Tests
- Chi-square and other statistical tests are performed using R scripts.


✅ **Tip:** Always run the pipeline in the designated **working directory** and ensure all paths in the configuration files reflect your environment.

---



## 📬 Contact

Created and maintained by Ouléye Sidibé

Questions or Feature requests ? Open an ![issue](https://github.com/OuleyeSidibe/IME_Rho_tet_Snakemake_Pipeline/issues/new) or contact ouleyehelena@gmail.com


## 📚 Citation
If you use  IME_Rho_tet Snakemake Pipeline in your work, please cite: .....



## 📝 License
This project is licensed under the GNU General Public License v3.0 (GPL-3.0).


