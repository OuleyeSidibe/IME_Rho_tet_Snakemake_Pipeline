## Run in qlogin
# qlogin -q long.q -pe thread 30

## Execute snakemake environnement
# conda activate snakemake-8.15.1

## Snk script execution in snakefile repository
# snakemake -s <<screening.smk>> -c30


## use conda environment by load profile.d/conda.sh
shell.executable("/bin/bash")
shell.prefix("source /usr/local/genome/Anaconda3/etc/profile.d/conda.sh;")

##librairies
import glob
import os


# ## Configure yaml file

configfile: "config.yaml"


""" Note : To launch snk file, positionning in workdir repository """


## Functions

# Define output files for rule2
def output_files(inputdir, threshold, prot_name):
    output=[]

    for prot in prot_name :
        file = os.path.join(f"{inputdir}/{prot}_db.fa")
        file_name = os.path.basename(file)
        out_name = f"{file_name.rsplit('.',1)[0]}_{threshold}"
        output.append(f"{inputdir}/{out_name}")
    return output




######################################################################---RULES---######################################################################
rule all:
    input:
        output_files(config['outputdir'], config['threshold'], config['prot_name']),
        f"{config['workdir']}/clstr_report_3.txt",
        f"{config['outputdir']}/name_unknowProt.txt",
        f"{config['workdir']}/RefseqFINAL2.csv",
        f"{config['workdir']}/RefseqFINAL3.csv"


## Rule1 : clusterised proteins database

rule clustering:
    input:
        f"{config['outputdir']}"

    output:
        output_files(config['outputdir'], config['threshold'], config['prot_name'])

    log:
        out = f"{config['workdir']}/clust.stdout",
        err = f"{config['workdir']}/clust.stderr"

    params:
        t = config['threshold']

    priority: 100

    shell: """
           conda activate cd-hit-4.8.1
           bash {config[workdir]}/clustering.sh {input} {params.t} 1> {log.out} 2> {log.err}
           conda deactivate
           """


## Rule2 : Processin and modify clusters

rule analyse_cluster:
    input:
        marq = f"{config['workdir']}/clust.stderr",
        i = f"{config['outputdir']}"

    output:
        o = f"{config['workdir']}/clstr_report_3.txt"

    params:
        p = config['prot_name']
    

    shell: """
           conda activate base
           touch {output.o}
           python3 {config[workdir]}/analyse_cluster.py -i {input.i} -t {config[refseq_t]} -p {params.p} -o {output.o} -th {config[threshold]}
           conda deactivate
           """


#rule3 : create query file of not referenced tet

rule create_query:
    input:
        marq = f"{config['workdir']}/clstr_report_3.txt",
        i1 = f"{config['outputdir']}"

    output:
        o1 = f"{config['outputdir']}/unknow_TET.fa"

    shell: """
          conda activate base
          touch {output.o1}
          python3 {config[workdir]}/create_query.py -i {input.i1}/TET_db_{config[threshold]}newRef.fa -o {output.o1}
          conda deactivate
          """


 ## Rule4 : blastp of unknows TET with Resfinder local DB

rule blastp:
    input:
        marq = f"{config['workdir']}/clstr_report_3.txt",
        db = f"{config['workdir']}",
        qu1 = f"{config['outputdir']}/unknow_TET.fa"

    output:
        f"{config['outputdir']}/Blastp_unknowProt.csv"

    shell: """
           conda activate blast-2.13.0
           touch {output}
           blastp -subject {input.db}/ResFinder_PPR_Moz.fa -query {input.qu1} -outfmt "6 delim=, qacc sacc qcovs pident" -out {output}
           conda deactivate
           """


# ## Rule5 : blastp analyse and define new name of unknows TET

rule balstp_analysis:
    input:
        f"{config['outputdir']}/Blastp_unknowProt.csv"

    output:
        o1 = f"{config['outputdir']}/name_unknowProt.txt"

    params:
        T1 = f"{config['outputdir']}/TET_db_{config['threshold']}newRef.fa",

    shell: """
           conda activate base
           touch {output.o1}
           python3 {config[workdir]}/blastp_analysisUnk.py -i {input} -o {output.o1} -c {params.T1} -w {config[outputdir]}
           conda deactivate
           """


# ## Rule6 : remove pseudogenes data in RefseqFINAL2 table

rule remove_pseudogenesData:
    input:
        i = f"{config['refseq_t']}",
        p = f"{config['workdir']}/report_3.stdout"

    output:
        o = f"{config['workdir']}/RefseqFINAL2.csv"

    shell: """
           conda activate base
           touch {output.o}
           python3 {config[workdir]}/remove_pseudoData.py -i {input.i} -o {output.o} -p {input.p}
           conda deactivate
           """


# ## Rule7 : add clusters data in Refseq table of clustering prot

rule add_clstrData:
    input:
        maq = f"{config['outputdir']}/name_unknowProt.txt",
        o = f"{config['outputdir']}"

    output:f"{config['workdir']}/RefseqFINAL3.csv"
    
    params:f"{config['workdir']}/RefseqFINAL2.csv"

    shell: """
           conda activate base
           python3 {config[workdir]}/add_clst_RefseqFinal.py -i {input.o} -o {output} -t {params} -g {config[groups]}
           conda deactivate
           """


""" Go to FIMO scripts """ 






