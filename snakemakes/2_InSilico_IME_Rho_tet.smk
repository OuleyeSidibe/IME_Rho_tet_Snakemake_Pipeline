## Run in Qlogin
# qlogin -q long.q -pe thread 30

## Execute snakemake environnement
# conda activate snakemake-8.15.1

## Snk script execution in snakefile repository
# snkemake -s <<createGroups.smk>> -C 30 -k


## use conda environment by load profile.d/conda.sh
shell.executable("/bin/bash")
shell.prefix("source /usr/local/genome/Anaconda3/etc/profile.d/conda.sh;")



# ## Configure yaml file
configfile: "config.yaml"



# ## Dictionnary of specie name with their annoted refseq file (GCF)

SPECIES=os.listdir(config['outputdir']['migale_data'])
dict_specie={}
for specie in SPECIES:
    # print(specie)
    specie_path = os.path.join(config['outputdir']['migale_data'], specie)
    list_refseq_annotation = os.listdir(specie_path)
    dict_specie[specie]=list_refseq_annotation




""" Note : To launch snk file, positionning in workdir repository """


######################################################################---RULES---######################################################################


rule all:
    input:  
        [f"{config['outputdir']['blastp_Groups_2']}/blastp/{specie}/{GCF_file}" for specie, annotations in dict_specie.items() for GCF_file in annotations],
        expand("{blastp_Groups_2}/groupes", blastp_Groups_2=config['outputdir']['blastp_Groups_2']),
        expand("{directory}/gb_files", directory=config['outputdir']['blastp_Groups_2']),
        f"{config['outputdir']['blastp_Groups_2']}/tBlastn/blastTronc.csv",
        f"{config['workdir']}/Genre_ab_Absolue_Relative.csv",
        f"{config['outputdir2_bis']}"

    


"""  Create a query file of proteins of interest (TetW, Relaxase, recombinase) in workdir """

# Rule1 : Proteins blast of query file on import migale data

rule blastp:
    input:
        qu = f"{config['workdir']}/query_p.fa",
        db = lambda wildcards: expand("{outputdir}/{specie}/{GCF_file}", outputdir=config['outputdir']['migale_data'], specie=wildcards.specie, GCF_file=wildcards.GCF_file)
        
    output: 
        b = f"{config['outputdir']['blastp_Groups_2']}/blastp/{{specie}}/{{GCF_file}}"

    threads: 1

        
    shell:  """
            conda activate blast-2.13.0 
            mkdir -p $(dirname {output.b}) 
            blastp -query {input.qu} -subject {input.db} -outfmt "6 delim=, qacc sacc qcovs pident evalue" -out {output.b}
            conda deactivate
            """



# Rule2 : blastp results analysis to identify RPP with/without relaxase and recombinase

rule blastp_analysis:
    input:
        out = expand("{stdout}/stdout_report_1_2.txt", stdout=config['outputdir']['import_Data_1'])
    output:
        o = directory(expand("{directory}/groupes", directory=config['outputdir']['blastp_Groups_2']))

    params:
        inp = f"{config['outputdir']['blastp_Groups_2']}/blastp",
        pi = config['min_identity'],
        pc = config['min_coverage'],
        tet = config['TET'],
        interv = config['CDS_interval']

    shell: """
           python3 {config[workdir]}/blastp_analysis.py -i {params.inp} -o {output.o} -pi {params.pi} -pc {params.pc} -t {params.tet} -int {params.interv} -a {input.out}
           """


# Rule3 : Summary table with metaData

rule refseq_table:
    input:
        i = expand("{directory}/groupes", directory=config['outputdir']['blastp_Groups_2']),
        gb = config['inputdir']
         

    output:
        o = directory(expand("{directory}/gb_files", directory=config['outputdir']['blastp_Groups_2'])),
        db = f"{config['outputdir']['blastp_Groups_2']}/db_troncSearch.txt"

    params:
        tet = config['TET'],
        mge = config ['MGE']

    shell: """
           conda activate base
           touch {output.db} 
           python3 {config[workdir]}/refseq.py -i {input.i} -gb {input.gb} -t {params.tet} -m {params.mge} -o {output.o} -db {output.db}
           conda deactivate
           """


""" Create a query file of TetW protein sequence in workdir """

## Rule4 : SEarch troncated tet(W) genes in IME_Rho_tet with tblastn tool and isolate them 

rule tblastn:
    input: 
        db = f"{config['outputdir']['blastp_Groups_2']}/db_troncSearch.txt",
        qu = f"{config['workdir']}/query_t.fa"

    output:f"{config['outputdir']['blastp_Groups_2']}/tBlastn/blastTronc.csv"

    shell : """
            conda activate blast-2.13.0
            mkdir -p $(dirname {output})
            tblastn -subject {input.db} -query {input.qu} -outfmt "6 delim=, qacc sacc qcovs pident evalue sstart send" -out {output} 
            conda deactivate
            """


# Rule5 : tblastn result analysis and summary table updating 

rule tblastn_analysis:
    input:
        i = expand("{directory}/groupes", directory=config['outputdir']['blastp_Groups_2']),
        b = f"{config['outputdir']['blastp_Groups_2']}/tBlastn/blastTronc.csv"
    output:
        o = f"{config['workdir']}/RefseqFINALE.csv"

    params:
        df = expand("{directory}/groupes/Refseq_2.csv", directory=config['outputdir']['blastp_Groups_2']) 

    shell: """
           python3 {config[workdir]}/tblastn_analysis.py -i {input.i} -t {params.df} -b {input.b} -o {output.o}
           """


# Rule6 :  File of relative abundances by genus per group - Absolute abundance by genus - Number of genera per group

rule abondance_genus:
    input:
        i = f"{config['workdir']}/RefseqFINALE.csv"

    output:
        f"{config['workdir']}/Genre_ab_Absolue_Relative.csv"

    params:
        r = f"{config['report_1']}",
        tax = f"{config['taxonomie_1']}"

    shell: """
           python3 {config[workdir]}/metaData_abund.py -o {config[workdir]} -t {input.i} -r {params.r} -tax {params.tax}
           """



## Rule7 : Create proteins database in each groups

 rule DB_cluster:
     input: 
         i = f"{config['outputdir']['blastp_Groups_2']}",
         f = f"{config['outputdir']['migale_data']}",
         r1 = f"{config['ref_file']['r1']}",
         r2 = f"{config['ref_file']['r2']}",
         r3 = f"{config['ref_file']['r3']}"

     output:
         directory(f"{config['outputdir2']}")
    
     params:
         p = config['proteins'],
         g = config['groups']
    
     shell: """
            conda activate base
            python3 {config[workdir]}/DB_cluster.py -i {input.i}/groupes -o {output} -f {input.f} -p {params.p} -g {params.g} -r1 {input.r1} -r2 {input.r2} -r3 {input.r3}
            conda deactivate
            """


""" go to 3_clustering snakefile """
