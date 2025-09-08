## Run in Qlogin
# qlogin -q long.q -pe thread 30

## Execute snakemake environnement
# conda activate snakemake-7.5.0

## Snk script execution in snakefile repository
# snkemake -S <<screening.smk>> -C 30 -k


## use conda environment by load profile.d/conda.sh
shell.executable("/bin/bash")
shell.prefix("source /usr/local/genome/Anaconda3/etc/profile.d/conda.sh;")



# ## Configure yaml file
configfile: "config.yaml"


""" Note : To launch snk file, positionning in workdir repository """



## Rule1 : export amino acide translated files of Bacillota and actinomycetota phyla from migale server """

rule export_data:
    input:
        i = config['workdir'] ,
        db = config['inputdir']
    output:
        o = directory(config['outputdir']),
        out = expand("{stdout}/stdout_report_1_2.txt", stdout=config['workdir'])

    shell:  """
            touch {output.out} | 
            python3 {config[workdir]}/Migale_Data.py -d {input.i} -o {output.o} -db {input.db} -out {output.out}
            rm {config[outputdir]}/.snakemake_timestamp
            """

""" go to 2_InSilico_IME_Rho_tet snakefile """
