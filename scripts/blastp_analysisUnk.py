## Librairies
import pandas as pd
from argparse import ArgumentParser
import sys, os
import re 
import random
from Bio import SeqIO




########################################################################################################################
# ## Notes

""" Filtrer les resultats blastp afin de définir un nouveau non pour les nouveau cariant Tet non connus de Resfinder """


## Load script :

# python3 /blastp_analysis.py -i <inputdir> -o <outputdir>  -c <tet new ref cluster>



#########################################################################################################################



# Analyseurs d'arguments
def config_parameters():
    parser=ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="input file for blastp result")
    parser.add_argument("-o", "--output", dest="outputdir", help="output file for report of bastp analyse")
    parser.add_argument("-c", "--cluster", dest="cluster_TET", help="fasta file for g_TET cluster")
    parser.add_argument("-w", "--workdir", dest="workdir", help="work directory of clustering")
    args=parser.parse_args()
    if len(sys.argv) < 4 :
        sys.exit("Warning : wrong number of argument")
    return args.inputdir, args.outputdir, args.cluster_TET, args.workdir


# generateur de nom des nouveaux tet
def generate_unique_TETnumber(start, end, generated_numbers):
    while True:
        num = random.randint(start, end-1)
        if num not in generated_numbers:
            generated_numbers.add(num)
            return num

    
inputdir, outputdir, cluster_TET, workdir = config_parameters()

# Read and named columns of blastp result file
read_blast = pd.read_csv(inputdir, header=None, index_col=[0])
read_blast.columns = ["Subject", "Coverage", "Identity"]

# Put unknow protein ID in list
query = read_blast.index.unique().tolist()
nb_query = len(query)




## Open outputfile for append 
with open(outputdir, "a") as name_file :
    name_file.write(f"{nb_query} Tetracycline(s) protein(s) not associated with Proteins referenced in ResFinder database !!\n\n\n")    
    
    # For each unknow prot blastp result, check if it's new variant or new tetracycline protein
    generated_numbers=set()
    for prot in query :
        
        match_prot=0
        new_prot=""
        str_query = re.split(r'\__|\_Tet', prot)
        id_query= str_query[1]
        print(id_query)
        

        data = read_blast.loc[prot,]
        for index, row in data.iterrows():

            if row["Coverage"] >= 95 and row["Identity"] >= 80 :
                name = row["Subject"]
                identity = row["Identity"]
                name_file.write(f"{index} is probably new variant of referenced prot ind Resfinder !\n")  
                name_file.write(f" ---> {index} = {name} identical with {identity}%\n")
                new_prot=f"{name}_newVar"
                name_file.write(f" ---> Tetracycline new variant = {new_prot}\n\n")
                name = index 
                match_prot+=1
                break
            
                
            
        if not match_prot :
            name_file.write(f"{index} is new variant of not referenced prot ind Resfinder, Soo create a name ;)\n")
            new_prot_number = generate_unique_TETnumber(70, 150, generated_numbers)
            new_prot = f"Tet{new_prot_number}_new"
            name_file.write(f" ---> Tetracycline new protein = {new_prot}\n\n")
            

    
        ## open cluster final file to define new proteine name
        """ ici recupére le numero du cluster"""
        with open(cluster_TET, "r") as tet_newref_fa, open(f"{workdir}/tmp.txt", "a") as tmp_fa:
            for record in SeqIO.parse(tet_newref_fa, 'fasta') :
                if id_query in record.id :
                    numb_clstr = record.id.rsplit("_",1)
                    # print("yessssssssssssssssssssssssssssssssss")
                    # print(id_query)
                    # print(record.id)
                    id_record = record.id
                    # if match_prot :
                    #     record.id = f"{new_prot}"
                    # else:
                    #     record.id = f"{new_prot}_{numb_clstr[1]}."
                    record.id = f"{new_prot}_{numb_clstr[1]}"
                    record.description=""
                    SeqIO.write(record, tmp_fa, 'fasta')
    
    with open(cluster_TET, "r") as tet_newref_fa, open(f"{workdir}/tmp.txt", "a") as tmp_fa:
        for record in SeqIO.parse(tet_newref_fa, 'fasta') :   
            if "WP" not in record.id:
                SeqIO.write(record, tmp_fa, 'fasta')
            # repalce original file by tmp file
        os.replace(f"{workdir}/tmp.txt", cluster_TET)
                
            
print("Process Done !!")
