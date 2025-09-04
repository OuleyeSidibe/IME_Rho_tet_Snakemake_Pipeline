## Librairies 
from argparse import ArgumentParser
import pandas as pd 
import os, sys, re 



###############################################################################################################################################################################################
# ## Notes

""" Remove pseudogenes data in RefseqFNALE2 """


## Load script :

# python3 /remove_pseudoData.py -i <inputdir> -o <outputdir> -p <pseudogenes file data>
##inputdir ="/home/osidibe/work/PPR_MGEproject/snakemakes/blastp_Groups_2/RefseqFINALE.csv"
##outputdir="/home/osidibe/work/PPR_MGEproject/snakemakes/clustering_by_prot_3_2025Article/RefseqFINALE2.csv"
##pseudo = "/home/osidibe/work/PPR_MGEproject/snakemakes/clustering_by_prot_3_2025Article/report_3.stdout"

################################################################################################################################################################################################


# Analyseurs d'arguments
def config_parameters():
    parser=ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="inputdir groups for protéins db")
    parser.add_argument("-o", "--output", dest="outputdir", help="outputdir for proteins groups of cluster")
    parser.add_argument("-p", "--pseudo", dest="pseudo", help="pseudogenes file data")
    args=parser.parse_args()
    if len(sys.argv) < 3 :
        sys.exit("Warning : wrong number of argument")
    return args.inputdir, args.outputdir, args.pseudo



inputdir, outputdir, pseudo = config_parameters()


table = pd.read_csv(inputdir)

# recuperer les accession des pseudo_génes
list_acc = []

with open (pseudo, "r") as pseudo_file :
    lines = pseudo_file.read().splitlines()

for line in lines :
    if ">lcl" in line :
        acc = re.split(r'\||\.1_prot', line)

        acc_gen = acc[1]
        list_acc.append(acc_gen)

# print(len(list_acc))
    

# Enlver les pseudo_genes de la table RefSeq
for pseudo in list_acc :
    
    for index, row in table.iterrows() :
        if row["Nucleotide_accession"] == pseudo :
            table = table.drop(index)
        
print(len(table))
table.to_csv(outputdir, index=False)