## Librairies 
from argparse import ArgumentParser
import pandas as pd 
import os, sys, re 



###############################################################################################
# ## Notes

""" Remove pseudogenes data in summary table """

## Load script :

# python3 /remove_pseudoData.py -i <inputdir> -o <outputdir> -p <pseudogenes file data>

################################################################################################


# Arguments
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

# pseudogene accessions
list_acc = []

with open (pseudo, "r") as pseudo_file :
    lines = pseudo_file.read().splitlines()

for line in lines :
    if ">lcl" in line :
        acc = re.split(r'\||\.1_prot', line)

        acc_gen = acc[1]
        list_acc.append(acc_gen)

# print(len(list_acc))


# remove pseudogenes data on summary table
for pseudo in list_acc :
    
    for index, row in table.iterrows() :
        if row["Nucleotide_accession"] == pseudo :
            table = table.drop(index)
        
print(len(table))
table.to_csv(outputdir, index=False)
