## Librairies
import pandas as pd 
import shutil
import os, sys
from argparse import ArgumentParser



###############################################################################################################
# ## Infos 

""" Analyse des resultats tblastn et mise à jour de la table avec les données de gène tet troncqués"""

# ## Lancement du script

# python3 -i <workdirectory>  -t <refseq table> -b <csv file tblastn> 

###############################################################################################################



# Analyseur d'argument
def config_parameters():
    parser= ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="work directory")
    parser.add_argument("-t", "--table", dest="df", help="refseq table")
    parser.add_argument("-b", "--blast", dest="tblastn", help="tblastn csv file result")
    parser.add_argument("-o", "--out", dest="outputTable", help="final table of refseq")
    options=parser.parse_args()
    if len(sys.argv)< 4:
        sys.exit("Warning : Wrong number of arguments")
    return options.inputdir, options.df, options.tblastn, options.outputTable


def main():
    
    inputdir, df, tblastn, outputTable= config_parameters()
    
    # lecture de la table de sortie tblantn
    read_blast = pd.read_csv(tblastn, header=None)
    

    #Considére uniquement les resultats blast avec un %identidé >= 50 et une couverture d'au moins 30% et recupérer leur numeros d'accessions
    read_blast = read_blast[read_blast[3] >= 50]
    new_file = read_blast[ read_blast[2]>= 30]
    last_file = new_file.drop_duplicates(subset=[1], keep='first')
    genomes = last_file[1].tolist()
    genome_list = []
    #retirer la version pour que les données collent à ceux du tableau generale
    for genome in genomes :
        genome = genome.split(".")
        genome_list.append(genome[0])
    # print(genome_list)


    read_table = pd.read_csv(df, index_col=[0])

    read_table["troncated_tet"] = None

    for index, row in read_table.iterrows() :
        if row["Nucleotide_accession"] in genome_list :
            groupe = row["groupe"]
            specie_name = row["Bacteria_specie"].replace(" ", "_")
            refseq_acc = row["Refseq_accession"]

            file_dir = os.path.join(inputdir + "/g_" + groupe, specie_name)
            file_path = os.path.join(file_dir, refseq_acc + ".csv")
            # Create a new repository 
            new_dir = os.path.join(inputdir + "/gTroncated", specie_name)
            os.makedirs(new_dir, exist_ok=True)
            # print(new_dir)
            
            shutil.move(file_path, new_dir)
            reposit = os.listdir(file_dir)
            if not reposit:
                os.rmdir(file_dir)

            read_table.loc[index, "troncated_tet"] = "yes"
            read_table.loc[index, "groupe"] = "Troncated"
            
        else :
            read_table.loc[index, "troncated_tet"] = "No"

        read_table.to_csv(outputTable)


# read_table.shape
print("Process done !!")
if __name__== "__main__" :
    main()