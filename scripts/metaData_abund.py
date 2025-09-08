## Libraries
from argparse import ArgumentParser
import pandas as pd
import os, re
import glob
import sys



################################################################################################################################################################################
# ## Infos 

""" abundance metadata tables for visualisation """

# ## Lancement script :

# python3 ~/metaData_abund.py -o "outputdir" -t "refseq_table finale obtain in step 2 of balstgroups" -r "report migale file obtain in step1" -tax "NCBI taxoResume_2 in step 1"


#################################################################################################################################################################################




# Analyseur d'arguments
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-o", "--out", dest="outputdir", help="output directory")
    parser.add_argument("-t", "--table", dest="table_ref", help="RefSeq finale table1")
    parser.add_argument("-r", "--report", dest="report_Migale", help="Report migale table")
    parser.add_argument("-tax", "--taxonomie", dest="taxonomie", help="taxonomie table")


    
    options=parser.parse_args()
    if len(sys.argv) < 4:
        sys.exit("Warning : wrong number of arguments")
    return options.outputdir, options.table_ref, options.report_Migale, options.taxonomie




outputdir, table_ref, report_Migale, taxonomie = config_parameters()




################################################## Fichier des abondances relatives par genre par groupe ########################################



# 1- Créer la table de l'abondance absolue et relative des genres bacteriens analysés dans les différents groupes (IME_Rho_tet /tet/autres)

table= pd.read_csv(table_ref, index_col=[-2])
groupe_name = table.index.unique()


dict_genre={}
for groupe in groupe_name:
    list_genre=[]

    table_groupe = table[table.index == groupe]

    for genre in table_groupe["Bacteria_specie"]:
            
        genre_name = genre.split(" ")

            
        list_genre.append(genre_name[0])

    table_groupe["Bacteria_specie"] = list_genre
    table_genre = table_groupe["Bacteria_specie"].value_counts()
    # print(groupe)
    for index, line in table_genre.items():
        if dict_genre.get(groupe) is not None :
                dict_genre[groupe] += [[index, line]]
        else:
                dict_genre[groupe] = [[index, line]]

    # print(dict_genre)

# Créer la table à partir du dictionnaire
df= pd.DataFrame(columns=["groupe", "genre", "Ab_relative"])
ab_rel={}
for key, value in dict_genre.items():
    for val in value :
        df2 = pd.DataFrame(data=[[key, val[0], val[1]]], columns=["groupe", "genre", "Ab_relative"])

        df = pd.concat([df, df2], ignore_index=True)
        
        # créer un dictionnaire des abondances relatives par genre
        if ab_rel.get(val[0]) is not None :
                ab_rel[val[0]] += [val[1]]
        else:
                ab_rel[val[0]] = [val[1]]


list_no_dupl_genre = df["genre"].drop_duplicates().tolist()

df.to_csv(f"{outputdir}/Genre_abondGroupe.csv", index=False)



# 2_ Créer la table du nombre de genre par groupe

df_2 = df.groupby("groupe").size().reset_index(name="nb_genre")

df_2.to_csv(f"{outputdir}/nbr_Genre_byGroupe.csv", index=False)






################################################## Fichier des abondances absolues et relatives par genre #####################################


# 3_ créer la table de l'abondance absolue et relative de chaque genre bacterien dans les données RefSeq

## Lire la table report_Migale et la table taxonomie
read_df = pd.read_csv(report_Migale, index_col=[0], sep ='\t')

taxo_df = pd.read_csv(taxonomie, index_col=[4])

## Retenir que les noms de genres dans les index
list_genre=[]

for index in read_df.index:
        index2 = index.split("_")
        index2 = index2[0]
        list_genre.append(index2)
set_genre = set(list_genre)

## liste des genres
list_genreF = list(set(list_genre))
# print(list_genreF)
        
dict_analysedGenus={}

for genre in list_genreF:
    if genre in list_no_dupl_genre :
        nb_genome=0
        genus = "" 
        # get absolute abundance
        for index, row in read_df.iterrows():
            if genre in index :
                nb_genome += row["Number of v1 annotation(s)"]
                genus = True
                
        # get relative abundance 
        for key, value in ab_rel.items():
            if genre == key :
                nb_relative_genome_Perc= (sum(value) * 100) / nb_genome
                    
                        
                
                
        #get family name and phylum
        if genus :
            for index, row in taxo_df.iterrows() :
                if index == genre:
                    family_name = row['family']
                    phylum = row["phylum"]
                    break

        dict_analysedGenus[genre]=nb_genome, nb_relative_genome_Perc, family_name, phylum
                        

## Create table for visualisation plot
table = pd.DataFrame(columns=["genre", "Ab_absolue", "Ab_relative", "family", "phylum"])

for key, value in dict_analysedGenus.items():
    df = pd.DataFrame(columns=["genre", "Ab_absolue", "Ab_relative", "family" ,"phylum"], data=[[key, value[0], value[1], value[2], value[3]]])
    table = pd.concat([table, df], ignore_index=True)


# table = table.sort_values(by=["Ab_absolue"], inplace=True)
table['family'] = table['family'].fillna('Unnammed')
table.to_csv(f"{outputdir}/Genre_ab_Absolue_Relative.csv", index=False)






################################################## Ajouter les informations genre, famille et phylum à la table RefseqFINALE  ########################################




r_input = pd.read_csv(table_ref)
r_input = r_input[(r_input["groupe"] != "Tet__MGE") & (r_input["groupe"] != "Troncated")]

## extraire les informations familles pour chaque genre
dict_GenFam={}

r_taxo = pd.read_csv(f"{outputdir}/Genre_ab_Absolue_Relative.csv", index_col=[0])
for index, row in r_taxo.iterrows():
    genre = index
    family = row["family"]
    phylum = row["phylum"]
    dict_GenFam[genre] = family, phylum 

dict_GenFam


# # créer les colonnes genres familles et phylum dans la table de resumé finale
genus=""
family=""
phylum=""

for index, row in r_input.iterrows() :
    genus_name = row["Bacteria_specie"].split(" ", 1)
    if genus_name[0] in dict_GenFam :
        values = dict_GenFam[genus_name[0]]
        family = values[0]
        phylum = values[1]

        
        r_input.loc[index, "genus"] = genus_name[0]
        r_input.loc[index, "family"] = family
        r_input.loc[index, "phylum"] = phylum



r_input.to_csv(table_ref, index=False)


print("Process Done !!")
