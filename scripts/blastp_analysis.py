## Librairies
import pandas as pd
from argparse import ArgumentParser
import glob
import os
import sys
import re
import fnmatch



###########################################################################################################################################################################
# ## Notes

""" Processing of blast outputs to identify ARGs close to relaxase and recombinase in Bacillota and Actinomycetota bacterial genomes  """

# ## Load script :

# python3 -i <inputdir with blast result> -o <outputdir> -pi <min of identity percentage> -pc <min of coverage percentage> -t <list of "tet protein" launch in blastp>
# -int <cds interval> -a <stdout.txt>

###########################################################################################################################################################################

# Arguments
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="input directory with blast data")
    parser.add_argument("-o", "--output", dest="outputdir", help="output directory for groupe")
    parser.add_argument("-pi", "--Pident", type=int, dest="Pidentity", help="min identity Percentage")
    parser.add_argument("-pc", "--Pcov", type=int, dest="Pcoverage", help="min coverage Pourcentage ")
    parser.add_argument("-t", "--tet", dest="tet_protein", help="tetracycline protein in query file")
    parser.add_argument("-int", "--interval", type=int, dest="interval", help="cds interval")
    parser.add_argument("-a", "--append", dest="stdout", help="stdout file")
    options = parser.parse_args()
    if len(sys.argv) < 7 :
        sys.exit("Warning : wrong number of arguments")
    return options.inputdir, options.outputdir, options.Pidentity, options.Pcoverage, options.tet_protein, options.interval, options.stdout 



## Functions 
dict_tetMGE={}
dict_tetOnly={}


# Function for extract proteine cds position
def prot_position(proteine, csv_file):
    prot = csv_file[csv_file.index == proteine]
    prot_position = prot["cds_position"].tolist()
    return prot_position
    

# Function for count number of genome for each specie
def nb_genome(path_dir):
    nb_genome=0
    path_rep = os.path.dirname(path_dir)
    for specie in os.listdir(path_rep):
        nb_genome+=len(os.listdir(os.path.join(path_rep, specie))) 
    nb_genome = nb_genome -1 
    return nb_genome


# Function for creted directory
def create_dir(outputdir, dir_name, specie_name):
    dir_path = os.path.join(outputdir, dir_name, specie_name)
    os.makedirs(dir_path, exist_ok=True)
    return dir_path


# Function for Tet proteines count in differents groupes
def tet_value(tet, dict_tet):
    if tet in dict_tet:
        dict_tet[tet] += 1
    else:
        dict_tet[tet] = 1
    return dict_tet


# Function to concate table if we have more than on MGE or tetMGE in genome
def concate_data(repository, table1, table2, tet_pos=None):
    
    if tet_pos is None :  
        if repository is None :
            repository = table1
        else :
            repository = pd.concat([repository, table1], axis=0)
    else :
        if repository is None : 
            repository = pd.concat([table1, table2[table2["cds_position"] == tet_pos]], axis=0)
        else :
            repository = pd.concat([repository, table1, table2[table2["cds_position"] == tet_pos]], axis=0)
    return repository



# Main
def main():

    inputdir, outputdir, Pidentity, Pcoverage, tet_protein, interval, stdout = config_parameters()

    dossier_blastp = glob.glob(f"{inputdir}/*/*")
    nb_file = 0
    # blastpF_path=""

    for file in dossier_blastp:
        print(file)
        tetMGE = None
        MGE = None
        tet_MGE = None
        nb_file += 1
        
        specie_name = os.path.basename(os.path.dirname(file))
        
        # print(specie_name)
        refseq_acc= os.path.basename(file)
        # print(refseq_acc)
        read_blastp= pd.read_csv(f"{inputdir}/{specie_name}/{refseq_acc}", header=None, usecols=[0,1,2,3], index_col=[0], names=["Query", "Subject", "Pcov", "Pident"])
 
        """ Filter on Pident and Pcov and stock data in new repository"""
        
        new_file = read_blastp[read_blastp["Pident"]>=Pidentity] 
        filter_file = new_file[new_file["Pcov"]>=Pcoverage]

        """ Extract CDS position """ 

        sbj_acces = filter_file["Subject"]

        list_cdsPosition=[]
        list_genome_acc=[]
        for subj in sbj_acces:
            position = subj.split("_")
            genome_acc = re.split(r'\||\.', subj)
            list_genome_acc.append(genome_acc[1])
            list_cdsPosition.append(position[-1])

        filter_file.loc[:, "cds_position"] = list_cdsPosition
        filter_file.loc[:, "genome_acc"] = list_genome_acc
        filter_file.drop("Pcov", axis=1, inplace=True)
        

        # Remove duplicate positions
        file_sort = filter_file.sort_values(by="cds_position", ascending=True)
        # table = file_sort.drop_duplicates(subset=["cds_position"], keep='first')

        table = file_sort 
        
        # Ne pas sauvegarder les fichiers vides
        if len(table) != 0 :
            blastpF_path = create_dir(outputdir, "Blastp_filter", specie_name)
            print(blastpF_path)
            table_bis = table
            table_bis.to_csv(f"{blastpF_path}/{refseq_acc}.csv")

            """ check closeness relaxase and recombinase genes in genomes """

            if (table_bis.index == "protein=relaxase").any() and (table_bis.index == "protein=recombinase").any():
                
                # CDS positions of relaxase, recombinase and tet

                for pos in prot_position("protein=relaxase", table_bis):
                    
                    # Create un intervale of 10 cds up and downstream of relaxase cds
                    cds_interval = range(-interval + int(pos), interval + int(pos))
                    

                    # check if recombinase cds is in relaxase intervalle cds
                    for pos_2 in prot_position("protein=recombinase", table_bis) :                    
                        if ((int(pos_2) > int(pos) and int(pos_2) <= cds_interval.stop -1)) or ((int(pos_2) < int(pos) and int(pos_2) >= cds_interval.start - 1)):
                            n_indexRel=table_bis.reset_index().index[table_bis["cds_position"] == pos]
                            n_indexRec=table_bis.reset_index().index[table_bis["cds_position"] == pos_2]
                            genomeRel = table_bis.iloc[n_indexRel,3].tolist()
                            genomeRec = table_bis.iloc[n_indexRec, 3].tolist()

                            if genomeRel[0] == genomeRec[0]:

                                rex_recom = pd.concat([table_bis[table_bis["cds_position"] == pos], table_bis[table_bis["cds_position"]== pos_2]], axis=0) 


                                # Check if tet genes belong to MGE
                                if (tet_protein == table_bis.index).any():
                                    
                                    tet_pos = prot_position(tet_protein, table_bis)                
                                    for postet in tet_pos:                                    

                                        # Si l'MGE complet associé au tet dans le brin sens ou antisens
                                        if ( ((int(pos_2) > int(pos) and int(pos_2) <= cds_interval.stop -1) and (int(postet) >= cds_interval.start and int(postet) < int(pos))) or ((int(pos_2) < int(pos) and int(pos_2) >= cds_interval.start - 1) and (int(postet) <= cds_interval.stop and int(postet) > int(pos))) ):                              
                                            
                                            # Check number of each Tet protein with element

                                            dict_tet_MGE = tet_value(tet_protein, dict_tetMGE)
                                            tetMGE = concate_data (tetMGE, rex_recom, table_bis, postet)

                                            TetMGE_path = create_dir(outputdir, "g_TetMGE", specie_name)
                                            # tetMGE = tetMGE.sort_values(by="cds_position", ascending=True)
                                            tetMGE.to_csv(f"{TetMGE_path}/{refseq_acc}.csv")

                                            # Drop MGE and tet in variable
                                            table_bis = table_bis[~table_bis["cds_position"].isin((postet, pos, pos_2))]
                                            break
                                    break

            

                ##  If tet in genomes but not with MGE check if it's in upstream of MGE
                
                # nb_groupeMGE = 0
                for pos_bis in prot_position("protein=relaxase", table_bis):
                    cds_interval = range(-interval + int(pos_bis), interval + int(pos_bis))
                    
                    nb_MGE=0
                    for pos_bis_2 in prot_position("protein=recombinase", table_bis):
                        if ((int(pos_bis_2) > int(pos_bis) and int(pos_bis_2) <= cds_interval.stop -1)) or ((int(pos_bis_2) < int(pos_bis) and int(pos_bis_2) >= cds_interval.start - 1)):
                            n_indexRel=table_bis.reset_index().index[table_bis["cds_position"] == pos_bis]
                            n_indexRec=table_bis.reset_index().index[table_bis["cds_position"] == pos_bis_2]
                            genomeRel = table_bis.iloc[n_indexRel,3].tolist()
                            genomeRec = table_bis.iloc[n_indexRec, 3].tolist()

                            if genomeRel[0] == genomeRec[0]:
                                rex_recom = pd.concat([table_bis[table_bis["cds_position"] == pos_bis], table_bis[table_bis["cds_position"]== pos_bis_2]], axis=0)

                                exit_loop = False
                                nb_tet = 0
                                nb_MGE +=1 
                                
                                if (tet_protein == table_bis.index).any():
                                    
                                    for position in prot_position(tet_protein, table_bis) :
                                        if ( ((int(pos_bis_2) > int(pos_bis) and int(pos_bis_2) <= cds_interval.stop -1) and (int(position) < cds_interval.start - 1)) or ((int(pos_bis_2) < int(pos_bis) and int(pos_bis_2) >= cds_interval.start - 1) and (int(position) > cds_interval.stop - 1))) :
                                            nb_tet+=1
                                            
                                            t_rex_recom = pd.concat([rex_recom, table_bis[table_bis["cds_position"] == position]], axis=0)
                                            Tet_path = create_dir(outputdir, "g_Tet__MGE", specie_name)
                                            t_rex_recom.to_csv(f"{Tet_path}/{refseq_acc}.csv")
                                            
                                            # Drop MGE and tet in variable
                                            table_bis = table_bis[~table_bis["cds_position"].isin((postet, pos_bis, pos_bis_2))]
                                            exit_loop = True
                                    
                                    if exit_loop :
                                        break
                                    
                                    
                                    # Create MGE groupe if no tet in genomes
                                    if nb_tet == 0 and nb_MGE == 1:
                                        tet_pos = None
                                        MGE_path = create_dir(outputdir, "g_MGE", specie_name)
                                        MGE = concate_data (MGE, rex_recom, table_bis)
                                        MGE.to_csv(f"{MGE_path}/{refseq_acc}.csv")
                                        table_bis = table_bis[~table_bis["cds_position"].isin((pos_bis, pos_bis_2))]
                                        break
                                        
                                else : 

                                    MGE_path = create_dir(outputdir, "g_MGE", specie_name)
                                    rex_recom.to_csv(f"{MGE_path}/{refseq_acc}.csv")
                                    table_bis = table_bis[~table_bis["cds_position"].isin((pos_bis, pos_bis_2))]
                                

            elif (table_bis.index != "protein=relaxase").any() and (table_bis.index != "protein=recombinase").any():
                table_bis_tet = pd.DataFrame()

    
                if tet_protein in table_bis.index :
                    table_bis_tet = pd.concat([table_bis_tet, table_bis[table_bis.index == tet_protein]], axis=0)
                    dict_tet_Only = tet_value(tet_protein, dict_tetOnly)

                onlyTET_path = create_dir(outputdir, "g_Tet", specie_name)
                table_bis_tet.to_csv(f"{onlyTET_path}/{refseq_acc}.csv")


            elif (table_bis.index == "protein=relaxase").all() or (table_bis.index != "protein=recombinase").all():
                relax_or_recom_path = create_dir(outputdir, "gOut", specie_name)
                table_bis.to_csv(f"{relax_or_recom_path}/{refseq_acc}.csv")


    # Number of genomes with/without tet genes with recombinase and relaxase    
    with open(stdout, "a") as stdout :
        stdout.write(f"\n\n=====================================--------Report blastp & groupes repartition--------=====================================\n\n")
        stdout.write(f"-- {nb_file} blastp launch genomes -- :  {nb_genome(blastpF_path)} retained genomes with threshold coverage ({Pcoverage}%) and identity ({Pidentity}%)\n\n")
        stdout.write(f"   --->>> {nb_genome(TetMGE_path)} MGE carring PPR coding gene : g_TetMGE \n") 
        stdout.write(f"   --->>> {nb_genome(MGE_path)} MGE not carring PPR coding gene : g_MGE\n")
        stdout.write(f"   --->>> {nb_genome(Tet_path)} MGE and PPR coding gene in different contig/scaffold : g_Tet__MGE\n")
        stdout.write(f"   --->>> {nb_genome(onlyTET_path)} only PPR coding gene : g_Tet\n")
        stdout.write(f"   --->>> {nb_genome(relax_or_recom_path)} out groupe : gOut")
        

print("Process Done !!")

if __name__ == "__main__":
    main()
