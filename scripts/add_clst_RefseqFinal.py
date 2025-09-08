## Libraries
from argparse import ArgumentParser
import pandas as pd
import os, re
import glob
import sys



################################################################################################################################################
# ## Infos 

""" Add cluster results in RefseqFinal table"""
# ## Lancement script :

# python3 ~/add_clstr_RefseqFinal.py -i "inputdir" -o "outputdire" -t "refseq_table for step blastp_groupe_2" -g <"MGE" "TetMGE" "Tet">


###############################################################################################################################################



# Analyseur d'arguments
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="input directory with clusters of proteines")
    parser.add_argument("-o", "--out", dest="outputdir", help="output directory")
    parser.add_argument("-t", "--table", dest="table", help="RefSeq finale table without pseudogenes data")
    parser.add_argument("-g", "--groupe", nargs="+", dest="groupes", help="list of groupe names")
    
    options=parser.parse_args()
    if len(sys.argv) < 4:
        sys.exit("Warning : wrong number of arguments")
    return options.inputdir, options.outputdir, options.table, options.groupes




# definir les proteines correspondantes suivant les groupes
def prot_by_grp(grp):
    
    proteins=["relaxase_id", "recombinase_id", "tet_id" ]
    
    if grp == "MGE":
        prot = [proteins[0], proteins[1]]
    elif grp == "TetMGE" :
        prot = [proteins[0], proteins[1], proteins[2]]
    else :
        prot = [proteins[2]]
        
    return prot






# Recupérer le numero de cluster des protéines correspondants
def get_clstr(proteins, grp, table, clstr_path):

    dict_clstr={}
    df = table[table["groupe"] ==  grp]    
    for prot in proteins :
        #define prot name
        print(prot)
        name = prot[:3].upper()

        list_i = df[prot].tolist()
        print(list_i)

        """ get cluster number """
        
        patern_i = os.path.join(clstr_path, f"{name}*.clstr")
        path = glob.glob(patern_i)
        for i in path:
            path_i = i
        
        for id in list_i :
            with open(path_i, "r") as clstr_file :
                lines = clstr_file.read().splitlines()
                
                position=None
                for i, line in enumerate(lines):
                    if id in line :
                        position = i
                        break
                    
                if position is not None :
                    cluster = False
                    
                    for i in range(position-1, -1,-1):
                        if "Cluster" in lines[i] :
                            line_cltr = lines[i].split()
                            nb_clstr = line_cltr[1]
                            nb_clstr = f"{name}_C{nb_clstr}"
                            cluster = True
                            break
            
            
            """ get reference cluster """
            
            patern_y = os.path.join(clstr_path, f"{name}*newRef.fa")
            path = glob.glob(patern_y)
            
            for y in path:
                path_y = y

            with open(path_y, "r") as ref_file :
                lines = ref_file.read().splitlines()
                # print(lines)
                for line in lines :
                    if f"C{line_cltr[1]}." in line :
                        # print("***")
                        ref_clstr = line.split(">")
                        # print(ref_clstr)
                        break
                        
                        
            if dict_clstr.get(prot) is not None :
                dict_clstr[prot] += [[id, nb_clstr, ref_clstr[1]]]
            else:
                dict_clstr[prot] = [[id, nb_clstr, ref_clstr[1]]]
                
                
    return dict_clstr



## main
def main():
    
    inputdir, outputdir, table, groupes = config_parameters()
    print(groupes)
    
    # Read refseq table
    read_table = pd.read_csv(table, index_col=[0])
    read_table = read_table.sort_values(by="groupe")


    
    for grp in groupes :
        
        grp = grp.split("_",1)
        grp = grp[1]
        proteins = prot_by_grp(grp)
        proteins = prot_by_grp(grp)

        dictClstr= get_clstr(proteins, grp, read_table, inputdir) 
        
        #remplir la table refseq des données de clustering
        for key, value in dictClstr.items() :
            for val in dictClstr[key] :
                col_name = key[:3].upper()
                col_clstr = f"{col_name}_clstr"
                col_ref = f"{col_name}_ref"
                
                for index, row in read_table.iterrows():
                    if row[key] == val[0] and row["groupe"] == grp:
                        
                        read_table.loc[index, col_clstr] = val[1]
                        read_table.loc[index, col_ref] = val[2]
                
    read_table.to_csv(outputdir, index=False)    

print("Process Done !!")

if __name__ == "__main__" :
    main()
