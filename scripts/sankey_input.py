# ## Librairies
from argparse import ArgumentParser
import pandas as pd 
import sys 



#######################################################################################
# ## Notes

""" create sankey input file for clustering by proteins """
# ## Load script :

# python3 ~/sankey_input.py -i "inputdir" -o "outputdire" -g  <"MGE" "TetMGE" "Tet">

#######################################################################################

# Arguments
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="Refseq table final")
    parser.add_argument("-o", "--out", dest="outputdir", help="outputdir of outputfile")
    parser.add_argument("-g", "--groupe", nargs="+", dest="groupes", help="list of groupe names")

    options=parser.parse_args()
    if len(sys.argv) < 3:
        sys.exit("Warning : wrong number of arguments")
    return options.inputdir, options.outputdir, options.groupes


# main
def main() :
    
    inputdir, outputdir, groupes = config_parameters()
    
    # refseq input table : choice of groupe 
    r_input = pd.read_csv(inputdir)
    r_input = r_input[(r_input["groupe"] != "Tet__MGE") & (r_input["groupe"] != "Troncated")]


    """ 1st groupbing by genus  """
        
    ## replace  (Nan) by "Na"
    df= r_input.fillna("Na")
    grouped_gg = df.groupby(by=["TET_ref","REL_clstr", "REC_clstr", "groupe", "genus"]).size().reset_index(name="nb_genome")
    ## replace "Na" by originate values
    grouped_gg = grouped_gg.replace("Na", pd.NA)


    """ 2nd grouping of the gender table according to the three proteins, summing the number of occurrences for each cluster link. """
    
    grouped_gg_ = grouped_gg.fillna("Na")
    grouped_ggref = grouped_gg_.groupby(by=["TET_ref","REL_clstr", "REC_clstr", "groupe"]).agg({"nb_genome" : 'sum'}).reset_index()
    grouped_ggref = grouped_ggref .replace("Na", pd.NA)
    grouped_ggref.to_csv(f"{outputdir}/grouped_tableRef.csv")

    # Open output file for writting
    with open(f"{outputdir}/sankey_input.txt", "a") as out_file :
        
        
        dicts={}
        for grp in groupes:
            out_file.write(f"---------------------{grp}---------------------\n")    
            
            # store each row of the table in a tuple       
            tup_clstr=[]
            grouped_ggref_ = grouped_ggref[grouped_ggref["groupe"] == grp]
            for index, row in grouped_ggref_.iterrows() :

                tup_clstr.append(tuple(row))

            # check tuple occurence of grouped_ggref in grouped_gg table to get genus and their occurence
            for tup in tup_clstr :
                nb_seq = 0

                grouped_gg_ = grouped_gg[grouped_gg["groupe"] == grp]
                for index, row in grouped_gg_.iterrows():

                    row_tup = tuple(row)

                    if tup[0:4] == row_tup[0:4]:
                        nb_seq += row["nb_genome"]

                        if dicts.get(tup) is not None :
                            dicts[tup] += [row_tup[-2:]]
                        else:
                            dicts[tup] = [row_tup[-2:]]
            


        # Filter na values
        for key_na, value in dicts.items():
            key = tuple(item for item in key_na if pd.notna(item))

            # create link betweeen prot nodes and group
            link = len(key) - 2
            nb_link = 0
            for i in range(0, link+1):
                if nb_link != link :
                    y = i +1
                    out_file.write(f"{key[i]} [{key[-1]}] {key[y]}\n")
                    nb_link +=1 
                    

print("Process Done !!")

if __name__ == "__main__" :
    main()

