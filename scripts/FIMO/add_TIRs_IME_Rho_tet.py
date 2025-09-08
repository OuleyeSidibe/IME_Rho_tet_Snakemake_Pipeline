## Librairies
import pandas as pd
from argparse import ArgumentParser


#######################################################################################
# ## Notes

""" Add TIRS of IME_Rho_tet data in summary table """

## Load script :

# python3 /add_TIRS_IME_Rho_tet.py -i <inputfile> -o <output table> -t <summary table>

#######################################################################################

# Arguments
def config_parameters():
    parser=ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputfile", help="TIRs RPP data")
    parser.add_argument("-o", "--output", dest="outputfile", help="summary table with TIRS RPP data")
    parser.add_argument("-t", "--table", dest="table", help="RefSeq finale table without pseudogenes data")
    args=parser.parse_args()
    if len(sys.argv) < 3 :
        sys.exit("Warning : wrong number of argument")
    return args.inputfile, args.outputfile, args.table


inputfile, outputfile, table = config_parameters()


#1 read refseq table to add information about TIRs
df = pd.read_csv(table, sep=",", index_col=0)
df = df[df["groupe"] == "TetMGE"]

# add new columns
df["TIR_amont"]= ""
df["TIR_amont_coord"]= ""
df["TIR_aval"]= ""  
df["TIR_aval_coord"]= ""
df["MGE_size"]=""

for index, row in df.iterrows():
    uniq_id = index
    # read the TIRs table to check uniq_id
    df_TIRs = pd.read_csv(inputfile, sep=",", index_col=0)
    #4 check if the uniq_id is in the TIRs table
    if uniq_id in df_TIRs.index and row["groupe"] == "TetMGE":

            
        #5 get the TIRs information
        tir_amont = df_TIRs.loc[index, "TIRamont_seq"]
        
        tir_aval = df_TIRs.loc[index, "TIRaval_seq"]
        
        if tir_amont is not None or tir_aval is not None :
            
            # print(df_TIRs.loc[index, "TIRamont_s"])
            tir_amont_s = df_TIRs.loc[index, "TIRamont_s"]
            tir_amont_e = df_TIRs.loc[index, "TIRamont_e"]
            if pd.notna(tir_amont_s) and pd.notna(tir_amont_e):
                tir_amont_coord = f"{int(tir_amont_s)}..{int(tir_amont_e)}"
                # print(tir_amont_coord)
            else :
                tir_amont_coord = ""
                
                
            tir_aval_s = df_TIRs.loc[index, "TIRaval_s"]
            tir_aval_e = df_TIRs.loc[index, "TIRaval_e"]
            if pd.notna(tir_aval_s) and pd.notna(tir_aval_e):
        
                tir_aval_coord = f"{int(tir_aval_s)}..{int(tir_aval_e)}"
            else :
                tir_aval_coord = ""
            

            #6 add the TIRs information to the dataframe
            df.at[index, "TIR_amont"] = tir_amont
            df.at[index, "TIR_aval"] = tir_aval
            df.at[index, "TIR_amont_coord"] = tir_amont_coord
            df.at[index, "TIR_aval_coord"] = tir_aval_coord
            
       

df.to_csv(outputfile, index=FALSE)
