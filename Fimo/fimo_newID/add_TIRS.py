import pandas as pd

table="/home/osidibe/work/PPR_MGEproject/Fimo/fimo_newID/RefseqFINAL3.csv"

table_TIRs_tet="/home/osidibe/work/PPR_MGEproject/Fimo/fimo_newID/parse_outTet.csv"
# table_TIRs_IME_Rho_tet="/home/osidibe/work/PPR_MGEproject/Fimo/fimo_newID/parse_outIME_Rho_tet.csv"
#1 read refseq table to add information about TIRs
df = pd.read_csv(table, sep=",", index_col=0)
df = df[df["groupe"] == "Tet"]

#2 ajouter les colonnes pour les TIRs
df["TIR_amont"]= ""
df["TIR_amont_coord"]= ""
df["TIR_aval"]= ""  
df["TIR_aval_coord"]= ""
df["MGE_size"]=""

for index, row in df.iterrows():
    uniq_id = index
    # read the TIRs table to check uniq_id
    df_TIRs = pd.read_csv(table_TIRs_tet, sep=",", index_col=0)
    #4 check if the uniq_id is in the TIRs table
    if uniq_id in df_TIRs.index and row["groupe"] == "Tet":

    # read the TIRs table to check uniq_ID
    # df_TIRs = pd.read_csv(table_TIRs_seb, sep=",", index_col=0)
    
    # for index_2, row_2 in df_TIRs.iterrows():
    #     #4 check if the uniq_id is in the TIRs table
        
    #     if index == index_2 and row["Nucleotide_accession"] == row_2["Nucleotide_accession"]  :
            
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
            
       

# print(df)
# df.to_csv("/home/osidibe/work/PPR_MGEproject/Fimo/fimo_newID/RefseqFINAL3_TIRS_IME_Rho_tet.csv")
df.to_csv("/home/osidibe/work/PPR_MGEproject/Fimo/fimo_newID/RefseqFINAL3_TIRS_Tet.csv")