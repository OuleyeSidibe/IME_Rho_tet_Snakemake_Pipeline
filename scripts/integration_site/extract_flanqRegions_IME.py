# Librairies
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# # from Bio.SeqRecord import SeqRecord
import os
import gzip

## paths
path="/home/osidibe/work/PPR_MGEproject/Fimo/fimo_newID/RefseqFINAL3_TIRS_IME_Rho_tet.csv"

# fa_dir = "/db/gb_bacteria/gb_bacteria_2025-03-16/flat"
fa_dir="/work_projet/isp-pgba/sebastien/tetPPR/Assemblies/g_TetMGE/"
outputdir = "/home/osidibe/work/PPR_MGEproject/integration_site/nucl_analyse_newID"
# lire la table et les colonnes necessaires
df= pd.read_csv(path, usecols=["Refseq_accession","Bacteria_specie","Nucleotide_accession","TIR_amont_coord","TIR_aval_coord","strand","groupe"], index_col=[2])


# parametre d'extension des palindrome en amont et en aval
extr_Pb = 150
max_lenght = (extr_Pb - 2)*2 # TSD

# test avec les 5 premiére lignes
# df = df.iloc[:5,]


# extraire la sequence fasta de l'accession nucléotidique
def download_fa_file(fa_dir, specie_name, refseq_acc, outputdir, groupe_name, nucl_acc, TIR_amont_coord,TIR_aval_coord,strand):
    fa_file = os.path.join(fa_dir, specie_name , refseq_acc + f"/{refseq_acc}_genomic.fna")          
    fa_dest_dir = os.path.join(outputdir + "/FA_files/", groupe_name, specie_name)
    if not os.path.exists(fa_dest_dir):
        os.makedirs(fa_dest_dir)

    fa_dest = os.path.join(fa_dest_dir + f"/{refseq_acc}_genomic.fna")
    
    try:
        with open(fa_file, "rt") as zip_in :
            
            mode ="a" if os.path.exists(fa_dest) else "w" 
            with open(fa_dest, f"{mode}t") as zip_out, open(f"{outputdir}/integ_fasta.fa", "a") as new_fa :
                
                #extraire uniquemment le fasta correspondant
                for record in SeqIO.parse(zip_in, 'fasta'):
                    if nucl_acc in record.id :
                        # get tir coord and #extract news seqs
                        if pd.notna(TIR_amont_coord) and pd.notna(TIR_aval_coord) :
                            Tir_amont = TIR_amont_coord.split("..")
                            Tir_aval = TIR_aval_coord.split("..")
                            
                            if strand == "+" :
                                SeqAmont_start = int(Tir_amont[0]) - extr_Pb
                                SeqAmont_stop = int(Tir_amont[0]) -2 ## 2pb of TSD
                                
                                SeqAval_start = int(Tir_aval[1])
                                SeqAval_stop = int(Tir_aval[1]) + extr_Pb
                                
                                new_seq_amont = record.seq[SeqAmont_start:SeqAmont_stop]
                                new_seq_aval = record.seq[SeqAval_start:SeqAval_stop]
                                
                                new_seq = new_seq_amont + new_seq_aval
                            else :
                                SeqAval_start = int(Tir_aval[0])- extr_Pb
                                SeqAval_stop = int(Tir_aval[0]) - 2 ## 2pb of TSD
                                
                                SeqAmont_start = int(Tir_amont[1])
                                SeqAmont_stop = int(Tir_amont[1]) + extr_Pb
                                
                                new_seq_aval = record.seq[SeqAval_start:SeqAval_stop]
                                new_seq_amont = record.seq[SeqAmont_start:SeqAmont_stop]
                                
                                new_seq = new_seq_aval + new_seq_amont
                                new_seq = new_seq.reverse_complement()
                            
                            
                            # check if seq is enough long
                            if len(new_seq) < max_lenght :
                                print(f"Seq amont or aval too short : {nucl_acc} amont: {len(new_seq_amont)}, aval: {len(new_seq_aval)}, strand: {strand}")
                                #do not take sequence
                                
                            else :
                                new_record = record
                                SeqIO.write(new_record, zip_out, "fasta")
                                new_fa.write(f">{nucl_acc}\n{new_seq}\n")
                    
                    
    except PermissionError:
        print(f"Permission denied : {fa_file}")   

    except FileNotFoundError:
        print(f"File not found : {fa_file}")

# Parcourir la table ligne par ligne
for index, row in df.iterrows():
    specie_name = row["Bacteria_specie"].replace(" ", "_")
    refseq_acc = row["Refseq_accession"]
    nucl_acc = index
    groupe_name = row['groupe']
    TIR_amont_coord = row['TIR_amont_coord']
    TIR_aval_coord =row['TIR_aval_coord']
    strand= row['strand']

    download_fa_file(fa_dir, specie_name, refseq_acc, outputdir, groupe_name, nucl_acc, TIR_amont_coord,TIR_aval_coord,strand)
    
