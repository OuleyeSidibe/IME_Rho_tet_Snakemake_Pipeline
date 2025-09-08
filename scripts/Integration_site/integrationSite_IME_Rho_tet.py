# Librairies
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import os
import gzip


##########################################################################################
# ## Notes

""" TIRs flanking regions extraction of IME_Rho_tet-carring genomes """

## Load script :

# python3 /integrationSite_IME_Rho_tet.py -i <inputdir> -o <outputdir> -f <fa_dir>

########################################################################################


# Arguments
def config_parameters():
    parser=ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputfile", help=" IME_Rho_tet table with TIRs information")
    parser.add_argument("-o", "--output", dest="outputdir", help="output directory")
    parser.add_argument("-f", "--fa", dest="fa_dir", help="assemblies directory")
    args=parser.parse_args()
    if len(sys.argv) < 3 :
        sys.exit("Warning : wrong number of argument")
    return args.inputfile, args.outputdir, args.fa_dir



inputfile, outputdir, fa_dir = config_parameters()

# read table
df= pd.read_csv(inputfile, usecols=["Refseq_accession","Bacteria_specie","Nucleotide_accession","TIR_amont_coord","TIR_aval_coord","strand","groupe"], index_col=[2])


# extraction parameters
extr_Pb = 150
max_lenght = (extr_Pb - 2)*2 # TSD


# Extract nucleotide sequence for accessions
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
                
                #extract fasta files
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
                            
                            
                            # check if sequence is enough long
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

# read table line by line
for index, row in df.iterrows():
    specie_name = row["Bacteria_specie"].replace(" ", "_")
    refseq_acc = row["Refseq_accession"]
    nucl_acc = index
    groupe_name = row['groupe']
    TIR_amont_coord = row['TIR_amont_coord']
    TIR_aval_coord =row['TIR_aval_coord']
    strand= row['strand']

    download_fa_file(fa_dir, specie_name, refseq_acc, outputdir, groupe_name, nucl_acc, TIR_amont_coord,TIR_aval_coord,strand)
    
