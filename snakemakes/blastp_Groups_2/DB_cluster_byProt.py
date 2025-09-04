## Librairies 
from argparse import ArgumentParser
from Bio import SeqIO
import pandas as pd 
import glob
import os, sys, re 
import fnmatch



###############################################################################################################################################################################################
# ## Notes

""" Create proteines database for clustering_by_prot step"""


## Load script :

# python3 /script/path/DB_cluster.py -i <inputdir> -o <outputdir>  -f >faa_dir> -p "TetW" "Tet32" "protein=relaxase" "protein=recombinase" -g <groups in order : TetMGE, MGE, Tet> -r <resfinder fasta file>

# """Les identifiants des fichiers fasta r1, r2, r3 doivent imperativemment commencaient par "Tet", "Rel" et "Rec" """ exple : >Tet(W)_1_DQ060146, >Relaxase_IME_RhoA2 ,>Recombinase_IME_RhoA2 

################################################################################################################################################################################################



# Analyseurs d'arguments
def config_parameters():
    parser=ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="inputdir groups for protéins db")
    parser.add_argument("-o", "--output", dest="outputdir", help="outputdir for proteins groups of cluster")
    parser.add_argument("-f", "--faa", dest="faa_dir", help="inputdir of faa_CDS files")
    parser.add_argument("-p", "--prot", nargs="+", dest="proteins", help="index proteins in data")
    parser.add_argument("-g", "--groups", nargs="+" ,dest="groups", help="groupes of prot to clusterise")
    parser.add_argument("-r1", "--resfinder", dest="tet_ref", help="fasta file of Tet PPR  and mozaÏc on resfinder DB")
    parser.add_argument("-r2", "--relaxase", dest="rel_ref", help="fasta file of relaxase protein in IME_RhoA2")
    parser.add_argument("-r3", "--recombinase", dest="rec_ref", help="fasta file of recombinase protein in IME_RhoA2")
    args=parser.parse_args()
    if len(sys.argv) < 8 :
        sys.exit("Warning : wrong number of argument")
    return args.inputdir, args.outputdir, args.faa_dir, args.proteins, args.groups, args.tet_ref, args.rel_ref, args.rec_ref


inputdir, outputdir, faa_dir, proteins, groups, tet_ref, rel_ref, rec_ref = config_parameters()


## Function to create database file of each proteins
def create_db (outputdir, line, prot):

    #create reposytory if necessary
    os.makedirs(outputdir, exist_ok=True)
    db_path = os.path.join(outputdir + f"/{prot}_db.fa")
    
    # Open for ID protein first write or append
    if not os.path.exists(db_path):
        with open(db_path, "w") as db_fasta :
            db_fasta.write(f"{line}\n")
    else :
        with open(db_path, "a") as db_fasta :
            db_fasta.write(f"{line}\n")
    return db_path



## Put fasta sequence in database file for each proteins
def fill_in_db (outputdir, faa_lines, prot_id, prot):
    # Add ID line of prot
    path=""
    start_fasta = None
    for i, line in enumerate(faa_lines):
        if prot_id in line :
            path = create_db(outputdir, line, prot)
            start_fasta = i
            break
        # Add sequence of prot
    if start_fasta is not None :
        for i in range(start_fasta +1, len(faa_lines)) :
            line = faa_lines[i]
            if not '>' in line :
                with open(path, "a") as db_fasta :
                    db_fasta.write(f"{line}\n")                    
            else :
                break

    return path


# extracte nucleotide accession
def nucl_accession(record_id):
    ID_ = record_id.rsplit("_", 1)
    ID_str = re.split(r'\||_prot', ID_[0]) 
    nucl_acc = ID_str[1]
    return nucl_acc




## ADD proteines reference in database file for each protein
def add_reference(prot_path, reference_file, outputdir):
    psd_id=""
    psd_seq=""
    psd_fa = ""
    nb_pseudo = 0
    nb_fa = 0 
    rm_nucl = []
    if prot_path :
        with open(reference_file, "r") as ref_fa, open(prot_path, "a") as prot_db:
            for record in SeqIO.parse(ref_fa, "fasta"):
                SeqIO.write(record, prot_db, "fasta")
            
        #report pseudo sequence in stdout file
        prot = os.path.basename(prot_path).split(".fa")
        # grp = os.path.basename(os.path.dirname(prot_path))
        with open (prot_path, "r") as prot_db2,  open(f"{outputdir}/report_3.stdout", "a") as fa_stdout :
            for record in SeqIO.parse(prot_db2, "fasta"):
                if record.id :
                    nb_fa +=1 
                    
                if "*" in record.seq :
                    rm_nucl.append(nucl_accession(record.id))

                    #ignored * in seq
                    nb_pseudo += 1
                    psd_id = f">{record.description}" 
                    psd_seq = record.seq
                    psd_fa += f"{psd_id}\n{psd_seq}\n"
                
                from Bio import SeqIO
# Delete psseudo gene in corresponding prot db
def delete_pseudo(path, rm_record):
    
    if path :
        with open (path, "r") as prot_db, open (f"{outputdir}tmp.txt", "w") as prot_db2 :
            
            
            for record in SeqIO.parse(prot_db, "fasta"):
                nb_nucl=0
                for nucl in rm_record :
                    if nucl not in record.id:
                        nb_nucl+=1
                if nb_nucl == len(rm_record) :
                    print(nb_nucl)
                    SeqIO.write(record, prot_db2, 'fasta')
            
            os.replace(f"{outputdir}tmp.txt", path)


## mean

for grp in groups :

    grp_prot = glob.glob(f"{inputdir}/{grp}/*/*")

    for file in grp_prot :
        REL_ID = []
        REC_ID = []
        TET_ID = []
        specie_name = os.path.basename(os.path.dirname(file))
        refseq_file = os.path.basename(file)
        path_file = os.path.join(inputdir, grp, specie_name, refseq_file)
        read_file = pd.read_csv(path_file, index_col=[0])

        # for each protein extract ID protein
        for index, row in read_file.iterrows() :
            if index == proteins[0] :
                TET_ID.append(row["Subject"])

            elif index == proteins[1] :
                REL_ID.append(row["Subject"])

            elif index == proteins[2] :
                REC_ID.append(row["Subject"])
    

        #Search ID proteine of protein in faa file
        refseq_acc = refseq_file.split(".csv")
        faa_path = os.path.join(faa_dir, specie_name, refseq_acc[0])
        
        with open(faa_path, "r") as faa_fasta :
            faa_lines = faa_fasta.read().splitlines()
            
            # for each protein fill in the create database file
            for tet_id in TET_ID :
                tet_path=fill_in_db(outputdir, faa_lines, tet_id, "TET")

            for rel_id in REL_ID :
                REL_path =fill_in_db(outputdir, faa_lines, rel_id, "REL")
                
            for rec_id in REC_ID :
                REC_path =fill_in_db(outputdir, faa_lines, rec_id, "REC")
                

# Add resfinder PPR reference in tet database and RhoA2-183 for relaxase and recombinase for clustering
TET_path, rm_TET = add_reference(tet_path, tet_ref, outputdir)
REL_path, rm_REL = add_reference(REL_path, rel_ref, outputdir)
REC_path, rm_REC = add_reference(REC_path, rec_ref, outputdir)

rm_record = rm_TET + rm_REL + rm_REC
rm_record = list(set(rm_record))
# print(rm_record)


# delete pseudo gene in corrresponding file
delete_pseudo(TET_path, rm_record)
delete_pseudo(REL_path, rm_record)
delete_pseudo(REC_path, rm_record)


print("Process Done !!")