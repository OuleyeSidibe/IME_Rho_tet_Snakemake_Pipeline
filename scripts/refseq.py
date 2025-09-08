## Librairies
from argparse import ArgumentParser
import pandas as pd
import glob
import os, sys
import re
import csv
import gzip
import random
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord




#############################################################################################################################################################################
# ## Notes 

""" Summary table of IME_Rho_tet metadata from Refseq/genbank"""

# ## Load script :

# python3 ~/refseq_table.py -i "inputdir" -gb "genbank_dir" -t "Tet(W)"  -m "protein=relaxase" "protein=recombinase" -o "refseq_directory" -db <db_troncSearch.txt>

############################################################################################################################################################################


# Arguments
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="input directory with groupes and gb files")
    parser.add_argument("-gb", "--genbank", dest="gb_dir", help="genbank files directory")
    parser.add_argument("-t", "--tet", dest="tet_protein", help="tet protein in balstp result")
    parser.add_argument("-m", "--MGE", nargs="+", dest="MGE_proteins", help="List of MGE proteins in blastp result")
    parser.add_argument("-o", "--output", dest="outputdir", help="outputdir of migale and result")
    parser.add_argument("-db", "--database", dest="db_troncSearch", help="fasta file of MGE groupe genomic sequence" )
    options=parser.parse_args()
    if len(sys.argv) < 6:
        sys.exit("Warning : wrong number of arguments")
    return options.inputdir, options.gb_dir, options.tet_protein, options.MGE_proteins, options.outputdir, options.db_troncSearch



# Generate unique id for each line of blast_table
def generate_id():
    id_list = []
    new_id = random.randint(100000, 999999)
    while new_id in id_list:
        new_id = random.randint(100000, 999999)
    id_list.append(new_id)
    return new_id



# Dowload Genbank files from migale 
def download_gbfile(gb_dir, specie_name, refseq_acc, outputdir, groupe_name):
    gb_file = os.path.join(gb_dir, specie_name + "/latest_assembly_versions", refseq_acc + f"/{refseq_acc}_genomic.gbff.gz")          
    gb_dest_dir = os.path.join(outputdir + "/GB_files/", groupe_name, specie_name)
    if not os.path.exists(gb_dest_dir):
        os.makedirs(gb_dest_dir)

    gb_dest = os.path.join(gb_dest_dir + f"/{refseq_acc}_genomic.gbff")
    with gzip.open(gb_file, "rb") as zip_in :
        with open(gb_dest, "wb") as zip_out :
            shutil.copyfileobj(zip_in, zip_out)
    return gb_dest



# Split subject of blastp result
def split_subject(subject, cds_position):
    cds_position = str(cds_position)
    cds_position_with_underscore = '_' + cds_position
    #S'assurer de l'echappement des caracteres speciaux et que la chaine se termine aprés la cds position with underescore
    pattern = re.escape(cds_position_with_underscore) + r'$'
    subject = re.sub(pattern, "", subject)
    subject_str = re.split(r'\||\.|prot_', subject)
    return subject_str


# Extracte tet_id and genome accession form subject of blastp result
def genome_acc(blast_table):
    # list_genomAcc=[]
    list_proteinId=[]

    for i in range(0, len(blast_table)):
        subject = blast_table.iloc[i, 0]
        cds_pos = blast_table.iloc[i, 2]
        subject_str = split_subject(subject, cds_pos)
        # genome_acc = subject_str[1]
        protein_id = subject_str[3]
        # list_genomAcc.append(genome_acc)
        list_proteinId.append(protein_id)
    # blast_table["genome_acc"] = list_genomAcc
    blast_table["protein_id"] = list_proteinId
    return blast_table


# For complete genome genbank_file define start_index
def prot_pos_index_gb(prot_position, nb_table, refseq_Annotation, genome_index):
    if nb_table > 1 and refseq_Annotation == "Complet genome" :
        start_index = prot_position
    else :
        start_index = genome_index
        
    return start_index



# Get proteins coordinates and strand 
def get_protein_position(prot_id, lines, start_index) :
    global prot_position
    prot_position_index = None
    #search protein_id in gb file
    for i, line in enumerate(lines[start_index:], start=start_index):
        if prot_id in line :
            # print(line, i)
            prot_position_index = i
            break

    
    # get prot position
    if prot_position_index is not None :
        cds = False
        for i in range(prot_position_index -1, -1, -1):
            if "CDS " in lines[i]: 
                CDS = (lines[i]).split()
                prot_position = CDS[1]
                cds = True
                break
        if not cds :
            for i in range(prot_position_index -1, -1, -1):
                if "gene " in lines[i]: 
                    GENE = (lines[i]).split()
                    prot_position = GENE[1]
                    break
        # print(prot_position)  
        #get protein location on strand
        if "complement" in prot_position :
            strand = "-"
        else:
            strand = "+"
            
        prot_position_str = [int(match.group()) for match in re.finditer(r'\d+', prot_position)]
        start = prot_position_str[0]
        stop = prot_position_str[1]
  
    
    return start, stop, strand, prot_position_index + 20



# Get the rignt coordianates according strand
def relax_start(strand, RELstart, RELstop, RECstart, RECstop):
    if strand == "+" :
        startREL = RELstart
        EGM_coord = str(RELstart) + '-' + str(RECstop)
    else :
        startREL = RELstop
        EGM_coord = str(RECstart) + '-' + str(RELstop)

    return startREL, EGM_coord


def split_gb(line):
    feat=""
    feature = re.findall(r'"(.*?)"', line)
    if len(feature) >= 1:
        feat = feature[0]
    return feat



# main
def main():

    inputdir, gb_dir, tet_protein, MGE_proteins, outputdir, db_troncSearch = config_parameters()

    # Create empty df
    final_table = pd.DataFrame(
        columns=["Unique_Id", "Refseq_accession", "Bacteria_specie", "strain" ,"Annotation", "Sample", "Host", "Country", "Nucleotide_accession", "Size(pb)", "relaxase_id", "relax_identity", "relaxase_coord", "recombinase_id", "recom_identity", "recombinase_coord",
                        "tet_id", "tet_identity", "tet_coord", "tet", "MGE_coord", "strand", "groupe"])


    # # open files of each groupe to get information for table
    # inputdir = "/home/osidibe/work/tetw_project"

    tet_rel_distance = []
    dict_Tet__MGE = {}
    dict_MGE = {}
   

    groupe_repo = glob.glob(f"{inputdir}/g_*")
    for groupe in groupe_repo :

        groupe_name = os.path.basename(groupe)
        print("******GROUPE**********")
        print(groupe_name)

        groupes = groupe_name.split("g_")
        files = glob.glob(f"{groupe}/*/*")
        
        for file in files :
            print(file)
            specie_name = os.path.basename(os.path.dirname(file))

            refseq_file = os.path.basename(file)
            refseq_acc = refseq_file.split(".csv")
            refseq_acc = refseq_acc[0]
            # print(specie_name)
            # print(refseq_acc)


            #Download gb_file from migale if not already dowload
            gb_path = os.path.join(gb_dir, groupe_name, specie_name + f"/{refseq_acc}_genomic.gbff")


            if os.path.exists(gb_path) :
                print("Genbank file already exist !")
            else :
                print("Download Genbank file from migale !")
                # gb_path = download_gbfile(gb_dir, specie_name, refseq_acc, outputdir, groupe_name)

            try : 
                gb_path = download_gbfile(gb_dir, specie_name, refseq_acc, outputdir, groupe_name)
                with open(gb_path, "r") as gb_file :
                    lines = gb_file.read().splitlines()

                # get genome acc and protein id
                read_file = pd.read_csv(file, index_col=[0])
                table = genome_acc(read_file)
                list_genom_acc = table["genome_acc"].unique()

                for genome in list_genom_acc :
                    # print(genome)
                    relax_coord = 0
                    relaxase_id = "-"
                    relax_identity = 0
                    recom_coord = 0
                    recombinase_id = "-"
                    recom_identity = 0
                    tet_coord = 0
                    tet = "-"
                    tet_identity = 0
                    tet_id = "-"
                    EGM_coord = 0
                    TETstart = 0
                    TETstop = 0
                    RELstart = 0
                    RELstop = 0
                    RECstop = 0
                    RECstart = 0
                    start_index = 0
                    genome_index = 0
                    strand=""
                    strand_bis=""
                    strain=None
                    sample=None
                    host=None
                    country=None
                    prot_position_TET=None
                    prot_position_REL=None
                    prot_position_REC=None
                    
                    
                    new_id = generate_id()
                    table_genome = table[table["genome_acc"] == genome]
                    table_genome= table_genome.sort_values(by="cds_position", ascending=True)
                    
                    
                    # get genomic acc size in gbff.gz file
                    line_counter = 1

                    for i, line in enumerate(lines):
                        if line_counter == 1 and genome in line:
                            genome_index = i
                            # print(genome_index)
                            gb_line = line
                            gb_str = gb_line.split()
                            genome_size = gb_str[2]
                            line_counter+=1                
                        
                        elif line_counter != 1 and 'KEYWORDS' in line:
                            line = line.replace(" ", "")
                            line = line.replace("KEYWORDS", "")
                            refseq_Annotation = line[:line.find("RefSeq")]
                            if refseq_Annotation == "" :
                                refseq_Annotation = "Complet genome"
                            line_counter+=1

                        elif line_counter != 1 and '/organism' in line:
                            specie = split_gb(line)

                        elif line_counter != 1 and '/strain' in line:
                            strain = split_gb(line)

                        elif line_counter != 1 and '/isolation_source' in line:
                            sample = split_gb(line)

                        elif line_counter != 1 and '/host' in line:
                            host = split_gb(line)
                    
                        elif line_counter != 1 and '/country' in line:
                            country = split_gb(line)
                            break

                        elif line_counter != 1 and '/note' in line:
                            break


                    table_list = []
                    list_pos = []
                    list_in = []

                    for index, row in table_genome.iterrows():
                        pos = row['cds_position']

                        if index not in list_in: 
                            list_in.append(index)
                            list_pos.append(pos)

                        else :                        
                            table_elt = table_genome[table_genome["cds_position"].isin(list_pos)]
                            table_list.append(table_elt)
                            list_in = []
                            list_pos = []
                            list_in.append(index)
                            list_pos.append(pos)

                    table_elt = table_genome[table_genome["cds_position"].isin(list_pos)]
                    table_list.append(table_elt)
                    
                    nb_table=0
                    for table_x in table_list :
                        
                        nb_table+=1
                        
                        for index, row in table_x.iterrows() :
                            # print(index, row)
                            position = row['cds_position']
                            
                            # Extract relaxase metadata
                            if index == MGE_proteins[0] :
                                relaxase_id = row["protein_id"]
                                relax_identity = row["Pident"]
                                start_index = prot_pos_index_gb(prot_position_REL, nb_table, refseq_Annotation, genome_index)
                                RELstart, RELstop, strand, prot_position_REL = get_protein_position(relaxase_id, lines, start_index)
                                strand_bis=strand
                                relax_coord = str(RELstart) + '..' + str(RELstop)

            
                            # Extract recombinase metadata
                            if index == MGE_proteins[1] :
                                recombinase_id = row["protein_id"]                      
                                recom_identity = row["Pident"]
                                start_index = prot_pos_index_gb(prot_position_REC, nb_table, refseq_Annotation, genome_index)
                                RECstart, RECstop, strand, prot_position_REC = get_protein_position(recombinase_id, lines, start_index)
                                recom_coord = str(RECstart) + '..' + str(RECstop)


                            # Extract tetracycline metadata
                            if index == tet_protein and row["genome_acc"] == genome :
                                tet_id = row["protein_id"]
                                tet = index
                                tet_identity = row["Pident"]
                                # print(nb_table)
                                # print(refseq_Annotation)
                                start_index = prot_pos_index_gb(prot_position_TET, nb_table, refseq_Annotation, genome_index)
                                # print(start_index)
                                TETstart, TETstop, strand, prot_position_TET = get_protein_position(tet_id, lines, start_index)
                                # print(TETstart, TETstop, strand, prot_position_TET , start_index)
                                tet_coord = str(TETstart) + '..' + str(TETstop)

                                
                            
                        # Calculate EGM coordinate and the maximal distance between TETstart and RELstart                    
                        if groupe_name == "g_TetMGE":
                            if strand == "+" :
                                EGM_coord = str(TETstart) + '..' + str(RECstop)
                                distance = int(RELstart) - int(TETstart)
                            else :
                                EGM_coord = str(RECstart) + '..' + str(TETstop)
                                distance = int(TETstop) - int(RELstop)
                            
                            tet_rel_distance.append(distance)
                            # if distance > 6000 :
                            #     # print(distance, specie_name, refseq_acc)
                            #     print(gb_path)
                                



                        # get distance between start genome and RELstart for groupe Tet_MGE
                        elif groupe_name == "g_Tet__MGE" :
                            
                            if RELstart != 0 and TETstart != 0 :
                                strand=strand_bis
                                startREL, EGM_coord = relax_start(strand, RELstart, RELstop, RECstart, RECstop)
                                genome_Tet__MGE = genome
                                dict_Tet__MGE[genome_Tet__MGE] = [[specie_name, refseq_acc, startREL, groupes[1]]]
                                
                            elif RELstart != 0 and TETstart == 0 :
                                startREL, EGM_coord = relax_start(strand, RELstart, RELstop, RECstart, RECstop)
                                tet_id = "-"
                                genome_Tet__MGE = genome
                                dict_Tet__MGE[genome_Tet__MGE] = [[specie_name, refseq_acc, startREL, groupes[1]]]
                                
                            elif RELstart == 0 and TETstart != 0 :
                                relax_identity = 0
                                recom_identity = 0
                                EGM_coord = str(TETstart) + '..' + str(TETstop)
                                
                                

                        # Check a possible troncated tet PPR gene in MGE group
                        elif groupe_name == "g_MGE" :
                            startREL, EGM_coord = relax_start(strand, RELstart, RELstop, RECstart, RECstop)

                            
                            if dict_MGE.get(genome) is not None :
                                dict_MGE[genome] += [[specie_name, refseq_acc, startREL, groupes[1]]]
                            else :
                                dict_MGE[genome] = [[specie_name, refseq_acc, startREL, groupes[1]]]
                            


                        elif groupe_name == "g_Tet" or tet != "-" :
                            EGM_coord = str(TETstart) + '..' + str(TETstop)
                                

                        result = pd.DataFrame(
                        columns=["Unique_Id", "Refseq_accession", "Bacteria_specie", "strain" ,"Annotation", "Sample", "Host", "Country", "Nucleotide_accession", "Size(pb)", "relaxase_id", "relax_identity", "relaxase_coord", "recombinase_id", "recom_identity", "recombinase_coord",
                        "tet_id", "tet_identity", "tet_coord", "tet", "MGE_coord", "strand", "groupe"], data = [[new_id, refseq_acc, specie_name, strain, refseq_Annotation, sample, host, country, genome, genome_size, relaxase_id, relax_identity, relax_coord, recombinase_id, recom_identity,
                            recom_coord, tet_id, tet_identity, tet_coord, tet, EGM_coord, strand, groupes[1]]])
                    
                        final_table = pd.concat([final_table, result], ignore_index=True)
                        final_table.to_csv(f"{inputdir}/Refseq_1.csv")
                        # print(final_table)
                    
            except FileNotFoundError:
                print(f"{gb_path} not found!!")
                pass
                
    max_distance = max(tet_rel_distance)
    print(max_distance)

    # print(dict_MGE)
    # print("****dict de base***")
    # print(dict_Tet__MGE)

    # ####################################################### Reorganize group in the summary table #############################################################
    grp2="MGE"
    grp3="Tet__MGE"

    #Fonction pour attribuer les bons groupes dans la table
    def arrange_tableGrp(dict_grp, max_distance, final_table, grp2, grp3):
        list_key = []
        for key, val in dict_grp.items() :
            for value in dict_grp[key] :
                if (value[3] == grp2 and int(value[2]) < max_distance) or (value[3] == grp3 and int(value[2]) > max_distance):
                    index_line = (final_table[final_table["Nucleotide_accession"] == key].index) 
                    

                    if value[3] == grp2 :                    
                        final_table.loc[index_line, "groupe"] = grp3  
                                        
                    else :
                        final_table.loc[index_line, "groupe"] = grp2
                        # print(index_line)
                        if (final_table.loc[index_line, "tet_coord"]).any() == 0 :
                            final_table = final_table.drop(index_line+1)
                        
                    list_key.append(key)
        
        return final_table, list_key
                    

    # list_key_MGE : list des genomes d'accessions à transferer du groupe MGE au groupe Tet_MGE
    table2, list_key_MGE = arrange_tableGrp(dict_MGE, max_distance, final_table, grp2, grp3)
    # list_key_TMGE : list des genomes d'accessions à transferer du groupe Tet_MGE au groupe MGE
    result_final, list_key_TMGE = arrange_tableGrp(dict_Tet__MGE, max_distance, table2, grp2, grp3)

    result_final.to_csv(f"{inputdir}/Refseq_2.csv")

#    print(list_key_MGE)
#    print("****list key***")
#    print(list_key_TMGE)

    # ######################################################## Reorganize group in repository #############################################################

    #Create new dictionnary of moved species or genomes
    def transfert_dict(dict_aceptor, dict_donor, list_key_donor):
        dict_aceptor={}
        for key in list_key_donor :
            if key in list(dict_donor):
                values = dict_donor[key]
                dict_aceptor[key] = values
                del dict_donor[key]

        return dict_aceptor

  
    to_T_MGE_dict = {}
    to_MGE_dict = {}
    to_T_MGE_dict = transfert_dict(to_T_MGE_dict, dict_MGE, list_key_MGE)
    to_MGE_dict = transfert_dict(to_MGE_dict, dict_Tet__MGE, list_key_TMGE)

    # print(to_T_MGE_dict)
    # print("********to_dict******************")
    # print(to_MGE_dict)

    # Function for creted directory
    def create_dir(inputdir, grp_name, specie_name):
        dest_path = os.path.join(inputdir + "/g_"+ grp_name, specie_name)
        os.makedirs(dest_path, exist_ok=True)
        return dest_path

    # Transfert genome in correct repository
    def right_grp_files(dict_grp):
        with open(f"{inputdir}/refseq_report_3.txt", "a") as log:
            log.write(f"{dict_grp}\n\n")
            nb_copy=0
            nb_move=0
            
            for key, val in dict_grp.items() :
                for value in dict_grp[key] :

                    # Moving file
                    file_path = os.path.join(inputdir + "/g_" + value[3], value[0], value[1] + ".csv")
                    read_filepath = pd.read_csv(file_path, index_col=[0])
                    table_filepath = genome_acc(read_filepath)
                    
                    
                    MGE_prot = 0
                    MGE_pos = []
                    
                    # Read moving file to export only concern proteins in acceptor group
                    for index, row in table_filepath.iterrows() :

                        if index in MGE_proteins and row["genome_acc"] == key :
                            MGE_pos.append(row["cds_position"])
                            MGE_prot += 1
                            

                        if MGE_prot == 2 :
                            old_file = pd.read_csv(file_path, index_col=[0])
                            
                            # MGE to move in new repository
                            move_file = old_file[old_file["cds_position"].isin(MGE_pos)]
                            # save moved file in tmp repository
                            source_dest = os.path.join(inputdir, value[1] + ".csv")

                            move_file.to_csv(source_dest)
                            # define destination path of moved file
                            if value[3] == grp2 :
                                dest_grp = grp3
                                dest_path = create_dir(inputdir, dest_grp, value[0])
                            elif value[3] == grp3 :
                                dest_grp = grp2
                                dest_path = create_dir(inputdir, dest_grp, value[0])
                            
                            
                            # move file in new repository groupe
                            # dest_pathF = os.path.join(dest_path, value[0] + "/")
                            file_dest = os.path.join(dest_path, value[1] + ".csv")
                            if os.path.exists(file_dest) :
                                with open(file_dest, "a", newline="") as dst_file, open(source_dest, "r") as src_file:
                                    dst_writer = csv.writer(dst_file)
                                    src_reader = csv.reader(src_file)
                                    #Copy src row in dst file en ignorant l'entête
                                    next(src_reader)                                  
                                    for row in src_reader :
                                        dst_writer.writerow(row)
                                nb_copy+=1
                                # delete tmp move file
                                os.remove(source_dest)
                            else : 
                                shutil.move(source_dest, dest_path)
                                print(dest_path)
                                nb_move+=1
                                
                            log.write(f"{value[0]}/{value[1]} from : {value[3]} to {dest_grp}\n")
                            
                            
                            
                            # Si le nouveau fichier ne comporte pas un autre EGM, le supprimer
                            new_file = old_file[~old_file["cds_position"].isin(MGE_pos)]
                            if len(new_file) < 2 :

                                os.remove(file_path)
                                log.write(f"{file_path} removed in origin group !\n")
                            else :
                                new_file.to_csv(f"{file_path}")
                                


                            # Delete empty repository
                            dir_path = os.path.join(inputdir + "/g_" + value[3], value[0])
                            reposit = os.listdir(dir_path)
                            if not reposit :
                                os.rmdir(dir_path)
                                log.write(f"{value[0]} removed because empty repository !\n")
                            MGE_prot = 0
                            log.write("----------------\n\n")
                            break
            log.write(f"{nb_copy} data copied in existing repository and {nb_move} data moved in new reposytories !!\n\n")
                    
                        
        return 
    
    # # Reorganise repository
    MGE_to_TET_MGE = right_grp_files(to_T_MGE_dict)
    Tet_MGE_to_MGE = right_grp_files(to_MGE_dict)
        


    # ######################################################## Search troncated tetracycline #############################################################


    # Extract contig or genome fasta sequence from genomic_fna file in migale db 
    tmp_file=os.path.join(inputdir + "/tmp")
    concat_dict = {**dict_Tet__MGE, **to_T_MGE_dict}
    nb_key=0

    for key, val in concat_dict.items():
        
        for value in concat_dict[key] :
            nb_key+=1
            fasta_sequence = ''
            genomicfna_file = os.path.join(gb_dir, value[0] + "/latest_assembly_versions", value[1] + f"/{value[1]}_genomic.fna.gz")

            with gzip.open(genomicfna_file, "rb") as genomic_fna_gz :
                lines = genomic_fna_gz.read().splitlines()
                for i, line in enumerate(lines) :
                    if key.encode('utf-8') in line :
                        # print(line)
                        fasta_sequence += line.decode('utf-8') + '\n'
                        startFasta = i+1
                        # print(startFasta)
                        break

                if startFasta is not None :
                    # print(startFasta)
                    for i in range(startFasta, len(lines)) :
                        line = lines[i]
                        
                        if not '>'.encode('utf-8') in line:                                          
                            fasta_sequence += line.decode("utf-8") + '\n'
                        else :
                            break


            # extract maximal distance sequence
            stop_fa = int(value[2])
            start_fa = stop_fa - max_distance
            if start_fa < 1 :
                start_fa = 1
                   
            with open(tmp_file, "w+") as tmp_fasta :
                tmp_fasta.write(fasta_sequence)
                
                tmp_fasta.flush()
                tmp_fasta.seek(0)
                    
                for record in SeqIO.parse(tmp_fasta, 'fasta') :
                    subseq=record.seq[start_fa-1:stop_fa]                                                                                                                             
                    # new_record=record[start_fa-1:stop_fa]
                    # new_record.id=f"{record.id}_{start_fa}_{stop_fa}"
                    # new_record.description=f"{record.description}_[{start_fa}:{stop_fa}]"
                    new_record= SeqRecord(
                        subseq,
                        description=f"{record.description}_[{start_fa}:{stop_fa}]",
                        id=""
                    )
            
                    with open (db_troncSearch, 'a') as db_file :
                        SeqIO.write(new_record, db_file, 'fasta')
        
print("Process Done !!")
if __name__ == "__main__" :
    main()
    
    
    
