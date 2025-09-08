
 ##Librairies
import pandas as pd 
import glob
import os, re, sys
from argparse import ArgumentParser
from Bio import SeqIO



################################################################################################################################################################
# ## Infos 

""" Sortir les chiffres de clustering par proteines (nb clsusters, seq, genres bacteriens) et définir les références Resfinder et RhoA2-183 comme référence des clusters """

# ## Lancement script :

# python3 ~/analyse_cluster.py -i "inputdir" -t "refseq_table" -p <prot_name> -o >outputdir> -th >threshold>



#################################################################################################################################################################



# Analyseur d'arguments
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="input directory with clusters of proteines")
    parser.add_argument("-t", "--table", dest="refseq_table", help="refseq summary table")
    parser.add_argument("-p", "--prot", nargs="+",dest="prot_name", help="proteine_name of cluster")
    parser.add_argument("-o", "--out", dest="outputdir", help="output file report")
    parser.add_argument("-th", "--theshold", dest="threshold", help="clustering threshold")
    options=parser.parse_args()
    if len(sys.argv) < 5:
        sys.exit("Warning : wrong number of arguments")
    return options.inputdir, options.refseq_table, options.prot_name, options.outputdir, options.threshold



# Extract reference ID of cluster
def cluster_ref(line):
    if "lcl" in line :
        ref= re.split(r'\.\.\.|\|', line)
        ref= ref[1]
    else :
        ref= re.split(r'\.\.\.|\>', line)
        ref= ref[1]
    return ref



# Extract protein and nucleotide accession
def nucl_accession(id_seq):
    ID_ = id_seq.rsplit("_", 1)
    ID_str = re.split(r'\||\.|prot_', ID_[0]) 
    nucl_acc = ID_str[0]
    prot_acc = ID_str[2]

    return nucl_acc, prot_acc



# Define genus bacteria of proteins cluster 
def bacteria(refseq_t, nucl_acc):
    df = pd.read_csv(refseq_t, index_col=[0])
    bact_genus= None
    
    for index, row in df.iterrows():
        if row["Nucleotide_accession"] == nucl_acc:
            bact_genus = df.loc[index, "Bacteria_specie"]
            
            bact_genus = bact_genus.split(" ")
            bact_genus = bact_genus[0]
    if bact_genus is None :
        print(nucl_acc)

    return bact_genus
            
 
 
## Main 
def main():
    
    inputdir, refseq_table, prot_name, outputdir, threshold = config_parameters()            


    with open(outputdir, "a") as stdout_file :
        

        stdout_file.write("===============================================================================================================================================================\n")
    
        for prot in prot_name :
            
            prot = f"{prot}_db_{threshold}"
            prot_dir = os.path.join(inputdir, prot)

            # change cdhit reference by resfinder or RhoA2-183 reference
            dict_clustSize = {}
            
            # outfile stat of clustering
            clstr_info = os.path.join(prot_dir + ".clstr")
            
            # outfile of clusters reference
            clstr_ref = prot_dir
            new_clstr_ref = os.path.join(clstr_ref + "newRef.fa")
            
            # fasta file of all sequences prot before clustering-
            prot_id = prot.split('_', 1)
            faa_db = os.path.join(inputdir, prot_id[0] + "_db.fa")



            if os.path.exists(clstr_info):

                # open cluster information file to calcule some parameters
                with open(clstr_info, "r") as out_clst :
                    nb_clust=0
                    nb_seq=0
                    

                    lines = out_clst.read().splitlines()
                    for i, line in enumerate(lines) :
    
                        size=0
                        seq_list=[]
                        genus_list=[]
                        genus_dict={}
                        
                        # define cluster number and number of sequence in cluster
                        if '>Cluster' in line :
                            cluster_name = line
                            cluster = line.split(" ")
                            nb_clust+=1
                            for i in range (i+1, len(lines)):
                                line = lines[i]
                                
                                # define number of sequence by cluster and cluster reference
                                if not '>Cluster' in line and not "*" in line:
                                    size+=1
                                    seq = cluster_ref(line)
                                    seq_list.append(seq)
                                elif not '>Cluster' in line and "*" in line :
                                    size+=1
                                    ref = cluster_ref(line)
                                    seq_list.append(ref)
                                else :
                                    break

                            # Make reference prot as same reference cluster if possible
                            seq2_list = sorted(seq_list, key=len)
                            nb_genome=0
                            nb_ref = 0
                            nb_x = 0
                            for x in seq2_list:
                                if x.startswith(prot_id[0].capitalize()) and nb_ref == 0:
                                    ref = x.replace("(", "").replace(")", "")
                                    nb_ref+=1
                                # define bacteria genus of each seq 
                                if not x.startswith(prot_id[0].capitalize()):
                                    nucl_acc, proteine = nucl_accession(x)
                                    bact_genus = bacteria(refseq_table, nucl_acc)
                                    nb_genome +=1 
                                    nb_x += 1
                                    
                                    if bact_genus not in genus_list :
                                        genus_dict[bact_genus] = 1
                                    else:
                                        genus_dict[bact_genus] += 1
                                    genus_list.append(bact_genus)
                                    

                            ## Create a new file of referenced prot as clusters reference
                            with open(new_clstr_ref, "a") as new_file, open(faa_db, 'r') as fa_db :
                                
                                for record in SeqIO.parse(fa_db, "fasta"):
                                    id_record = record.id.replace("(", "").replace(")", "")
                                                                            
                                    if ref in record.id and "protein" in record.description:

                                        ID_spl = record.description.split(" ")
                                        ID_rpl = ID_spl[0].replace("lcl|", "")
                                        genome, prot_acc = nucl_accession(ID_rpl)
                                        id2 = prot_acc

                                        for x in ID_spl :


                                            if "gene=" in x and prot_id[0] == prot_name[0] :                                            
                                                seq_id1 = re.split(r"=|\]", x)


                                            elif ("protein=" in x and prot_id[0] != prot_name[0]) or ("protein=" in x and prot_id[0] == prot_name[0]):
                                                seq_id1 = x.split("=")
                                                

                                        id1 = seq_id1[1].replace("(", "").replace(")", "")
                                        if ("rec" in id1 or "rel" in id1) or ("Rec" in id1 or "Rel" in id1 ):
                                            id1 = id1[:3]
                                        

                                        id_fin = f"{id1}__{id2}_C{cluster[1]}."
                                        seq_fa = record.seq
                                            # print(id)

                                    elif ref in id_record:
                                        id_fin = f"{id_record}_C{cluster[1]}."
                                        
                                        seq_fa = record.seq
                                if nb_genome != 0 :
                                        
                                    new_file.write(f">{id_fin}\n{seq_fa}\n")
                            dict_clustSize[cluster_name] = [size, id_fin, nb_genome, genus_dict]

                        else:
                            nb_seq +=1

                        
                stdout_file.write(f"{prot_id[0]}_{threshold} : {nb_seq} Sequence(s), {nb_clust} Cluster(s)\n\n")
                for key, value in dict_clustSize.items():
                    if value[2] != 0 :
                        stdout_file.write(f"{key} : {value[0]} Sequence(s), {value[1]}, {value[2]} Genome(s) : {value[3]}\n")
                # print("\n")
                # print("\n")
                stdout_file.write("\n--------------------------------------------------------------------------------------------------------------\n\n")

                        

print("Process Done !!")
if __name__ == "__main__" :
    main()
