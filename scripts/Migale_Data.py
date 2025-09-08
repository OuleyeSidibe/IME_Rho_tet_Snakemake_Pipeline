# Librairies
import pandas as pd
import sys, os
import fnmatch
import shutil
import glob
import gzip
from argparse import ArgumentParser


##############################################################################################################################################################################
# ## Notes:

""" Download "translated_cds.faa" files in migale database for bacillota and Actinomycetota """

# ## Load script :

# python3 /Migale_data.py -d <directory> -o <output directory> -db /db/gb_bacteria/current/flat/ -out /home/osidibe/work/PPR_MGEproject/snakemakes/screening_1/stdout_report

##############################################################################################################################################################################

## Argument
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-d", "--directory", dest="directory", help="directory")
    parser.add_argument("-o", "--output", dest="outputdir", help="output directory")
    parser.add_argument("-db", '--database', dest="db_repository", help="migale database of genbank bacteria")
    parser.add_argument("-out", "--out", dest="stdout_report", help="stdout report")
    options= parser.parse_args()
    if len(sys.argv) < 4:
        sys.exit("Warning : wrong number of arguments")

    return options.directory, options.outputdir, options.db_repository, options.stdout_report


# main
def main():
    directory, outputdir, db_repository, stdout_report = config_parameters()
    
    os.makedirs(outputdir, exist_ok=True)

    ## Download/open and modified taxonomy file from NCBI
    with open (f"{directory}/NCBI_taxonomie.csv", "r") as taxo_file:
        read_taxo = taxo_file.read()
    caracter_1 = "\t|\t"
    caracter_2 = "\t|"
    replace_caracter1 = "\t"
    replace_caracter2 = ""
    replace_file = read_taxo.replace(caracter_1, replace_caracter1).replace(caracter_2, replace_caracter2)

    with open (f"{directory}/NCBI_taxoResume_1.csv", "w") as taxo_replace:
        taxo_replace.write(replace_file)
    

    ## Only keep taxonomy data from bacillota (frimicutes) et Actinomycetota (actinobacterie)
    read_taxo = pd.read_csv(f"{directory}/NCBI_taxoResume_1.csv", sep="\t", low_memory=False)
    read_taxo = read_taxo.replace('\t|','')
    read_taxo.columns = ["tax_id", "tax_name", "species" , "genus", "family", "order", "class", "phylum", "kingdom", "super_kingdom"]
    read_taxo = read_taxo[(read_taxo["phylum"] == "Bacillota") | (read_taxo["phylum"] == "Actinomycetota")] 
    read_taxo["tax_name"] = read_taxo["tax_name"].str.replace(" ", "_")
    read_taxo.to_csv(f"{directory}/NCBI_taxoResume_2.csv")


    # ## Create à list of NCBI species
    taxe_name = read_taxo["tax_name"].tolist()
    taxe_name_fin =[]
    nb_delete=0
    annoted=0
    annoted_v1=0
    not_annoted=0
    not_migale=0
    in_migale=0
    for name in taxe_name:
        
        ## delete unclassified, uncultured and unidentifies species in NCBI data
        if fnmatch.fnmatch(name, "unclassified_*") or fnmatch.fnmatch(name, "uncultured_*") or fnmatch.fnmatch(name, "unidentified_*"):
            nb_delete+=1
        else :
            taxe_name_fin.append(name)         
    
    # print(taxe_name_fin)
    # ## Exporte uncompressed translated_cds.faa files of NCBI species from migale
    annotedSpecies_list=[]
    phylum=""
    table = pd.DataFrame(columns=["Species", "Phylum", "Number of v1 annotation(s)"]) 
    for species in os.listdir(db_repository):
        
        if species in taxe_name_fin:
            in_migale+=1
                
            GCF_path = os.path.join(db_repository, species + "/latest_assembly_versions")
            
            if os.path.exists(GCF_path):
                
                # file_pattern = "*_translated_cds.faa.gz"
                # GCF_subfile = os.path.join(GCF_path, "/GCF_*.1*", file_pattern)
                CDS_dest = os.path.join(outputdir, species)                
                if not os.path.exists(CDS_dest):
                    os.makedirs(CDS_dest)
                # annotedSpecies_list.append(species)
                for index, row in read_taxo.iterrows():
                    if row['tax_name'] == species : 
                        phylum = row["phylum"]
                        break
                # GCF_dir = glob.glob(f"{GCF_path}/GCF_*.1*/*_translated_cds.faa.gz")
                GCF_dir = glob.glob(f"{GCF_path}/GCF_*.1*")
                for repository in GCF_dir :
                    annoted+=1
                    
                    if ".fasta" not in repository or ".gz" not in repository :
                    
                        GCF_file=glob.glob(f"{repository}/*_translated_cds.faa.gz")
                    
                        for cds in GCF_file :
                            
                            
                            # Uncompress translated_cds.faa.gz in CDS reposytory
                            with gzip.open(cds, "rb") as zip_in:
                                file_accesion = os.path.basename(os.path.dirname(cds))
                                unzip_dest = os.path.join(CDS_dest, file_accesion)

                                with open(unzip_dest, "wb") as zip_out:
                                    shutil.copyfileobj(zip_in, zip_out)
                    else:
                        print(repository)
                    
                # Table of v1 annoted genome number for each specie
                nb_annotation = len(os.listdir(CDS_dest))
                if nb_annotation != 0 :
                    annoted_v1+=1
                    tmp = pd.DataFrame(data=[[species, phylum, nb_annotation]], columns=["Species", "Phylum", "Number of v1 annotation(s)"])                
                    table = pd.concat([table, tmp])

            else :
                not_annoted+=1
            
        else:
            not_migale+=1    
                          
    # nb_species_analysed = len(annotedSpecies_list)
    def genom_by_phylum(table, phylum):
        table_bis = table[table["Phylum"] == phylum]
        nb_annot = table_bis["Number of v1 annotation(s)"].sum()
        return nb_annot
        
    
    # # Number of species in gb_bacteria of migale
    with open(stdout_report, "w") as stdout:
        stdout.write(f"\n=====================================--------Report NCBI & Migale data--------=====================================\n\n")
        stdout.write(f"{len(taxe_name)} NCBI species of Bacillota and Actinomycetota taxonomy.\n\n")
        stdout.write(f"--{nb_delete} excluded data\n")
        stdout.write(f"--{len(taxe_name_fin)} retained data\n\n")
        stdout.write(f"{len(os.listdir(db_repository))} bacteria species in Migale\n")
        
        stdout.write(f"--{in_migale} migale species of Bacillota and Actinomycetota taxonomy\n")
        stdout.write(f"   -> {not_annoted} annotation in progress\n")
        stdout.write(f"   -> {annoted} annoted species data with {annoted_v1} retained species  with version.1 annotation\n\n")
        stdout.write(f"--{not_migale} others bacteria species reported in migale\n\n")
        

        stdout.write(f"** {annoted_v1} retained species with {table['Number of v1 annotation(s)'].sum()} annotations data ** : \n ")
        stdout.write(f" -- {genom_by_phylum(table, 'Bacillota')} Bacillota annotations\n -- {genom_by_phylum(table, 'Actinomycetota')} Actinomycetota annotations")


    table.to_csv(f"{directory}/dataReport_1.csv", sep="\t", index=False)


print("Process Done !!")

if __name__ == "__main__" :
    main()
    ##
    
