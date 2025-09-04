from Bio import Entrez, SeqIO
import pandas as pd


file = "/home/osidibe/work/PPR_MGEproject/integration_site/nucl_analyse/blastn_out_filtered_sorted_97.csv"
filtered_df_sort = pd.read_csv(file, sep=',', index_col=1)

# Ajout des colonnes de résultat
filtered_df_sort["refseqGB_acc"]=""
filtered_df_sort['CDS_product'] = ""
filtered_df_sort['CDS_locus_tag'] = ""
filtered_df_sort['CDS_protein_id'] = ""
filtered_df_sort['CDS_location'] = ""

Entrez.email = "ouleye.sidibe@inrae.fr"

integration_offset = 150
tolerance = 10

list_not_found_NZ = []
records_cache = {}

for index, row in filtered_df_sort.iterrows():
    
    accession = f"NZ_{index}"
    
    s_start = row['s_start']
    strand = row['strand']

    
    # voir si pour une même accession il y a plusieurs CDS qui couvrent le site
    if row['CDS_product'] != "" or pd.isna(row["CDS_product"]):
        print(f"⚠️ Attention : plusieurs CDS couvrent le site {int_site} pour {accession}.")  
        print(f"⚠️ CDS déjà trouvé : {filtered_df_sort.at[index, 'CDS_locus_tag']} ({filtered_df_sort.at[index, 'CDS_product']})")
            
    
    # Correction du type de strandbastp du
    if accession not in records_cache:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
            record = SeqIO.read(handle, "genbank")
            records_cache[accession] = record
            handle.close()
        except Exception as e:
            print(f"❌ Erreur pour {accession} : {e}")
            list_not_found_NZ.append(index) 
            continue

    record = records_cache[accession]


    # Calcul de la position du site d'intégration
    int_site = s_start + integration_offset if strand == "plus" else s_start - integration_offset

    # Recherche du gène couvrant la position du site d'intégration
    found = False
    for feature in record.features:
        
        if feature.type in ["CDS"]:
            start = int(feature.location.start)
            end = int(feature.location.end)

            if start - tolerance <= int_site <= end + tolerance:
                
                
                # Récupération des informations du CDS
                product = feature.qualifiers.get("product", ["N/A"])[0]
                locus = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                location_str = f"{start}..{end}"

                filtered_df_sort.at[index, 'refseqGB_acc'] = accession
                filtered_df_sort.at[index, 'CDS_product'] = product
                filtered_df_sort.at[index, 'CDS_locus_tag'] = locus
                filtered_df_sort.at[index, 'CDS_protein_id'] = protein_id
                filtered_df_sort.at[index, 'CDS_location'] = location_str

                print(f"✅ Site {int_site} couvert par : {start}..{end} : ({product})")
                found = True
                break

    if not found:
        print("❌ Aucun CDS trouvé pour cette position.")
        list_not_found_NZ.append(index)

# Export du fichier modifié
output_path = "/home/osidibe/work/PPR_MGEproject/integration_site/nucl_analyse/blastn_out_filtered_sorted_refseqGB97.csv"
filtered_df_sort.to_csv(output_path)



print("############################## same script testing refseq accesion with 'NC' #########################")

list_not_found_NC = []
filtered_df_sort = pd.read_csv(output_path, sep=',', index_col=0)
for index, row in filtered_df_sort.iterrows():
    
    if index in list_not_found_NZ and pd.isna(row["CDS_product"]) or row['CDS_product'] == "":

    
        accession = f"NC_{index}"
        
        s_start = row['s_start']
        strand = row['strand']

        
        # Correction du type de strandbastp du
        if accession not in records_cache:
            try:
                # handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
                record = SeqIO.read(handle, "genbank")
                records_cache[accession] = record
                handle.close()
            except Exception as e:
                print(f"❌ Erreur pour {accession} : {e}")
                list_not_found_NC.append(index) 
                continue

        record = records_cache[accession]


        # Calcul de la position du site d'intégration
        int_site = s_start + integration_offset if strand == "plus" else s_start - integration_offset
        print(f"🔍 Position estimée du site : {int_site} (+/- {tolerance}pb)")
        

        # Recherche du gène couvrant la position du site d'intégration
        found = False
        for feature in record.features:
            
            if feature.type in ["CDS"]:
                start = int(feature.location.start)
                end = int(feature.location.end)

                if start - tolerance <= int_site <= end + tolerance:

                    # Récupération des informations du CDS
                    product = feature.qualifiers.get("product", ["N/A"])[0]
                    locus = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                    location_str = f"{start}..{end}"

                    filtered_df_sort.at[index, 'refseqGB_acc'] = accession
                    filtered_df_sort.at[index, 'CDS_product'] = product
                    filtered_df_sort.at[index, 'CDS_locus_tag'] = locus
                    filtered_df_sort.at[index, 'CDS_protein_id'] = protein_id
                    filtered_df_sort.at[index, 'CDS_location'] = location_str

                    print(f"✅ Site {int_site} couvert par : {start}..{end} : ({product})")
                    found = True
                    break

        if not found:
            print("❌ Aucun CDS trouvé pour cette position.")
            list_not_found_NC.append(index)

# Export du fichier modifié
output_path = "/home/osidibe/work/PPR_MGEproject/integration_site/nucl_analyse/blastn_out_filtered_sorted_refseqGB97.csv"
filtered_df_sort.to_csv(output_path)



print("############################## same script testing gb accesion with not found accesion in refseq ########################")

list_not_found_NZ_NC_gb = []
filtered_df_sort = pd.read_csv(output_path, sep=',', index_col=0)
for index, row in filtered_df_sort.iterrows():
    
    if index in list_not_found_NC and pd.isna(row["CDS_product"]) or row['CDS_product'] == "":

        accession = index
        print(f"🔍 Recherche pour l'accession : {accession}")
        
        s_start = row['s_start']
        strand = row['strand']

        
        # Correction du type de strandbastp du
        if accession not in records_cache:
            try:
                # handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
                record = SeqIO.read(handle, "genbank")
                records_cache[accession] = record
                handle.close()
            except Exception as e:
                print(f"❌ Erreur pour {accession} : {e}")
                list_not_found_NZ_NC_gb.append(index) 
                continue

        record = records_cache[accession]


        # Calcul de la position du site d'intégration
        int_site = s_start + integration_offset if strand == "plus" else s_start - integration_offset
        print(f"🔍 Position estimée du site : {int_site} (+/- {tolerance}pb)")
        

        # Recherche du gène couvrant la position du site d'intégration
        found = False
        for feature in record.features:
            
            if feature.type in ["CDS"]:
                start = int(feature.location.start)
                end = int(feature.location.end)

                if start - tolerance <= int_site <= end + tolerance:
                    
                    # Récupération des informations du CDS
                    product = feature.qualifiers.get("product", ["N/A"])[0]
                    locus = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                    location_str = f"{start}..{end}"

                    filtered_df_sort.at[index, 'refseqGB_acc'] = accession
                    filtered_df_sort.at[index, 'CDS_product'] = product
                    filtered_df_sort.at[index, 'CDS_locus_tag'] = locus
                    filtered_df_sort.at[index, 'CDS_protein_id'] = protein_id
                    filtered_df_sort.at[index, 'CDS_location'] = location_str

                    print(f"✅ Site {int_site} couvert par : {start}..{end} : ({product})")
                    found = True
                    break

        if not found:
            print("❌ Aucun CDS trouvé pour cette position.")
            list_not_found_NZ_NC_gb.append(index)

# Export du fichier modifié
output_path = "/home/osidibe/work/PPR_MGEproject/integration_site/nucl_analyse/blastn_out_filtered_sorted_refseqGB97.csv"
filtered_df_sort.to_csv(output_path)
print(f"CDS not found : {list_not_found_NZ_NC_gb}")
