###############################################
# Extract source and host metadata for the ~65K genomes of the original dataset

# 1. on MIGALE -> extract metadata for each genome
./1_collect_sourcehost_forall.sh > sourcehost.out 2> missing_assemblies.out

# 2. categorize according to identified keywords
cat sourcehost.out| ./2_categorize.sh > sourcehost.categorized.out

# 3. return distribution gut / non gut / unknown
cat sourcehost.categorized.out | cut -f 3 | sort | uniq -c

# 4. return distribution host type
cat sourcehost.categorized.out | grep -P "\tgut" | cut -f 4 | sort | uniq -c


###############################################
# Extract source and host metadata for the IME_Rho_tet-carrying genomes

# 1. extract source and host for each genome
cat  ../RefseqFINAL5_TIRS_IME_Rho_tet.csv | awk 'BEGIN{FS=","}{print $6"\t"$7"\t"$1 }'> IME_Rho_tet_sourcehost.out

# 2. categorize according to identified keywords
cat IME_Rho_tet_sourcehost.out| ./2_categorize.sh > IME_Rho_tet_sourcehost.categorized.out

# 3. Merge manually the source_host and the original genome table. 


###############################################
# Extract source and host metadata for the tet-carrying genomes

# 1. extract source and host for each genome
cat  ../RefseqFINAL5_TIRS_tet.csv | awk 'BEGIN{FS=","}{print $6"\t"$7"\t"$1 }'> tet_sourcehost.out

# 2. categorize according to identified keywords
cat tet_sourcehost.out| ./2_categorize.sh > tet_sourcehost.categorized.out

# 3. Merge manually the source_host and the original genome table. 

