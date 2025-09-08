#!/bin/bash

#################################################################################################################################
# ## Notes 

"""  Search for TIR motifs in all genomes selected by the variable FILTER, using the fimo algorithm from the meme suite """

# Load script : 
#  1 - launch script twice, one for normal species name and then for unamed species 'sp. XXX'
#  2 - same script is launched for RPP (tet) and IME_Rho_tet (TetMGE) replacing in FILTER argument

################################################################################################################################


# check numbers of arguments
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <directory> <inputdir> <outputdir> <refseqTable>"
  exit 1
fi

directory=$1    
inputdir=$2     
outputdir=$3   
refseqTable=$4 


# filter data about groupe and analyse table line by line
FILTER="Tet"
IFS=$'\n'


# conda activate meme-5.5.7

#search in csv file line which contains filter, extract interest columns and stock data in vals
# 1st launch
for vals in `grep $FILTER $refseqTable | cut -d "," -f 2,3,9,21,22`; do

# 2nd launch
# for vals in `grep $FILTER $refseqTable | grep " sp. " | cut -d "," -f 2,3,9,21,22`; do # process only unamed species: 'sp. XXX' 

 #extract accession genome and species name
 ass=`echo $vals | cut -d ',' -f 1`
 # 1st launch
 species=`echo $vals | cut -d ',' -f 2| cut -d " " -f 1,2 | tr " " "_"`  
 
 # 2nd launch
#  species=`echo $vals | cut -d ',' -f 2| tr " " "_"`  # without the 'cut -d " " -f 1,2' to process unamed species: 'sp. XXX' 
 

# Launch alone (not with fimo )if necessary t dezip assembli files fisrt if necessary

#   if [ -e "$inputdir/g_Tet/$species/$ass/${ass}_genomic.fna.gz" ]; then
#   echo dezip process of $ass in $species
#   gzip -d "$inputdir/g_Tet/$species/$ass/${ass}_genomic.fna.gz"

#   else
#     if [ -e "$inputdir/g_Tet/$species/$ass/${ass}_genomic.fna" ]; then
#     echo $ass in $species is already dezip




 ##### then run the fimo loop

 #check if genomic sequence files of genomes available in migale data
 if [ -e "$inputdir/g_Tet/$species/$ass/${ass}_genomic.fna" ] ; then
  echo process $ass in $species
  fimo --text --thresh 1e-9 --no-qvalue --oc $outputdir/results_Tet3 $directory/TIR_model.meme $inputdir/g_Tet/$species/$ass/${ass}_genomic.fna > $outputdir/results_Tet3/$ass.tsv 
  
 else
  if ! [[ -e "$inputdir/g_Tet/$species/$ass/${ass}_genomic.fna" ]];then

   echo fna for $ass does not exist in $species >> $outputdir/run_fimo.out


  fi
 fi
done 
