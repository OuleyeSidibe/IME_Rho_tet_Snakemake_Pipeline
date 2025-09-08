#!/bin/bash


inputdir=$1
threshold=$2

# conda activate cd-hit-4.8.1

for prot in "$inputdir"/*.fa; do

    file_name=$(basename "$prot")
    out_name=${file_name%.*}_$threshold

    # clusterisation at thresholt param each protein db file
    cd-hit -i $prot -o $inputdir/$out_name -c $threshold -d 0 -sc 1 -g 1
    
    

done
    

