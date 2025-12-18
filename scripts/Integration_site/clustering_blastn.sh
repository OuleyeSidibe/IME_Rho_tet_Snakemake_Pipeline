#!/bin/bash
# This script is used to run BLASTP on a set of protein sequences against the NCBI non-redundant protein database (nr).
# It uses the BLAST+ toolkit installed via conda.

inputfile="path/of/reconstruct/file.fa"
output_clstr="repository/of/clustering_out"
output_blast="/home/osidibe/work/PPR_MGEproject/integration_site/1000G/blastn_out.txt"


#activate conda environment if necessary

cd-hit -i $inputfile -o $output_clstr -c 0.9 -d 0 -sc 1 -g 1  


csplit -z -f seq_ -b "%03d.fasta" $output_clstr '/^>/' '{*}' 
for f in seq_*.fasta; do
  blastn -query "$f" -remote -db nr -evalue 1e-5 -max_target_seqs 10 -outfmt "6 qacc sacc qcovs pident evalue sstart send sstrand stitle" >> $output_blast
done

