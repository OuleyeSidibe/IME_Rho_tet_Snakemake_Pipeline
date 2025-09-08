## Librairies
from argparse import ArgumentParser
from Bio import SeqIO
import sys



################################################################################################################################################################
# ## Notes

""" Create query file of unreferenced tet proteins clusters """

# ## Load script :

# python3 ~/create_query.py -i "inputfile" -o "outputfile"

#################################################################################################################################################################



# Arguments
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputdir", help="input directory with clusters of proteines")
    parser.add_argument("-o", "--out", dest="outputdir", help="output file report")
    options=parser.parse_args()
    if len(sys.argv) < 2:
        sys.exit("Warning : wrong number of arguments")
    return options.inputdir, options.outputdir


inputdir, outputdir = config_parameters()

# Create query file for blastp analysis
def create_query(input_file, output_file):
    with open (input_file, "r") as newref, open(output_file, "a") as query_file :
        for record in SeqIO.parse(newref, "fasta"):
            if "WP_" in record.id :
                SeqIO.write(record, query_file, "fasta")
    return query_file


create_query(inputdir, outputdir)

