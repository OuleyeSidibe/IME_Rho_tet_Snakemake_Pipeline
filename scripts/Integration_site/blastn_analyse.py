# librairies
import pandas as pd
from argparse import ArgumentParser

##############################################################################################
# ## Notes

""" blastn result sorting and filtering  """

## Load script :

# python3 /blastn_analyse.py -i <inputfile> -o <outputfile> 

##############################################################################################

# Arguments
def config_parameters():
    parser=ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputfile", help=" blastn result ")
    parser.add_argument("-o", "--output", dest="outputfile", help="outputfile of sorted and filtered table")
    args=parser.parse_args()
    if len(sys.argv) < 2 :
        sys.exit("Warning : wrong number of argument")
    return args.inputfile, args.outputfile


inputfile, outputfile = config_parameters()



# Read the CSV file
df = pd.read_csv(inputfile, sep='\t', header=None, on_bad_lines='skip')
df.columns = ["query_id", "subject_id", "coverage", "identity", "evalue", "s_start", "s_end", "strand", "description"]

#sort df by e_value
df.sort_values(by='evalue', ascending=True, inplace=True)
df.reset_index(drop=True, inplace=True)

#threshold for coverage and identity
coverage_threshold = 97
identity_threshold = 90

# Filter the DataFrame based on the thresholds
filtered_df = df[(df['coverage'] >= coverage_threshold) & (df['identity'] >= identity_threshold)]


filtered_df_sort = filtered_df.drop_duplicates(subset='query_id', keep='first')

# Save the sorted and filtered DataFrame to a new CSV file
filtered_df_sort.to_csv(outputfile, index=False)           


