import pandas as pd

file="/home/osidibe/work/PPR_MGEproject/integration_site/nucl_analyse/blastn_out.csv"

# Read the CSV file

df = pd.read_csv(file, sep='\t', header=None, on_bad_lines='skip')

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
filtered_df_sort.to_csv("/home/osidibe/work/PPR_MGEproject/integration_site/nucl_analyse/blastn_out_filtered_sorted_97.csv", index=False)           

# the proportion of the filtered DataFrame compared to the original DataFrame
# proportion = len(filtered_df) / len(df) * 100
# print(f"Proportion of filtered DataFrame: {proportion:.2f}%")



