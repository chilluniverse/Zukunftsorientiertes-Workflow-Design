#!/usr/bin/env python3

#- Dependencies
import csv
import argparse
import pandas as pd

#- CLI Argument Parser
parser = argparse.ArgumentParser(description='Speperate genes to clusters')
parser.add_argument('--clusters', dest='clusters', help='Path of the clusters file')
parser.add_argument('--path', dest='output_path', default='cluster')
args = parser.parse_args()

# Define a dictionary to store the group values
group_dict = {}

# Read the input tsv file
with open(args.clusters, 'r') as input_file:
    reader = csv.reader(input_file, delimiter=',')

    # Iterate through each row in the tsv file
    for row in reader:
        # Get the first column value (group value) and second column value (id)
        id_value = row[0]
        group_value = row[1]

        if group_value != 'x':
            # If the group value is not already in the dictionary, add it
            if group_value not in group_dict:
                group_dict[group_value] = []

            # Append the id value to the corresponding group value in the dictionary
            group_dict[group_value].append(id_value)


max_length = max([len(item) for item in group_dict.values()])

for gen, value in group_dict.items():
    if len(value) < max_length:
        group_dict[gen] = value + [''] * (max_length - len(value))

clusters = pd.DataFrame(group_dict)

job_file = open("cluster.job", "a")
# Iterate through each column in the dataframe
for column in clusters.columns:
    # Create a new dataframe with the current column as the data and the column name as the header
    new_df = pd.DataFrame(clusters[column].values, columns=[column])
    
    # Write the new dataframe to a new csv file with the column name as the filename
    new_df.to_csv(args.output_path + f'/c{column}.txt', index=False, header=None)
    
    # Write .job file for metascape container
    job_file.write("""{"input":"/data/cluster/c"""+column+""".txt", "output":"/data/results/c"""+column+"""", "single":true, "source_tax_id": 7227, "target_tax_id":7227}\n""")

job_file.close()

