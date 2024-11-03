#!/usr/bin/env python3

import os
import argparse

# Argument Parser
def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Add an optional argument
    parser.add_argument('-p', '--path', default='',  help='path to folder', type=str)
    parser.add_argument('-g', '--genes', help='Name or Path to gene list', default="peaks/annotation_genes_only.txt", type=str)
    parser.add_argument('-a', '--annotation', help='Filename of Annotations', default="annotation_filtered.csv", type=str)
    parser.add_argument('-c', '--count', help='threshold at which genes should be kept', default=10, type=int)
    parser.add_argument('-o', '--outname', help='filename of output')

    # Parse the arguments
    args = parser.parse_args()

    return args

#> Import all gene names
#? read in all genes given in text file and return it as a dict; gene names := key
def get_gene_names(filepath):
    dict = {}                                       # create dictionary dict

    with open(filepath, 'r') as file:               # Open the file
        lines = file.readlines()                    # and save each line in list

    for line in lines:
        line = line.strip()                         # strip whitespaces/new line chars
        if line and line.startswith("FBgn"):        # If the line /= empty and is a gene, add it to the dictionary (init with 0 as count) and list
            dict[line] = 0

    return dict

#> Count all apperances of the genes in all annoations
def count_genes(folder, annotation, gene_count):
    # extract the value of the 12th column of a tsv file?
    # Open the annotation-file and read its content
    with open(annotation, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                values = line.split('\t')
                if len(values) >= 12:
                    gene = values[11]
                    if gene in gene_count:
                        gene_count[gene] += 1
    return gene_count

def filter_apperance(count, gene_count):
    # Discard key and value from dictionary if value is smaller then value of variable "count"?
    dict = {key: value for key, value in gene_count.items() if value >= count}
    return dict

def write_genes(path, gene_count, outname):
    # Write keys of dictionary to file?
    file = open(path + outname, "w")
    for key in gene_count.keys():
        file.write(key +"\n") #+ "\t" + str(gene_count[key]) + "\n")
    file.close()

#! Main method
def main():
    args = parse_arguments()
    gene_list_file = args.path + args.genes
    gene_counts = get_gene_names(gene_list_file)
    gene_counts = count_genes(args.path, args.annotation, gene_counts)
    gene_counts = filter_apperance(args.count, gene_counts)

    write_genes(args.path, gene_counts, args.outname)

    print(len(gene_counts))

if __name__ == '__main__':
    main()