#!/usr/bin/env python3

#- Dependencies
import argparse
from pathlib import Path
import os
import csv

#- CLI Argument Parser
parser = argparse.ArgumentParser(description='Compile mapping results to matrix')
parser.add_argument('--files', dest='files', help='Path of *.tsv files containing read count data')
parser.add_argument('--name', dest='name', default= 'countData', help='Name of the output file; default: countData')
parser.add_argument('--path', dest='path', default= "", help='Output Path')
args = parser.parse_args()

#? Convert Nextflow Input ("[/path/to/file1, /path/to/file2]")to Array
file_paths_str = args.files.strip('[]')
file_paths = [path.strip() for path in file_paths_str.split(',')]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!            FUNCTIONS            !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#* split filename and extension
#? helper function used in function 'search_file()'
def get_filename_and_extension(file_path):
    file = Path(file_path).name
    filename, file_extension = file.split(os.extsep,1)
    
    return filename, file_extension

#- 1. Check Array Input and convert to dictionary
# @inputs:       files = [/path/to/file1, /path/to/file2]
# @output: result_dict = {'filename':'/path/to/file1'}
def search_file(files):
    result_dict = {}
    
    for file in files:
        filename, file_extension = get_filename_and_extension(file)
        if file_extension=="sorted.bam.tsv":
            result_dict[filename] = file
    return result_dict

#* read tsv file and write to a dictionary
#? helper function used in function 'get_file_values()'
def read_tsv_to_dict(file_path):
    result_dict = {}

    with open(file_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        for row in tsv_reader:
            result_dict[row[0]] = row[2:]
        
    return result_dict

#- 2. extract values from files and write to dictionary
# @inputs: read_count_files = {'filename':'/path/to/file1'}
# @output: result_dict
def get_file_values(read_count_files):
    result_dict = {}
    for file in read_count_files:
        result_dict[file] = read_tsv_to_dict(read_count_files[file])

    return result_dict

#* put data in a single data structure
def combine_counts(read_count_data):
    
    #? get all genes / gene-names
    all_genes = set([])
    print(read_count_data)
    for runs in read_count_data:
        run = read_count_data[runs]
        
        for gen in run:
            all_genes.add(gen)
    
    
    result_dict = {}
    for gen in all_genes:
        if result_dict.get(gen) != None:   
                tmp_list = result_dict[gen]
        else: tmp_list = []
        
        for run in read_count_data:
            run_data = read_count_data[run]
            
            if gen in run_data: 
                tmp_list.append(run_data[gen][0])
            else: tmp_list.append('0')
        result_dict[gen] = tmp_list

    return result_dict

#- 3. Build the Matrix and generate output file
# @inputs: read_count_dat: = {file:data}
# @                   path = 'path/to/folder' 
# @          out_file_name = 'filename'
# @output: 
def build_matrix(read_count_data, path, out_file_name):
    #? write header "\tRUN\t\RUN..."
    header = ""
    for run in read_count_data: header = header + "\t" + str(run)

    #? combine all rows from all files
    rows = combine_counts(read_count_data)                              

    with open(path+out_file_name+".tsv", 'w') as file:

        file.write(header+'\n')
        
        for run in rows:
            tmp_row = run
            values = rows[run]
            for item in values:
                tmp_row = tmp_row + '\t' + item
            file.write(tmp_row+'\n')

    return


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!         FUNCTION CALLS          !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

read_count_files = search_file(file_paths)

read_count_data = get_file_values(read_count_files)

build_matrix(read_count_data,args.path, args.name)