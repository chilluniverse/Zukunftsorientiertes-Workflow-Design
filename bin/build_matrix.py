#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import csv
import sys

parser = argparse.ArgumentParser(description='Compile mapping results to matrix')
parser.add_argument('--files', dest='files', nargs='+', help='Path of *.tsv files containing read count data')
parser.add_argument('--path', dest='path', default= "", help='Output Path')
args = parser.parse_args()


def get_filename_and_extension(file_path):
    file = Path(file_path).name
    filename, file_extension = file.split(os.extsep,1)
    
    return filename, file_extension

#! 1 - DONE
def search_file(files):
    result_dict = {}
    
    for file in files:
        filename, file_extension = get_filename_and_extension(file)
        if file_extension=="sorted.bam.tsv":
            result_dict[filename] = file
    
    return result_dict

def read_tsv_to_dict(file_path):
    result_dict = {}

    with open(file_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        
        for row in tsv_reader:
            result_dict[row[0]] = row[1:4]
            
    return result_dict

#! 2 
def get_file_values(read_count_files):
    result_dict = {}
    for file in read_count_files:
        result_dict[file] = read_tsv_to_dict(read_count_files[file])
    return result_dict
        
def combine_counts(read_count_data):
    result_dict = {}
    for file_name in read_count_data:
        file = read_count_data[file_name]
        for gen in file:
            # gen ist key
            tmp_values = file[gen]
            if result_dict.get(gen) != None:
                tmp_list = result_dict[gen]
            else:
                tmp_list = []
            tmp_list.append(tmp_values[1])
            result_dict[gen] = tmp_list
    
    return result_dict

#! 3
def build_matrix(read_count_data, path):
    
    out_file_name = 'countData'
    header = ""
    
    tmp_length = -1
       
    for file in read_count_data:
        file_length = len(read_count_data[file].values())
        if tmp_length == -1: tmp_length = file_length
        elif tmp_length != file_length: 
            print("files where not mapped to same genome! (mismatch of gene count)")
            sys.exit(1)
            
        # out_file_name = out_file_name + '_' + str(file)
        header = header + "\t" + str(file)
    
    rows = combine_counts(read_count_data)
    
    print ('writing "'+out_file_name+'.tsv"')

    with open(path+out_file_name+".tsv", 'w') as file:        

        file.write(header+'\n')
        
        for key in rows:
            tmp_row = key
            values = rows[key]
            for item in values:
                tmp_row = tmp_row + '\t' + item
            file.write(tmp_row+'\n')
            
    print ("finished ;-)")
    return



read_count_files = search_file(args.files)

read_count_data = get_file_values(read_count_files)

build_matrix(read_count_data,args.path)