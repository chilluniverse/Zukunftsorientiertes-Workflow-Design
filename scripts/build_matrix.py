
import argparse
from pathlib import Path
import os
import csv
import sys

parser = argparse.ArgumentParser(description='Compile mapping results to matrix')
parser.add_argument('--path', dest='path', help='Path of folder and/or subfolder containing SRA Runs')
args = parser.parse_args()


def get_filename_and_extension(file_path):
    file = Path(file_path).name
    filename, file_extension = file.split(os.extsep,1)
    
    return filename, file_extension

def search_file(folder_path):
    result_dict = {}
    
    for root, files in os.walk(folder_path):
        for file in files:
            filename, file_extension = get_filename_and_extension(file)
            if file_extension=="sorted.bam.txt":
                result_dict[filename] = os.path.join(root, file)
    
    return result_dict

def read_tsv_to_dict(file_path):
    result_dict = {}

    with open(file_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        
        for row in tsv_reader:
            result_dict[row[0]] = row[1:4]
            
    return result_dict

def get_run_values(run_paths):
    result_dict = {}
    for run in run_paths:
        result_dict[run] = read_tsv_to_dict(run_paths[run])
    return result_dict
        
def combine_counts(run_values):
    result_dict = {}
    for run_name in run_values:
        run = run_values[run_name]
        for gen in run:
            # gen ist key
            tmp_values = run[gen]
            if result_dict.get(gen) != None:
                tmp_list = result_dict[gen]
            else:
                tmp_list = []
            tmp_list.append(tmp_values[1])
            result_dict[gen] = tmp_list
    
    return result_dict

def build_matrix(run_values, path):
    
    out_file_name = 'countData'
    header = ""
    
    tmp_length = -1
       
    for run in run_values:
        run_length = len(run_values[run].values())
        if tmp_length == -1: tmp_length = run_length
        elif tmp_length != run_length: 
            print("runs where not mapped to same genome! (mismatch of gene count)")
            sys.exit(1)
            
        # out_file_name = out_file_name + '_' + str(run)
        header = header + "\t" + str(run)
    
    rows = combine_counts(run_values)
    
    print ('writing "'+out_file_name+'.tsv"')

    with open(path+"/"+out_file_name+".tsv", 'w') as file:        

        file.write(header+'\n')
        
        for key in rows:
            tmp_row = key
            values = rows[key]
            for item in values:
                tmp_row = tmp_row + '\t' + item
            file.write(tmp_row+'\n')
            
    print ("finished ;-)")
    return



run_paths = search_file(args.path)
run_values = get_run_values(run_paths)

build_matrix(run_values,args.path)