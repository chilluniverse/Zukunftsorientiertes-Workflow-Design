#!/usr/bin/env python3

import csv
import argparse

#> Argument Parser
def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Add an optional argument
    parser.add_argument('-p', '--path', default="", help='path to folder', type=str)
    parser.add_argument('-a', '--annotation', help='Filename of Annotations', type=str)
    parser.add_argument('-i', '--input', help='Filename of Input bed file', type=str)
    parser.add_argument('-o', '--output', help='Filename of filtered bed file', type=str)

    # Parse the arguments
    args = parser.parse_args()

    return args

#> Read in filtered annotations an write all PeakIDs into an array
def get_PeakID(file_name):
    PeakID = []
    with open(file_name, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        next(csvreader) # skip header row
        for row in csvreader:
            PeakID.append(row[0].split("-")[0])
    
    return set(PeakID)

#> Read in *.narrowPeak file and filter out all lines without valid annotation
def filter_Peaks(file_name, PeakID):
    peaks = []
    with open(file_name, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            if(row[3] in PeakID):
                peaks.append(row)
    return peaks

#> write new *.narrowPeak file
def write_Peaks(file_name, peaks):
    with open(file_name, 'w') as file:
        for subarray in peaks:
            for item in subarray:
                file.write(str(item))
                file.write('\t')
            file.write('\n')
            
            
#- main method        
def main():
    args = parse_arguments()
    
    annotation_file = args.path + args.annotation
    peak_file = args.path + args.input
    filtered_peaks = args.path +args.output
    
    PeakID = get_PeakID(annotation_file)
    peaks = filter_Peaks(peak_file, PeakID)
    write_Peaks(filtered_peaks, peaks)
    print(len(peaks))

if __name__ == '__main__':
    main()