#!/bin/python3
#!/usr/bin/env python3

#############################################################################
### Author: Karolina Bojkowska (karolina.bojkowska@unil.ch)
###
### Uses two arguments:
### path = path to the /../Unaligned/GTF/FastqC folder
### metric = 'metric name from the FastqC output' as one string
#############################################################################

#############################################################################
# Import modules I
#############################################################################

import argparse
from argparse import RawTextHelpFormatter
import pandas
import os
import glob

#############################################################################
# Load the arguments passed to the python script
#############################################################################

parser=argparse.ArgumentParser(
    description='''
***********************************************************************************************************
*** This script retrieves fastqc metrics for 
*** and saves it to a tab-separated file.
*** Usage : python 
*** Usage example: python get_sra_metadata.py /path/to/GTFQC/FastqC '>>Per seuence quality' '/path/my_file.txt'
***********************************************************************************************************''',
    formatter_class=RawTextHelpFormatter)

parser.add_argument('path', type = str)
parser.add_argument('metric',type = str)
parser.add_argument('out_file',type = str)

args=parser.parse_args()

#print(args)

path = ''.join(args.path)
metric = ''.join(args.metric)
out_file = ''.join(args.out_file)

print('****** Processing for fastqc folder ', path,' and metric ', metric, ' ******')

#############################################################################
# 
#############################################################################


dir_pattern = '*_fastqc'

fqc_output_file_name = 'fastqc_data.txt'

final_path=os.path.join(path, dir_pattern, fqc_output_file_name)

print(final_path)

def get_sample_id( myPath ):
    """
    Gets the fastq id from the path provided.
    """

    sampleID = '_'.join(myPath.split(sep = "/")[-2].split(sep = '_')[:-1])

    return sampleID

fqc_files_list = glob.glob(final_path)

#print(fqc_files_list)


end_module_string = '>>END_MODULE'


with open(out_file, 'w+') as outFile:
    for f in fqc_files_list:
        sample_id = get_sample_id( f )
#        print(sample_id)
        with open (f , 'r') as inFile:
            for line in inFile:
                if metric in line and line.startswith('>>'):
                    for line in inFile:
                        if line.startswith('#'):
                            new_header = [line.strip(), 'Sample ID', '\n']
                            if not '\t'.join(new_header) in outFile.read():
                                outFile.write(('\t'.join(new_header)))
                            else:
                                continue
                        else:
                            if line.strip() == end_module_string:
                                break
                            new_line=[line.strip(), sample_id, '\n']
                            outFile.write(('\t'.join(new_line)))

print('****** Processing completed ******')

