
#!/usr/bin/env python3

#############################################################################
### Author: Karolina Bojkowska (karolina.bojkowska@unil.ch)
### Date: 07/11/2022
###
### This script cleans the Undetermined.fastq.gz removing redundant reads that are present in all project Unaligned fastq.gz within the DEMUX folder.
###
### Requires that these files are present per Unaligned_PROJECT_ID in the DEMUX folder:
### Undetermined_R1.fastq.gz per project
### Unaligned_R1.fastq.gz
###
### Uses as argument the master_demux.config file present in the DEMUX folder.
#############################################################################

#############################################################################
# Import modules I
#############################################################################

import argparse
from argparse import RawTextHelpFormatter

#############################################################################
# Load the arguments passed to the python script
#############################################################################

parser=argparse.ArgumentParser(
    description='''*********************************************************************
*** This script sorts undetermined for bcl2fastq pipeline output. ***
*********************************************************************

*** Removes all sample-attributed reads from Undetermined.fastq.gz reads per DEMUX declared in the config_file.
*** Processes R1, R2, I1, I2 and R3 separately.
*** Outputs a new Undetermined_clean.fastq.gz file per DEMUX contaning only non-attributed reads (real Undetermined).
*** Usage example: python clean_undetermined.py /work/../DEMUX/221107_JABBA_0000_AHABABA/master_demux.conf 

*** Before using on DCSR, activate the appropriate conda environment:
    module load gcc miniconda3
    conda activate demux_dev''',
epilog="""*** All is well that ends well.""", formatter_class=RawTextHelpFormatter )
parser.add_argument('config_file',
help='''Config file master_demux.conf that must be present in the DEMUX_RUN directory.
It defines all demultiplexings done for the Sequencing Run. Example master_demux.conf :

[DEMUX]
unaligned_suffix=BCL1
OverrideCycles=Y151;I10N6;I10;Y151
BarcodeMismatchesIndex1=1
sample_project=SpiderConeSnail_MR,CoPo_JM,PhixSeqC_LG

[DEMUX]
unaligned_suffix=BCL2
OverrideCycles=Y151;I8Y8;I8N2;Y151
BarcodeMismatchesIndex1=1
sample_project=TimemaChip_TS
''')
args=parser.parse_args()

my_conf = ''.join(args.config_file)

#############################################################################
# Import modules II
#############################################################################

import gzip
import sys
import time
import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import fnmatch
from datetime import datetime as dt
#import random

#############################################################################
# Functions
#############################################################################
def extractUnalignedPaths( file ):
    """
    Process a master_demux.conf config file to extract the Unaligned Suffix and make a list of demuxed projects.
    Returns a list of projects per Unaligned folder.
    Usage example: processConfig( file = /work/.../DEMUX/221014_JABBA113_AXXXXX/master_demux.conf )
    """
    # root path
    r_path='/'.join(file.rsplit(sep='/')[:-1])
    # read config file
    with open(file, 'r') as f:
        printMsg("Processing the following config file :")
        myD={}
        for lane in f:
            if not lane.isspace():
                if lane.strip().split(sep="=")[0] == "unaligned_suffix":
                    myKey=lane.strip().split(sep="=")[1]
                else:
                    if lane.strip().split(sep="=")[0] == "sample_project":
                        myItem=lane.strip().split(sep="=")[1]
                        myD[myKey]=[i for i in myItem.split(sep=",")]
#        print(myD)
        myList=[]
        for i in myD.keys():
            for a in myD[i]:
                mypath=r_path+'/Unaligned'+i
#                projectID='_'.join(a.split(sep='_')[-2:])
                proj_path=mypath+'/'+a
                myList.append(proj_path)
    print(file, "\n", sep="")
    printMsg("Unaligned project paths :")
    print('\n'.join(myList), '\n')
    return(myList)

def extractReadCoordinates( header_string ):
    """
    Process a read header string to extract the read coordinates. Returns read coordinates separated by ':'.
    Usage example: extractReadCoordinates( string = 'my_header_string' )
    """
    ss = ':'.join(header_string.split(None)[0].split(sep=':')[3:])
    return(ss)

def makeFileList( folder_path_list ):
    """
    Makes a list of files presnt in the folder paths within the input "list of folder_paths".
    Usage example: clean_fastq( folder_path_list = ['/home/myPath1', '/home/myPath2'] )
    """
    files_list=list()
    # make a list of files from a list of folder paths
    for path in folder_path_list:
        for i in fnmatch.filter(os.listdir(path), pat='*gz'):
            files_list.append(os.path.join(path,i))
    return(files_list)

def printMsg ( text ):
    print("[INFO] [",dt.now().strftime('%Y/%m/%d %H:%M:%S').strip("\s"),"] ", text, "\n", sep="")

def makeSet2Eliminate( folder_path_list, lane_pattern, output_file ):
    """
    Makes a set with fastq headers from the fastq files present in the list folder_paths that match the wanted read pattern (R1, R2, R3 or I1, I2).
    Usage example: clean_fastq( folder_path_list = ['/home/myUnalignedPath1', '/home/myUnalignedPath2'] , lane_pattern = 'L001' , read_pattern = '*_R1_*' , output_file = "read_coordinates2remove.txt" )
    """
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        printMsg("Headers file \n\n",output_file,"\n\nalready exists and is non-empty. Loading data into unwanted headers set.")
        unwanted_set=set()
        with open(output_file, 'r') as hf:
            for i in hf:
                unwanted_set.add(i.strip())
        #        print(unwanted_set)
        printMsg("Done importing headers set from file.")
    else:
        with open(output_file, "w") as headers_file:
            pass
        printMsg("Empty headers file ",output_file," created.")
        
        printMsg("Processing the following project fastq files from lane ", lane_pattern, " to extract headers of reads to be eliminated from the Undetermined fastq files.")
        process_unal_fq_files_xtract_headers(output_file=output_file, L_read1_file_list=folder_path_list)
        printMsg("Done making headers file for ", lane_pattern,".")
        printMsg("Importing headers file for ", lane_pattern,".")
        with open(output_file, "r") as headers_file:
            unwanted_set=set()
            for ll in headers_file:
                unwanted_set.add(ll.strip())
        printMsg("Done making Headers set of reads to be eliminated for ", lane_pattern, ".")
    return(unwanted_set)

def process_unal_fq_files_xtract_headers(output_file, L_read1_file_list):
    for fq in L_read1_file_list:
        print(fq)
        with gzip.open(fq,"rt") as head_handle, open(output_file, "a") as headers_file:
            for i in FastqGeneralIterator(head_handle):
                r_cor=extractReadCoordinates(i[0])
 #               line_found = any(r_cor in line for line in headers_file)
  #              if not line_found:
   #                 headers_file.seek(0, os.SEEK_END)
                headers_file.write(r_cor+'\n')

def cleanUndetermined( headers_set, undetermined_fq, out_file ):
    """
    Uses the set with fastq headers provided as ARG1 to clean the fastq file goven as ARG2 to output file in ARG3.
    Writes the final reads to the out_file.
    Usage example: cleanUndetermined( headers_set = {'', '', ''}, undetermined_fq = '/home/Undaligned_PROJECT_1/Undetermined/Undetermined_L001_R1_001.fastq.gz',
    out_file = '/home/Undaligned_PROJECT_1/Undetermined/Undetermined_clean.fastq.gz' )
    """
#    final_set=set()
    with gzip.open(undetermined_fq,"rt") as undet_handle, gzip.open(out_file,"wt") as out_handle:
        # go over Undetermined fastq headers read positions; pick those that do not match the headers_set and write to out_file
        for j in FastqGeneralIterator(undet_handle):
            z=extractReadCoordinates(j[0])
            #print(j, '\n',z)
            # test if read position present
            if not z in headers_set:
                #print(z)
                out_handle.write("@%s\n%s\n+\n%s\n" % (j))
    #return(print("[INFO] [",dt.now().strftime('%Y/%m/%d %H:%M:%S').strip("\s"),"] CleanUndetermined completed.","\n", sep=""))

def clean_fastq( config_file ):
    """
    Removes reads from the target Undetermiend FQ.GZ per Demux declared in the config_file for all Undetermined reads (R1, R2, I1, I2, R3) separately.
    Uses the fastq.gz files per Unaligned Project. Outputs a new, cleaned out_FQ.GZ file per .
    Usage example: clean_fastq( config_file = master_demux.conf )
    """
    st=time.time()
    my_project_list=extractUnalignedPaths(config_file)
    #print(my_project_list)
    file_list=makeFileList(my_project_list)
    read_pattern='*_R1_00?.fastq.gz'

    read1_file_list=fnmatch.filter(file_list, read_pattern)
    
    lanes=( 'L001', 'L002', 'L003', 'L004' )

    for L in lanes:
        printMsg(".............. Processing lane ",L," ..............")
        L_read1_file_list=fnmatch.filter(read1_file_list, '*_'+L+'_*' )

        if L_read1_file_list:
            headers_file='/'.join(config_file.split(sep="/")[:-1])+'/'+L+'read_coordinates_to_eliminate.txt'
            headers2remove = makeSet2Eliminate(folder_path_list = L_read1_file_list, lane_pattern=L, output_file = headers_file)
    
            fq_list=set()
            for path in my_project_list:
                undet_FQ_path='/'.join(path.split(sep="/")[:-1])+'/Undetermined'
            #    print(undet_FQ_path)
                myfqs=(fnmatch.filter(os.listdir(undet_FQ_path), '*_'+L+'_*gz')) 
                for i in myfqs:
                    undet_FQ=undet_FQ_path+'/'+i
                    fq_list.add(undet_FQ) 
            printMsg("Undetermined fastqs :\n\n", '\n'.join(fq_list))     
            make_clean_undetermined(headers2remove, fq_list)
        else:
            printMsg("Lane ", L, " not present.")
            pass
   
#    with open("headers_set.txt", "w") as of:
 #       of.write("\n".join(headers2remove))
    et = time.time()
    # get the execution time
    elapsed_time = et - st
    return(printMsg("Cleaning Undetermined completed. Execution time: ", round(elapsed_time/60,1), ' minutes.'))

def make_clean_undetermined(headers2remove, fq_list):
    for undet_FQ in fq_list: 
        if os.path.isfile(undet_FQ):
            printMsg("Starting fastq cleaning module for \n\n", undet_FQ)
            new_Undet_path='/'.join(undet_FQ.split("/")[:-1])+'_clean' #os.path.join(undet_FQ_path,'Undetermined')
            if not os.path.exists(new_Undet_path):
                os.mkdir(new_Undet_path,mode=0o777)
                printMsg("Directory ", new_Undet_path, " created.")
            else:
                printMsg("Directory ", new_Undet_path, " already exists.")

            file_name=undet_FQ.split(sep="/")[-1].replace('Undetermined_', 'Undetermined_clean_')
            out_file_FQ=os.path.join(new_Undet_path, file_name)
            # clean the Undetermined fastq
            if os.path.isfile(out_file_FQ): 
                printMsg("Output file ",out_file_FQ," already present.")      
                pass
            else:
                cleanUndetermined(headers2remove,undet_FQ,out_file_FQ)
                printMsg("Output file created: ",out_file_FQ)
            printMsg("Fastq cleaning module for ", undet_FQ," completed.")
        else:
            print('File does not exist; skipping ',undet_FQ)
            pass


#############################################################################
#                                      MAIN
#############################################################################

printMsg("............................ Starting python module ............................")
clean_fastq(config_file = my_conf)
printMsg("............................ End of python module ............................")

#############################################################################
#                                    End
#############################################################################


