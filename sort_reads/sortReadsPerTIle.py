
#!/usr/bin/env python3

#############################################################################
### Author: Karolina Bojkowska (karolina.bojkowska@unil.ch)
### Date: 06/01/2023
###
### This script extracts read coordinates for all fastq files from Unaligned projects. 
### It then sorts the read coordinates per tile and writes to file per tile.
### Requires that these files are present per Unaligned_PROJECT_ID in the DEMUX folder:
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
*** This script extracts read coordinates and writes to file per tile for all fastq files from Unaligned projects. ***
*********************************************************************

*** Outputs a folder with a file per tile with read coordiantes associated to Unaligned projects.
*** Usage example: python sortReadsPerTIle.py /work/../DEMUX/221107_JABBA_0000_AHABABA/master_demux.conf 

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
import time
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import fnmatch
from datetime import datetime as dt
import sys


#############################################################################
# Variables
#############################################################################
# dev: define separator
separ='/'#os.sep

#############################################################################
# Functions
#############################################################################
# print information message
def printMsg ( text ):
    """
    Prints a message in a pre-defined format for log purpose. If more than one string is used, 
    concatenate strings with '+' as in Usage example below.
    Usage example: printMsg( text = "This is the message for processing file "+myFile )
    """
    print("[INFO] [",dt.now().strftime('%Y/%m/%d %H:%M:%S').strip("\s"),"] ", text, "\n", sep="")

def rootPath(file):
        # root path
    return(separ.join(file.rsplit(sep=separ)[:-1]))

def extractPaths( file ):
    """
    Process a master_demux.conf config file to extract the Unaligned Suffix and make a list of demuxed projects.
    Returns a tuple with 
    - the list of projects per Unaligned folder as first element and the 
    - absolute path of the RUN_FOLDER as the second element.
    Usage example: 
        myUnaligned_list=extractPaths( file = /work/.../DEMUX/221014_JABBA113_AXXXXX/master_demux.conf )[0]
        RUN_FOLDER_PATH=myUnaligned_list=extractPaths( file = /work/.../DEMUX/221014_JABBA113_AXXXXX/master_demux.conf )[1]
    """
    # check that file exists and is readable
    if not os.access(file, mode=os.R_OK):
        sys.exit("[FATAL] File "+file+" does not exist or is not readable.")
    else:
        pass
    r_path=rootPath(file)    
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
        myList=[]
        for i in myD.keys():
            for a in myD[i]:
                mypath=r_path+separ+'Unaligned'+i
                proj_path=mypath+'/'+a
                myList.append(proj_path)
    print(file, "\n", sep="")
    printMsg("Unaligned project paths :")
    print('\n'.join(myList), '\n')
    return(myList)

def makeFQlist( folder_path_list ):
    """
    Makes a list of FASTQ.GZ files presnt in the folder paths within the input "folder_path_list".
    Usage example: makeFQlist( folder_path_list = ['/home/myPath1', '/home/myPath2'] )
    """
    files_list=list()

    # make a list of files from a list of folder paths
    for path in folder_path_list:
        for i in fnmatch.filter(os.listdir(path), pat='*fastq.gz'):
            files_list.append(path+separ+i)#os.path.join(path,i))
    return(files_list)

def select_files_pattern( file_list , pattern ):
    """
    Selects only pattern - matching files from the provided "file_list".
    Usage example: select_files_pattern( file_list = ['/home/myPath1/myFile_R1_001.fastq.gz', '/home/myPath1/myFile_R2_001.fastq.gz'] ,
    pattern = "*_R1_001*) - will select files matching *_R1_001* only.
    """
    matching_list=fnmatch.filter(file_list , pattern)  
    return(matching_list)

def extractReadCoordinates( header_string ):
    """
    Process a read header string to extract the read coordinates. Returns read coordinates separated by ':'.
    Usage example: extractReadCoordinates( string = 'my_header_string' )
    """
    ss = ':'.join(header_string.split(None)[0].split(sep=':')[3:])
    return(ss)

def extractTile( coord ):
    """
    Process a read header read coordinate string to extract the tile. Returns tile name separated by '_:'.
    Usage example: extractTile( string = 'my_read_coordinates' )
    """
    tt = '_'.join(coord.split(":")[:2])
    return(tt)

def processUnalFqXtractReadCoords( file_list , read_coord_file ):
    """
    Process list of fastq files to extract read coordinates and write to file.
    Usage example: processUnalFqXtractReadCoords( file_list = ['/home/myPath1/myFile_R1_001.fastq.gz', '/home/myPath1/myFile_R2_001.fastq.gz'] ,
    read_coord_file = "/path/to/read_coordinates_file.txt" )
    """
    tiles_set=set()
    for fq in file_list:
        print(fq)
        with gzip.open(fq,"rt") as head_handle, open(read_coord_file, "a+") as headers_file:
            for i in FastqGeneralIterator(head_handle):
                tiles_set.add(extractTile(extractReadCoordinates(i[0])))
                headers_file.write(extractReadCoordinates(i[0])+'\n')
    return(tiles_set, printMsg("Read coordiantes file ready: "+'\n'+read_coord_file))

def makeReadCoorFile( read_coord_file , file_list ):
    """

    """
    if not os.path.exists(read_coord_file):
        printMsg("Extracting read coordinates and writing to file "+read_coord_file)
        return( processUnalFqXtractReadCoords(read_coord_file=read_coord_file, file_list=file_list) )
    else:
        printMsg("File "+read_coord_file+" already present.")
        printMsg("Making tile names set.")
        t_set=set(extractTile(i) for i in open(read_coord_file, "r")) 
        printMsg("Tile names set ready.")
        return(t_set)
    
        
def sortReadCoorsPerTile( conf_file ):

    # set variables
    ROOT_PATH=rootPath(conf_file)
    tmp_directory_name='perTileReadsTMP'
    
    unal_paths=extractPaths(conf_file)

    # make list of R1 fastq files to be used for extracting read coordiantes
    read1_file_list=select_files_pattern(file_list = makeFQlist(unal_paths), pattern ="*_R1_001.*")

    lanes=( 'L001', 'L002', 'L003', 'L004' )

    for lane in lanes: 
        printMsg(".............. Processing lane "+lane+" ..............")

        tmp_dir=ROOT_PATH+separ+lane+tmp_directory_name
            
        #make temporary folder for read_coordinates per tile
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
            printMsg("Temporary directory for read coordinates per tile created.")
        else:
            printMsg("Temporary directory for read coordinates per tile already exists.")
        
        # filter file list for lane match
        L_read1_file_list=fnmatch.filter(read1_file_list, '*_'+lane+'_*' )
        
        # define file name 
        read_coord_file=tmp_dir+separ+lane+'readsCoordiates.txt'

        if L_read1_file_list:
        # make read coordiantes file and specify tiles names
            tiles=makeReadCoorFile( file_list=L_read1_file_list, read_coord_file = read_coord_file)
        # crete tile files
            for i in tiles:
                if not os.path.exists(tmp_dir+separ+i):
                    printMsg("Initializing tile files for read coordinates sorting.")
                    with open(tmp_dir+separ+i, "w"):
                        pass
                else:
                    printMsg("Read coordinates per tile already sorted.")
                    return None
            
            #print(files)
            with open(read_coord_file, "r") as rcf:
                for lane in rcf:
                    #print(lane)
                    tt=extractTile(lane.strip())
                    file=tmp_dir+separ+tt
                    with open(file, "a") as f:
                        f.write(lane)
            printMsg("Sorting read coordinates per tile for lane "+lane+" completed.")
        else:
            printMsg("Lane "+lane+" not present.")
            pass



#############################################################################
#                                      MAIN
#############################################################################

printMsg("............................ Starting python module ............................")

sortReadCoorsPerTile(conf_file = my_conf)

printMsg("............................ End of python module ............................")

#############################################################################
#                                    End
#############################################################################


