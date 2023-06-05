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
from collections import OrderedDict
import pandas as pd
#############################################################################
# Load the arguments passed to the python script
#############################################################################

parser=argparse.ArgumentParser(
    description='''
***********************************************************************************************************
*** This script reformats the feature annotation table from Ensembl
*** and saves it to a tab-separated file.
*** Usage : python
*** Usage example: python make_annotation_table.py inFile.txt outFile.txt 
***********************************************************************************************************''',
    formatter_class=RawTextHelpFormatter)

parser.add_argument('inFile', type = str)
parser.add_argument('outFile',type = str)

args=parser.parse_args()

#print(args)

inF = ''.join(args.inFile)
outF = ''.join(args.outFile)

print('****** Making annotation table for ',inF,' ******')

#############################################################################
#
#############################################################################
myDict = OrderedDict()
myDict["GeneID"] = [ "Gene_Type", "Gene_Name" , "Description", "assembly", "Chromosome", "start", "end", "strand", "old_locus_tag", "product_accession", "feature_interval_length", "product_length" ]


with open(inF, "r") as file, open(outF, "w+") as out:
    header = file.readline()
    
    for line in file:
        myLine = line.strip().split(sep="\t")
#        print("current line: ", myLine)
        if myLine[0] == "gene":
            start1 = myLine[7]
            end1 = myLine[8]
            geneId1 = myLine[16]
            geneType = myLine[1]
            geneName = myLine[14]
            featLen = myLine[17]
            strand = myLine[9]
            if len(myLine) == 20:
                oldTag = myLine[19].split(sep="=")[-1] 
            elif len(myLine) == 19:
                oldTag = ""
            asem = myLine[2]
            chrom = myLine[6]
            intLen=myLine[15]
            print(geneId1,geneType)
            with open(inF, "r") as ff:
                for l in ff:
                    myLine2 = l.strip().split(sep="\t")
     #               print("next line :", myLine2)
                    start2 = myLine2[7]
                    end2 = myLine2[8]
                    geneId2 = myLine2[16]
    
                    if start1 == start2 and end1 == end2 and geneId1 == geneId2:
                        desc = myLine2[13]
                        if len(myLine2) == 19:
                            prodLen = myLine2[18]
                        elif len(myLine2) == 18:
                            prodLen = ""
                        elif len(myLine2) > 19:
                            prodLen = myLine2[18]
                        prodAcc = myLine2[10]
                        print(geneId1,geneType, geneId2, desc)
                        myDict[geneId1] = [ geneType, geneName, desc, asem,  chrom, start1, end1, strand, oldTag, prodAcc, featLen, prodLen ]
                    else:
                       continue 
        else:
            continue
        
  #  print(myDict)

    df = pd.DataFrame.from_dict(myDict)
   # print(df)
    df = df.transpose()
    df.to_csv(outF, sep="\t", index=True, header = False)

