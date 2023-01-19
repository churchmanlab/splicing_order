"""

Author: Karine Choquet

Date: March 3, 2020

This script will determine the splicing status of individual introns and reads from nanopore sequencing data

Usage: python get_multi_introns_isoforms_df.directRNA.py splicing_dict multi_introns_df multi_introns_counts

splicing_dict: dictionary (key: read) of dictionaries (key: transcript, value: list of introns and splicing statuses)


"""

import sys
import numpy as np
import pandas as pd
import pysam
from collections import Counter

import re
import math

import pybedtools
from pybedtools import BedTool

#####################################

# CONFIG


# BED file containing reads that intersect introns and the info from both
intersect_df = sys.argv[1]

# Output files
out_splicing_info = sys.argv[2]
out_splicing_dict = sys.argv[3]

######################################



# function to create a dataframe with splicing information for
# every read that spans an intron in the dataset
def get_splicing_info(intersect_df, min_overlap):
    
    df = open(intersect_df, 'r')

    # prepare a list for splice calls
    spliceCalls = []

    # set variables for parsing the cigar string
    pattern = re.compile('([MIDNSHPX=])')
    Consumes_Query = ["M", "I", "S", "=", "X"]
    Consumes_Reference = ["M", "D", "N", "=", "X"]    

    # loop through all read-intron intersects
    for line in df:

        split_line = line.split("\t")
        df_count = int(split_line[-1])

        # ignore reads that do not overlap intron by minimum threshold
        if (df_count < min_overlap):
            continue

        # record the start and ends of reads
        # record the start and ends of reads
        aln_start = int(split_line[1]) # record the start of the read
        aln_end = int(split_line[2]) # record the end of the read
        aln_name = split_line[3]
        intron_chr = split_line[7]
        intron_start = int(split_line[8]) # record the end of the intron
        intron_end = int(split_line[9]) # record the end of the intron
        aln_strand = split_line[5] # strand of the read
        gene_strand = split_line[-2]
        cigar_aln = split_line[6] # cigar string
        name_gene = split_line[-4]
        intron_count = split_line[-3]

        # parse cigar string into a list of tuples for easy parsing
        Sep_Values = Sep_Values = pattern.split(cigar_aln)[:-1]
        CigarPairs = list((Sep_Values[n:n+2] for n in range(0, len(Sep_Values), 2)))  

        # set up variables for measuring the length of cigar string operators
        CigarOp_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        start_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        end_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        intron_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        currentloc = aln_start

        # go through list of cigar strings and grab splicing information
        for cigar_Entry in CigarPairs:

            op_Length = int(cigar_Entry[0]) # get length of cigar operator
            cigarOp = cigar_Entry[1] # get type of cigar operator  
            CigarOp_counts[cigarOp] += op_Length # add the cigar operator length to the counts dictionary
            cigarOp_start=currentloc # get the starting coordinate of the cigar operator

            if (cigarOp in Consumes_Reference):
                currentloc=currentloc+op_Length # add the cigar operator length to the current location coordinate 

            cigarOp_end=currentloc # get the ending coordinate of the cigar operator

            # gather information if the portion of the cigar string spans the designated intron start
            if (cigarOp_start<intron_start-min_overlap and cigarOp_end>=intron_start-min_overlap):
                if (cigarOp_end>=intron_start+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(intron_start-min_overlap)+1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=intron_start-min_overlap and cigarOp_end<intron_start+min_overlap):
                count=op_Length
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary       

            elif (cigarOp_start<intron_start+min_overlap and cigarOp_end>=intron_start+min_overlap):
                if (cigarOp_start<=intron_start-min_overlap):
                    count=min_overlap*2
                else:
                    count=(intron_start+min_overlap)-cigarOp_start-1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string is within the intron
            if (cigarOp_start<intron_start and cigarOp_end>=intron_start):
                if (cigarOp_end>=intron_end):
                    count=intron_end-intron_start
                else:
                    count=cigarOp_end-intron_start
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=intron_start and cigarOp_end<intron_end):
                count=op_Length
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<intron_end and cigarOp_end>=intron_end):
                if (cigarOp_start<=intron_start):
                    count=intron_end-intron_start
                else:
                    count=intron_end-cigarOp_start
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string spans the designated intron end
            if (cigarOp_start<intron_end-min_overlap and cigarOp_end>=intron_end-min_overlap):
                if (cigarOp_end>=intron_end+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(intron_end-min_overlap)
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=intron_end-min_overlap and cigarOp_end<intron_end+min_overlap):
                count=op_Length
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<intron_end+min_overlap and cigarOp_end>=intron_end+min_overlap):
                if (cigarOp_start<=intron_end-min_overlap):
                    count=min_overlap*2
                else:
                    count=(intron_end+min_overlap)-cigarOp_start
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

        # get length of the aligned portion of this read from cigar string
        aligned_read_length = CigarOp_counts['M']+CigarOp_counts['D']

       # get 5'SS and 3'SS counts as determined by gene strand
        strand = gene_strand
        if (strand == '+'):
            aln_start = aln_start  # record the start of the read
            aln_end = aln_end # record the end of the read
            intron_5SS_counts = start_counts # record the cigar string counts over the 5'SS
            intron_3SS_counts = end_counts # record the cigar string counts over the 3'SS
            read_overlap = intron_end - (aln_end - aligned_read_length + min_overlap)
            
        if (strand == '-'):
            aln_start = aln_end # record the start of the read
            aln_end = aln_start  # record the end of the read
            intron_5SS_counts = end_counts # record the cigar string counts over the 5'SS
            intron_3SS_counts = start_counts # record the cigar string counts over the 3'SS  
            read_overlap = (aln_end + aligned_read_length - min_overlap) - intron_start
            
        # annotate splicing status based on CIGAR string information around splice sites
        splice='UND'

        if (intron_5SS_counts['N']==0 and intron_3SS_counts['N']==0):
            if (intron_3SS_counts['M']+intron_3SS_counts['D']==min_overlap*2):
                if (intron_3SS_counts['M']>min_overlap):
                    splice = 'NO'

        if (intron_5SS_counts['N']>0 and intron_5SS_counts['N']<min_overlap*2):
            if (intron_3SS_counts['N']>0 and intron_3SS_counts['N']<min_overlap*2):
                splice = 'YES'
            if (intron_3SS_counts['N']==min_overlap*2):
                splice = 'SKP3SS'

        if (intron_5SS_counts['M']==0 and intron_3SS_counts['M']==0):
            if (intron_3SS_counts['N']==min_overlap*2):
                    splice = 'SKP'

        if (intron_5SS_counts['N']==min_overlap*2):
            if (intron_3SS_counts['N']>0 and intron_3SS_counts['N']<min_overlap*2):
                splice = 'SKP5SS'


        # annotate splicing status based on CIGAR string information within the intron 
        if (splice == 'YES'):
            if (float(intron_end-intron_start) > 0.0):
                ratio = float(intron_counts['N'])/float(intron_end-intron_start)
                difference = abs(intron_counts['N']-(intron_end-intron_start))

                # if read is spliced, between 90-100% of the intron has to be spliced 
                # and no more than 100 nucleotides within the intron can be matching the intron sequence
                if( ratio < 0.8 or ratio > 1.1 or difference > 50):
                    splice='UND'
            if (float(intron_end-intron_start) == 0.0):
                splice='UND'

        if (splice == 'NO'):
            if (float(intron_end-intron_start) > 0.0):
                ratio = float(intron_counts['M'])/(float(intron_counts['M'])+float(intron_counts['N'])+float(intron_counts['D'])+1)

                # if read is unspliced, at least 75% of the read has to match (CIGAR=M) the intron sequence
                if(intron_counts['M'] < min_overlap/2 or ratio < 0.75):
                    splice='UND'
            
            if (float(intron_end-intron_start) == 0.0):
                splice='UND'

        if (splice == 'SKP_3SS') or (splice == 'SKP_5SS'):
            if (float(intron_end-intron_start) > 0.0):
                ratio = float(intron_counts['N'])/float(intron_end-intron_start)
                difference = abs(intron_counts['N']-(intron_end-intron_start))

                # if read is spliced, between 90-110% of the intron has to be spliced
                # and no more than 100 nucleotides within the intron can be matching the intron sequence
                if( ratio < 0.8 or ratio > 1.1 or difference > 50): # changed this for DDX39A
                    splice='UND'
            if (float(intron_end-intron_start) == 0.0):
                splice='UND'


        # save read, intron, and splicing information
        spliceCalls.append([aln_name,intron_chr,intron_start,intron_end,gene_strand,name_gene,intron_count,read_overlap,splice])

    spliceCalls_df = pd.DataFrame(spliceCalls)
    spliceCalls_df.columns = ["read_name","chrom","intron_start","intron_end","strand","gene_name","intron_count","read_overlap","splice_status"]

    return spliceCalls_df


# every read that spans an intron in the dataset
def get_read_junctions_dictionary(splice_df):

    read_junctions = {}

    for i in range(0,splice_df.shape[0]):       

        # define the read name
        read_name = splice_df['read_name'].iloc[i]
        gene_name = splice_df['gene_name'].iloc[i]
        chrom = splice_df['chrom'].iloc[i]
        intron_start = splice_df['intron_start'].iloc[i]
        intron_end = splice_df['intron_end'].iloc[i]
        intron_count = int(splice_df['intron_count'].iloc[i])
        strand = splice_df['strand'].iloc[i]
        read_overlap = splice_df['read_overlap'].iloc[i]
        splice_status = splice_df['splice_status'].iloc[i]

        # check if read name is in the dictionary, if not save it
        if read_name not in read_junctions.keys():

            # make a new dictionary for the gene and add intron info to it
            read_junctions[read_name] = {}
            read_junctions[read_name][gene_name] = [[chrom, intron_start, intron_end, intron_count, strand, read_overlap, splice_status]]

        # check if read name is in the dictionary, if it is proceed to gene information
        elif read_name in read_junctions.keys():

            # if gene_name is not already in read dictionary, 
            # make a new dictionary for the gene and add intron info to it
            if gene_name not in read_junctions[read_name].keys():
                read_junctions[read_name][gene_name] = [[chrom, intron_start, intron_end, intron_count, strand, read_overlap, splice_status]]

            # if gene_name is already in read dictionary, add new intron info to it
            elif gene_name in read_junctions[read_name].keys():
                read_junctions[read_name][gene_name].append([chrom, intron_start, intron_end, intron_count, strand, read_overlap, splice_status])

    return read_junctions



##############################################

# set all variables for analysis
min_overlap = 25


# get splicing information for every read that spans an intron
print("getting splicing info...")
splice_info = get_splicing_info(intersect_df,min_overlap)
splice_info = splice_info.sort_values(by=['chrom','intron_start','intron_end']).reset_index(drop=True)
splice_info.to_csv(out_splicing_info, sep='\t', index=False, header=True)

splice_dictionary = get_read_junctions_dictionary(splice_info)
out_sub_dict = out_splicing_dict + '.npy'
np.save(out_sub_dict, splice_dictionary)

