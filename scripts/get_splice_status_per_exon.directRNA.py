"""

Author: Karine Choquet

Date: March 3, 2020

This script will determine the splicing status of individual exons and reads from nanopore sequencing data

Usage: python get_multi_exons_isoforms_df.directRNA.py splicing_dict multi_exons_df multi_exons_counts

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


# BED file containing reads that intersect exons and the info from both
intersect_df = sys.argv[1]

# Output files
out_splicing_info = sys.argv[2]
out_splicing_dict = sys.argv[3]

######################################



# function to create a dataframe with splicing information for
# every read that spans an exon in the dataset
def get_splicing_info_exons(intersect_df, min_overlap):
    
    df = open(intersect_df, 'r')

    # prepare a list for splice calls
    spliceCalls = []

    # set variables for parsing the cigar string
    pattern = re.compile('([MIDNSHPX=])')
    Consumes_Query = ["M", "I", "S", "=", "X"]
    Consumes_Reference = ["M", "D", "N", "=", "X"]    

    # loop through all read-exon intersects
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
        exon_chr = split_line[7]
        exon_start = int(split_line[8]) # record the end of the intron
        exon_end = int(split_line[9]) # record the end of the intron
        aln_strand = split_line[5] # strand of the read
        gene_strand = split_line[-2]
        cigar_aln = split_line[6] # cigar string
        name_gene = "NM_" + split_line[-4].split("_")[1]
        exon_count = int(split_line[-4].split("_")[3])

        # parse cigar string into a list of tuples for easy parsing
        Sep_Values = pattern.split(cigar_aln)[:-1]
        CigarPairs = list((Sep_Values[n:n+2] for n in range(0, len(Sep_Values), 2)))  

        # set up variables for measuring the length of cigar string operators
        CigarOp_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        start_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        end_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        exon_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
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

            # gather information if the portion of the cigar string spans the designated exon start
            if (cigarOp_start<exon_start-min_overlap and cigarOp_end>=exon_start-min_overlap):
                if (cigarOp_end>=exon_start+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(exon_start-min_overlap)+1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=exon_start-min_overlap and cigarOp_end<exon_start+min_overlap):
                count=op_Length
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary       

            elif (cigarOp_start<exon_start+min_overlap and cigarOp_end>=exon_start+min_overlap):
                if (cigarOp_start<=exon_start-min_overlap):
                    count=min_overlap*2
                else:
                    count=(exon_start+min_overlap)-cigarOp_start-1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string is within the exon
            if (cigarOp_start<exon_start and cigarOp_end>=exon_start):
                if (cigarOp_end>=exon_end):
                    count=exon_end-exon_start
                else:
                    count=cigarOp_end-exon_start
                exon_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=exon_start and cigarOp_end<exon_end):
                count=op_Length
                exon_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<exon_end and cigarOp_end>=exon_end):
                if (cigarOp_start<=exon_start):
                    count=exon_end-exon_start
                else:
                    count=exon_end-cigarOp_start
                exon_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string spans the designated exon end
            if (cigarOp_start<exon_end-min_overlap and cigarOp_end>=exon_end-min_overlap):
                if (cigarOp_end>=exon_end+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(exon_end-min_overlap)
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=exon_end-min_overlap and cigarOp_end<exon_end+min_overlap):
                count=op_Length
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<exon_end+min_overlap and cigarOp_end>=exon_end+min_overlap):
                if (cigarOp_start<=exon_end-min_overlap):
                    count=min_overlap*2
                else:
                    count=(exon_end+min_overlap)-cigarOp_start
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

        # get length of the aligned portion of this read from cigar string
        aligned_read_length = CigarOp_counts['M']+CigarOp_counts['D']

        # get 5'SS and 3'SS counts as determined by gene strand
        strand = gene_strand
        if (strand == '+'):
            exon_start_counts = start_counts # record the cigar string counts over the 5'SS
            exon_end_counts = end_counts # record the cigar string counts over the 3'SS
            
        if (strand == '-'):
            exon_start_counts = end_counts # record the cigar string counts over the 5'SS
            exon_end_counts = start_counts # record the cigar string counts over the 3'SS  
            
            
        # annotate splicing status based on CIGAR string information around splice sites
        splice='UNDETERMINED'

        if (exon_start_counts['N']==0 and exon_end_counts['N']==0):
            if (exon_end_counts['M']+exon_end_counts['D']==min_overlap*2) and (exon_start_counts['M']+exon_start_counts['D']==min_overlap*2):
                if (exon_end_counts['M']>min_overlap) and (exon_start_counts['M']>min_overlap):
                    # this means the exon is included and surrounded by two unspliced exons, but this could still result in skipping afterwards
                    splice = 'UNSPLICED'
        
        if (exon_start_counts['M']>0 and exon_start_counts['M']<min_overlap*2):
            if (exon_end_counts['M']>0 and exon_end_counts['M']<min_overlap*2):
                # this means the exon is included and surrounded by two spliced exons
                splice = 'INCL'
            if (exon_end_counts['N']==0 and (exon_end_counts['M']+exon_end_counts['D']==min_overlap*2)):
                if (exon_end_counts['M']>min_overlap):
                    # this means the exon is included and the next exon is unspliced
                    splice = 'INCL'
            if (exon_end_counts['N']==0 and (exon_end_counts['M']+exon_end_counts['D']==0)):
                splice = 'INCLLAST' # this means the exon is (partially) included but the read ends within the exon
        
        if (exon_end_counts['M']>0 and exon_end_counts['M']<min_overlap*2):
            if (exon_start_counts['N']==0 and (exon_start_counts['M']+exon_start_counts['D']==min_overlap*2)):
                if (exon_start_counts['M']>min_overlap):
                    # this means the exon is included and the previous exon is unspliced
                    splice = 'INCL'
            if (exon_start_counts['N']==0 and (exon_start_counts['M']+exon_start_counts['D']==0)):
                splice = 'INCLFIRST' # this means the exon is (partially) included but the read ends within the exon
            
            
        if (exon_start_counts['M']==0 and exon_end_counts['M']==0):
            if (exon_end_counts['N']==min_overlap*2) and (exon_start_counts['N']==min_overlap*2):
                # this would be an excluded exon
                splice = 'EXCL'

        # annotate splicing status based on CIGAR string information within the exon 
        if (splice == 'EXCL'):
            if (float(exon_end-exon_start) > 0.0):
                ratio = float(exon_counts['N'])/float(exon_end-exon_start)
                difference = abs(exon_counts['N']-(exon_end-exon_start))

                # if exon is excluded, between 90-100% of the exon has to be excluded 
                # and no more than 10 nucleotides within the exon can be matching the exon sequence
                if( ratio < 0.9 or ratio > 1.1 or difference > 10):
                    splice='UNDETERMINED'
            if (float(exon_end-exon_start) == 0.0):
                splice='UNDETERMINED'

        if (splice == 'INCL') or (splice == "UNSPLICED"):
            if (float(exon_end-exon_start) > 0.0):
                ratio = float(exon_counts['M'])/(float(exon_counts['M'])+float(exon_counts['N'])+float(exon_counts['D'])+1)

                # if exon is included, at least 75% of the read has to match (CIGAR=M) the exon sequence
                if(exon_counts['M'] < min_overlap/2 or ratio < 0.75):
                    splice='UNDETERMINED'
            
            if (float(exon_end-exon_start) == 0.0):
                splice='UNDETERMINED'

        # save read, exon, and splicing information
        spliceCalls.append([aln_name,exon_chr,exon_start,exon_end,gene_strand,name_gene,exon_count,splice])

    spliceCalls_df = pd.DataFrame(spliceCalls)
    spliceCalls_df.columns = ["read_name","chrom","exon_start","exon_end","strand","gene_name","exon_count","splice_status"]

    return spliceCalls_df



# every read that spans an exon in the dataset
def get_read_junctions_exon_dictionary(splice_df):

    read_junctions = {}

    for i in range(0,splice_df.shape[0]):       

        # define the read name
        read_name = splice_df['read_name'].iloc[i]
        gene_name = splice_df['gene_name'].iloc[i]
        chrom = splice_df['chrom'].iloc[i]
        exon_start = splice_df['exon_start'].iloc[i]
        exon_end = splice_df['exon_end'].iloc[i]
        exon_count = splice_df['exon_count'].iloc[i]
        strand = splice_df['strand'].iloc[i]
        splice_status = splice_df['splice_status'].iloc[i]

        # check if read name is in the dictionary, if not save it
        if read_name not in read_junctions.keys():

            # make a new dictionary for the gene and add exon info to it
            read_junctions[read_name] = {}
            read_junctions[read_name][gene_name] = [[chrom, exon_start, exon_end, exon_count, strand, splice_status]]

        # check if read name is in the dictionary, if it is proceed to gene information
        elif read_name in read_junctions.keys():

            # if gene_name is not already in read dictionary, 
            # make a new dictionary for the gene and add exon info to it
            if gene_name not in read_junctions[read_name].keys():
                read_junctions[read_name][gene_name] = [[chrom, exon_start, exon_end, exon_count, strand, splice_status]]

            # if gene_name is already in read dictionary, add new exon info to it
            elif gene_name in read_junctions[read_name].keys():
                read_junctions[read_name][gene_name].append([chrom, exon_start, exon_end, exon_count, strand, splice_status])

    return read_junctions



##############################################

# set all variables for analysis
min_overlap = 25


# get splicing information for every read that spans an exon
print("getting splicing info...")
splice_info = get_splicing_info_exons(intersect_df,min_overlap)

# apply some filters and sort the dataframe
splice_info = splice_info.sort_values(by=['chrom','exon_start','exon_end']).reset_index(drop=True)
splice_info_nodup = splice_info.drop_duplicates(subset=['read_name','chrom','exon_start','exon_end','strand']).reset_index(drop=True)
splice_info_nodup.to_csv(out_splicing_info, sep='\t', index=False, header=True)

# get dictionary with all exon junctions that a read spans
# do for each chromosome, otherwise it runs into memory problems
chrom_list = splice_info['chrom'].drop_duplicates().tolist()

for c in chrom_list:
    splice_info_sub = splice_info[splice_info['chrom'] == c].reset_index(drop=True)
    print("making splicing dictionary for chr", c)
    splice_dictionary = get_read_junctions_exon_dictionary(splice_info_sub)
    out_sub_dict = out_splicing_dict + '_chr' + str(c) + '.npy'
    np.save(out_sub_dict, splice_dictionary)

