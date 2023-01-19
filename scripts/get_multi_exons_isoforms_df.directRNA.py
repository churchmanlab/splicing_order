"""

Author: Karine Choquet

Date: March 19, 2021

This script will concatenate the splicing status of all exons within a transcript for each read.

Usage: python get_multi_introns_isoforms_df.directRNA.py splicing_dict multi_introns_df multi_introns_count

splicing_dict: dictionary (key: read) of dictionaries (key: transcript, value: list of exons and splicing statuses)
multi_exons_df: splicing status of each read per transcript, one read per line. Splicing status of each exon is separated by an underscore
multi_exons_count: number of reads per splicing status and per transcript


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


# Splicing dictionary for that chromosome
splicing_dict_file = sys.argv[1]

# Output files
out_multi_exons_df = sys.argv[2]
out_multi_exons_counts = sys.argv[3]

######################################


def get_all_multi_exons(read_junctions, exon_min):

    multi_exons_list = []

    for read in read_junctions.keys():

        # make a set for all exon pairs within a read
        # this will avoid duplicate pairs being called due to alternative splicing
        uniq_splice_pattern = set()

        # loop through all genes that has exons that a read maps to
        for gene in read_junctions[read].keys():

            # only go through genes that have 2 or more exons
            if len(read_junctions[read][gene]) >= exon_min:

                # characterize the number of spliced and unspliced exons in the read
                strand = [row[4] for row in read_junctions[read][gene]][0]
                splice_status = [row[5] for row in read_junctions[read][gene]]
                exon_numbers = [row[3] for row in read_junctions[read][gene]]
                splice_status_join = '_'.join(splice_status)
                exon_numbers_join = '_'.join([str(i) for i in exon_numbers])
                status_count = Counter(splice_status)
                multi_exons_list.append([read,gene,strand,exon_numbers_join,splice_status_join])
                
    multi_exons_df = pd.DataFrame(multi_exons_list)
    
    if len(multi_exons_df)>0:
        multi_exons_df.columns = ['read','gene','strand','exon_numbers','splice_status']
    else:
        multi_exons_df.columns = pd.DataFrame(['read','gene','strand','exon_numbers','splice_status'])

    return multi_exons_df


##############################################

# set all variables for analysis
min_overlap = 25

splice_dictionary = np.load(splicing_dict_file, encoding='latin1', allow_pickle=True).item()
# Get multi-exons dataset
multi_exons_df = get_all_multi_exons(splice_dictionary, 2)
multi_exons_df.to_csv(out_multi_exons_df, sep='\t', index=False, header=False)

# Summarize and count
if len(multi_exons_df)>0:
    multi_exons_counts = pd.DataFrame(multi_exons_df.groupby(['gene','strand','exon_numbers','splice_status'])['read'].count()).reset_index()
    multi_exons_counts.columns = ['gene','strand','exon_numbers','splice_status','count']
    multi_exons_counts.to_csv(out_multi_exons_counts, sep='\t', index=False, header=False)

