"""

Author: Karine Choquet

Date: Feb 16, 2021

This script will concatenate the splicing status of all introns within a transcript for each read.

Usage: python get_multi_introns_isoforms_df.directRNA.py splicing_dict multi_introns_df multi_introns_count

splicing_dict: dictionary (key: read) of dictionaries (key: transcript, value: list of introns and splicing statuses)
multi_introns_df: splicing status of each read per transcript, one read per line. Splicing status of each intron is separated by an underscore
multi_introns_count: number of reads per splicing status and per transcript


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
out_multi_introns_df = sys.argv[2]
out_multi_introns_counts = sys.argv[3]

######################################


def get_all_multi_introns(read_junctions, intron_min):

    multi_introns_list = []

    for read in read_junctions.keys():

        # make a set for all intron pairs within a read
        # this will avoid duplicate pairs being called due to alternative splicing
        uniq_splice_pattern = set()

        # loop through all genes that has introns that a read maps to
        for gene in read_junctions[read].keys():

            # only go through genes that have 2 or more introns
            if len(read_junctions[read][gene]) >= intron_min:

                # characterize the number of spliced and unspliced introns in the read
                strand = [row[4] for row in read_junctions[read][gene]][0]
                splice_status = [row[6] for row in read_junctions[read][gene]]
                intron_numbers = [row[3] for row in read_junctions[read][gene]]
                splice_status_join = '_'.join(splice_status)
                intron_numbers_join = '_'.join([str(i) for i in intron_numbers])
                status_count = Counter(splice_status)
                multi_introns_list.append([read,gene,strand,intron_numbers_join,splice_status_join])
                
    multi_introns_df = pd.DataFrame(multi_introns_list)
    
    if len(multi_introns_df)>0:
        multi_introns_df.columns = ['read','gene','strand','intron_numbers','splice_status']
    else:
        multi_introns_df.columns = pd.DataFrame(['read','gene','strand','intron_numbers','splice_status'])

    return multi_introns_df




##############################################

# set all variables for analysis
min_overlap = 25

splice_dictionary = np.load(splicing_dict_file, encoding='latin1', allow_pickle=True).item()
# Get multi-introns dataset
multi_introns_df = get_all_multi_introns(splice_dictionary, 2)
multi_introns_df.to_csv(out_multi_introns_df, sep='\t', index=False, header=True)


# Summarize and count
if len(multi_introns_df)>0:
    multi_introns_counts = pd.DataFrame(multi_introns_df.groupby(['gene','strand','intron_numbers','splice_status'])['read'].count()).reset_index()
    multi_introns_counts.columns = ['gene','strand','intron_numbers','splice_status','count']
    multi_introns_counts.to_csv(out_multi_introns_counts, sep='\t', index=False, header=True)

