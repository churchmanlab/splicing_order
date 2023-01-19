"""

Author: Karine Choquet

Date: March 11, 2021

This script will determine the splicing status of every possible intron/exon pair within nanopore reads in candidate genes

Usage: python get_splice_status_per_intron_exon_pair_directRNA.py ctrl_multi_introns_df U2_multi_introns_df ctrl_multi_exons_df U2_multi_exons_df out_intron_exon_pairs candidate_genes

*_multi_introns_df and *_multi_exons_df are the output files from get_splice_status_introns_direct_RNA_nanopore_seq.sh and get_splice_status_exons_direct_RNA_nanopore_seq.sh
candidate_genes is a list of genes in which both SE and RI were detected by short-read RNA-seq

"""

import sys
import numpy as np
import pandas as pd

import scipy
from scipy import stats

import itertools

from statsmodels.stats.multitest import multipletests

from more_itertools import consecutive_groups

#####################################

# CONFIG

print("Reading input files")

# Splicing dataframes for each sample
ctrl_multi_introns_df = pd.read_table(sys.argv[1])
U2_multi_introns_df = pd.read_table(sys.argv[2])
ctrl_multi_exons_df = pd.read_table(sys.argv[3])
U2_multi_exons_df = pd.read_table(sys.argv[4])

# Output files
out_intron_exon_pairs_counts = sys.argv[5]

# FLAIR significant isoforms or other candidate genes
sig_genes = pd.read_table(sys.argv[6], header=None)
sig_genes.columns = ['gene']
sig_genes_list = sig_genes['gene'].tolist()
print(sig_genes_list)

# Merge the intron and exon dataframes by read name
print("Merging intron and exon files")
U2_multi_both_df = U2_multi_introns_df.merge(U2_multi_exons_df, on=['read','gene','strand'])
ctrl_multi_both_df = ctrl_multi_introns_df.merge(ctrl_multi_exons_df, on=['read','gene','strand'])

U2_multi_both_df.columns = ['read','gene','strand','intron_numbers','intron_statuses','exon_numbers','exon_statuses']
ctrl_multi_both_df.columns = ['read','gene','strand','intron_numbers','intron_statuses','exon_numbers','exon_statuses']


######################################

# retrieve splicing status for each possible pair of introns
def get_splice_status_per_intron_exon_pair(multi_both_df):
    pairs_list = []

    for i in range(len(multi_both_df)):
        readname = multi_both_df.loc[i]['read']
        gene_id = multi_both_df.loc[i]['gene']
        strand = multi_both_df.loc[i]['strand']
    
        intron_numbers = multi_both_df.loc[i]['intron_numbers']
        intron_splice_status = multi_both_df.loc[i]['intron_statuses']
        exon_numbers = multi_both_df.loc[i]['exon_numbers']
        exon_splice_status = multi_both_df.loc[i]['exon_statuses']
    
        intron_list = [int(a) for a in intron_numbers.split("_")]
        intron_status_list = [b for b in intron_splice_status.split("_")]
        exon_list = [int(a) for a in exon_numbers.split("_")]
        exon_status_list = [b for b in exon_splice_status.split("_")]
    
        n_introns_SKP = len([c for c in intron_status_list if "SKP" in c]) # SKIP reads where all junctions are skipped
    
        # Go back to the beginning of the loop if all the introns in the read are SKP
        if n_introns_SKP == len(intron_status_list):
            continue
    
        # Iterate through every possible intron/exon combination and retrieve the splicing statuses
        for pair in itertools.product(intron_list, exon_list):
            my_int = pair[0]
            my_ex = pair[1]
            if my_int != my_ex and my_int+1 != my_ex: # no need to look at combinations that overlap
                index_int = intron_list.index(my_int)
                index_ex = exon_list.index(my_ex)
                int_splice = intron_status_list[index_int]
                ex_splice = exon_status_list[index_ex]
                
                # Add a special status for the exons that are first or last in a read, since they are undetermined 
                # only because the read doesn't span both the start and end
                if ex_splice == "UNDETERMINED" and index_ex == 0:
                    ex_splice = "UND-START-READ"
                if ex_splice == "UNDETERMINED" and index_ex == exon_list.index(exon_list[-1]):
                    ex_splice = "UND-END-READ"
        
                pair_id = gene_id + "," + str(my_int) + "," + str(my_ex)
                pair_splice_status = int_splice + "_" + ex_splice
        
                pairs_list.append([readname, pair_id, pair_splice_status])
        
    pairs_df = pd.DataFrame(pairs_list)
    pairs_df.columns = ['read','pair_id','splice_status']

    pairs_counts = pd.DataFrame(pairs_df.groupby(['pair_id','splice_status'])['read'].count()).reset_index()
    
    return pairs_counts


# Compare alternative counts and reference counts (YES_INCL) in control and U2
# Even if the intron is later than the exon in the gene, the intron splice status is always first in the splice status
def fishers_exact_test_per_row(row, df_yes_incl):
    pair_id = row['pair_id']
    splice_status = row['splice_status']
    count_ctrl_alt = row['count_ctrl']
    count_U2_alt = row['count_U2']
    
    # Retrieve the corresponding YES_INCL count
    try:
        count_ctrl_ref = int(df_yes_incl.loc[df_yes_incl['pair_id'] == pair_id]['count_ctrl'])
    except TypeError:
        count_ctrl_ref = 0
        
    try:
        count_U2_ref = int(df_yes_incl.loc[df_yes_incl['pair_id'] == pair_id]['count_U2'])
    except TypeError:
        count_U2_ref = 0
    
    p = stats.fisher_exact([[count_ctrl_ref, count_U2_ref], [count_ctrl_alt, count_U2_alt]])[1]
        
    return(p)


def test_intron_exon_pairs(HeLa_intron_exon_pairs_counts):
    # For each intron pair that showed significant differences in the U2 and control, test for co-occurrence of the two
    # events in the U2 KD

    results_list = []

    for i in range(len(HeLa_intron_exon_pairs_counts)):
        pair_id = HeLa_intron_exon_pairs_counts.loc[i]['pair_id']
        splice_status = HeLa_intron_exon_pairs_counts.loc[i]['splice_status']
        count_U2 = HeLa_intron_exon_pairs_counts.loc[i]['count_U2']
        count_ctrl = HeLa_intron_exon_pairs_counts.loc[i]['count_ctrl']
        pvalue = HeLa_intron_exon_pairs_counts.loc[i]['pvalue']
        int_splice = HeLa_intron_exon_pairs_counts.loc[i]['int_splice']
        ex_splice = HeLa_intron_exon_pairs_counts.loc[i]['ex_splice']
    
        if pvalue < 0.05:
            if count_U2 > 10:
        
                # Test co-occurrence of SE/RI in U2 KD:
                if int_splice == "NO" and ex_splice == "EXCL":
            
                    # Test dependence of RI on SE
                    try:
                        int_INCL_U2 = int(HeLa_intron_exon_pairs_counts.loc[(HeLa_intron_exon_pairs_counts['pair_id'] == pair_id) &
                                                                (HeLa_intron_exon_pairs_counts['int_splice'] == int_splice) & 
                                                                (HeLa_intron_exon_pairs_counts['ex_splice'] == "INCL")]['count_U2'].sum())
                    except TypeError:
                        int_INCL_U2 = 0
                    
                    try:
                        int_INCL_ctrl = int(HeLa_intron_exon_pairs_counts.loc[(HeLa_intron_exon_pairs_counts['pair_id'] == pair_id) &
                                                                (HeLa_intron_exon_pairs_counts['int_splice'] == int_splice) & 
                                                                (HeLa_intron_exon_pairs_counts['ex_splice'] == "INCL")]['count_ctrl'].sum())
                    except TypeError:
                        int_INCL_ctrl = 0
            
                
            
                    # Test dependence of SE on RI in U2 KD:
                    try:
                        YES_ex_U2 = int(HeLa_intron_exon_pairs_counts.loc[(HeLa_intron_exon_pairs_counts['pair_id'] == pair_id) &
                                                                (HeLa_intron_exon_pairs_counts['int_splice'] == "YES") & 
                                                                (HeLa_intron_exon_pairs_counts['ex_splice'] == ex_splice)]['count_U2'].sum())
                    except TypeError:
                        YES_ex_U2 = 0
                    
                    try:
                        YES_ex_ctrl = int(HeLa_intron_exon_pairs_counts.loc[(HeLa_intron_exon_pairs_counts['pair_id'] == pair_id) &
                                                                (HeLa_intron_exon_pairs_counts['int_splice'] == "YES") & 
                                                                (HeLa_intron_exon_pairs_counts['ex_splice'] == ex_splice)]['count_ctrl'].sum())
                    except TypeError:
                        YES_ex_ctrl = 0
            

                    # Retrieve the number of normally spliced reads
                    try:
                        YES_INCL_U2 = int(HeLa_intron_exon_pairs_counts.loc[(HeLa_intron_exon_pairs_counts['pair_id'] == pair_id) &
                                                                (HeLa_intron_exon_pairs_counts['int_splice'] == "YES") &
                                                                (HeLa_intron_exon_pairs_counts['ex_splice'] == "INCL")]['count_U2'].sum())
                    except TypeError:
                        YES_INCL_U2 = 0

                    try:
                        YES_INCL_ctrl = int(HeLa_intron_exon_pairs_counts.loc[(HeLa_intron_exon_pairs_counts['pair_id'] == pair_id) &
                                                                (HeLa_intron_exon_pairs_counts['int_splice'] == "YES") &
                                                                (HeLa_intron_exon_pairs_counts['ex_splice'] == "INCL")]['count_ctrl'].sum())
                    except TypeError:
                        YES_INCL_ctrl = 0

            
                    results_list.append([pair_id, splice_status, "SE_RI", int_INCL_U2, YES_ex_U2,
                                          int_INCL_ctrl, YES_ex_ctrl, YES_INCL_U2,YES_INCL_ctrl])
                
            
        
    if len(results_list) >= 1:
        results_df = pd.DataFrame(results_list)
        results_df.columns = ['pair_id','splice_status','event_types','NO_INCL_U2','YES_EXCL_U2',
                                   'NO_INCL_ctrl','YES_EXCL_ctrl','YES_INCL_U2','YES_INCL_ctrl']

        
        return(results_df)
        
    
    


def compare_intron_exon_pairs_between_conditions(U2_test_df, ctrl_test_df):
    # Retrieve splicing status for each intron pair
    U2_intron_exon_pairs_counts = get_splice_status_per_intron_exon_pair(U2_test_df)
    ctrl_intron_exon_pairs_counts = get_splice_status_per_intron_exon_pair(ctrl_test_df)
    
    # Merge the counts from the two conditions
    HeLa_intron_exon_pairs_counts = pd.merge(ctrl_intron_exon_pairs_counts, U2_intron_exon_pairs_counts, on=['pair_id','splice_status'], how='outer').fillna(0)
    HeLa_intron_exon_pairs_counts.columns = ['pair_id','splice_status','count_ctrl','count_U2']

    # Make a dataframe that contains the counts for spliced/included intron/exon pairs
    # Even when the exon is upstream of the intron in the gene, the order of the splice status is always intron-exon
    # So the status wanted will always be YES_INCL
    HeLa_df_yes_incl = HeLa_intron_exon_pairs_counts[HeLa_intron_exon_pairs_counts['splice_status']=='YES_INCL'].reset_index(drop=True)
        
    # Apply Fisher's exact test to detect differences between the two conditions, followed by FDR correction
    HeLa_intron_exon_pairs_counts['pvalue'] = HeLa_intron_exon_pairs_counts.apply(lambda row: fishers_exact_test_per_row(row, HeLa_df_yes_incl), axis=1)
    #HeLa_intron_exon_pairs_counts['FDR'] = multipletests(HeLa_intron_exon_pairs_counts['pvalue'], alpha=0.05, method='fdr_bh')[1]

    # Make a separate column for the splicing status of each intron
    HeLa_intron_exon_pairs_counts['int_splice'] = HeLa_intron_exon_pairs_counts['splice_status'].str.split("_").str[0]
    HeLa_intron_exon_pairs_counts['ex_splice'] = HeLa_intron_exon_pairs_counts['splice_status'].str.split("_").str[1]
    
    # Apply the binomial test for each intron pair
    results_df = test_intron_exon_pairs(HeLa_intron_exon_pairs_counts)
    
    # Merge back with the original dataframe that contains the counts
    try:
        HeLa_intron_exon_pairs_counts_binom = HeLa_intron_exon_pairs_counts.merge(results_df, on=['pair_id','splice_status'])
    
        return(HeLa_intron_exon_pairs_counts_binom)
    
    except TypeError:
        pass
    

#####################################


# Retrieve genes with significant differences between control and U2 KD

for gene in sig_genes_list:

    print(gene)   
 
    U2_test_df = U2_multi_both_df[U2_multi_both_df['gene']==gene].reset_index(drop=True)
    ctrl_test_df = ctrl_multi_both_df[ctrl_multi_both_df['gene']==gene].reset_index(drop=True)
    
    if len(U2_test_df) > 1 and len(ctrl_test_df) > 1:
        U2_vs_ctrl_df = compare_intron_exon_pairs_between_conditions(U2_test_df, ctrl_test_df)
        
        if U2_vs_ctrl_df is not None:
            try:
                results_df = pd.concat([results_df,U2_vs_ctrl_df]).reset_index(drop=True)
        
            except NameError:
                results_df = U2_vs_ctrl_df



# Write to file
results_df.to_csv(out_intron_exon_pairs_counts, sep="\t", header=True, index=False) 
    


