"""

Author: Karine Choquet

Date: June 4, 2021

This script will identify possible splicing paths in IFRD2 and assign them a score based on the counts from chromatin RNA sequencing, allowing for SE events

Usage: ython get_splicing_paths.cDNA-PCR.with_SE.py gene_names_df intron_df ctrl_1_multi_introns_counts ctrl_2_multi_introns_counts U2_1_multi_introns_counts U2_2_multi_introns_counts splicing_paths_with_SE

gene_names_df is an annotation file containing transcript ids (NM_) in one column and the gene name in the other
intron_df is a BED file with coordinates of each intron of interest
*_multi_intron_counts are output files from get_splice_status_introns_cDNA_PCR_nanopore_seq.sh containing the splicing statuses of each intron and each read


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

import itertools
from more_itertools import consecutive_groups

#####################################

# CONFIG

# genes names info
gene_names_df = pd.read_table(sys.argv[1], header=None)
gene_names_df.columns = ['gene_name','gene_id']

# get intron info to get the true intron numbers for transcripts on the minus strand
gene_id_list = ['NM_006764.4']

hg38_intron_info = pd.read_table(sys.argv[2])
hg38_intron_info_sub = hg38_intron_info[hg38_intron_info['gene'].isin(gene_id_list)].reset_index(drop=True)


# Multi-intron counts from each sample
sample1 = pd.read_table(sys.argv[3])
sample2 = pd.read_table(sys.argv[4])
sample3 = pd.read_table(sys.argv[5])
sample4 = pd.read_table(sys.argv[6])

sample1['target'] = "ctrl_kd_1"
sample2['target'] = "ctrl_kd_2"
sample3['target'] = "U2_kd_1"
sample4['target'] = "U2_kd_2"

all_multi_introns_df = pd.concat([sample1,sample2,sample3,sample4]).reset_index(drop=True)
all_multi_introns_df = all_multi_introns_df[~all_multi_introns_df['splice_status'].str.contains("UND")].reset_index(drop=True)

all_multi_introns_df['gene_id'] = all_multi_introns_df['gene'].str.split(".").str[0]
all_multi_introns_df = all_multi_introns_df.merge(gene_names_df, on='gene_id')
all_multi_introns_df = all_multi_introns_df.drop('gene_id', axis=1)

hg38_intron_info_sub = hg38_intron_info_sub.merge(all_multi_introns_df[['gene','gene_name']],on='gene')

# Output file
out_splicing_paths = sys.argv[-1]

######################################

def get_score_per_path_pair_v2(multi_introns_df, total_introns, introns_of_interest, n_introns_gene, strand):
    

    # Make a dictionary with the patterns and the counts and another with the number of introns spliced and the counts
    pattern_dict = {}
    results_list = []
    spliced_counts_dict = {}
    
    introns_of_interest_list_tmp = introns_of_interest.split("_")
    if strand == "-":
        introns_of_interest_list = [str(n_introns_gene - int(i)) for i in introns_of_interest_list_tmp]
        introns_of_interest_fix = "_".join(introns_of_interest_list)
    elif strand == "+":
        introns_of_interest_list = [str(int(i)+1) for i in introns_of_interest_list_tmp]
        introns_of_interest_fix = "_".join(introns_of_interest_list)
    
    # Initiate pattern_dict
    pattern_dict = {}
    for n in range(len(introns_of_interest_list)+1):
        pattern_dict[n] = {}
        
    
    # Iterate to get counts for isoforms
    for row in range(len(multi_introns_df)):
        gene_name = multi_introns_df.loc[row]['gene_name']
        splice_status_temp = multi_introns_df.loc[row]['splice_status']
        intron_numbers_list_tmp = multi_introns_df.loc[row]['intron_numbers'].split("_")
        if strand == "-":
            intron_numbers_list = [str(n_introns_gene - int(i)) for i in intron_numbers_list_tmp]
            intron_numbers = "_".join(intron_numbers_list)
        elif strand == "+":
            intron_numbers_list = [str(int(i)+1) for i in intron_numbers_list_tmp]
            intron_numbers = "_".join(intron_numbers_list)
        pattern_count = multi_introns_df.loc[row]['count']
        sample_name = multi_introns_df.loc[row]['target']
        splice_status_list_temp = splice_status_temp.split("_")
        splice_status_list = ["SKP" if "SKP" in a else a for a in splice_status_list_temp]
        splice_status = "_".join(splice_status_list)
        if len(splice_status_list) == total_introns and introns_of_interest_fix == intron_numbers:
            skipped_count = Counter(splice_status_list)['SKP']
            skipped_count_3SS = Counter(splice_status_list_temp)['SKP3SS']
            skipped_count_5SS = Counter(splice_status_list_temp)['SKP5SS']
            spliced_count = Counter(splice_status_list)['YES']
            unspliced_count = Counter(splice_status_list)['NO']
            if skipped_count_3SS == skipped_count_5SS and skipped_count_3SS == 1: # this ensures that we are looking at truly skipped exons on both sides, and not some artefacts
                true_skipped_count = skipped_count_3SS
                level = spliced_count + true_skipped_count
                if skipped_count < total_introns:
                    if splice_status not in pattern_dict[level].keys():
                        pattern_dict[level][splice_status] = pattern_count
                    elif splice_status in pattern_dict[level].keys():
                        pattern_dict[level][splice_status] += pattern_count
                    if level not in spliced_counts_dict.keys():
                        spliced_counts_dict[level] = pattern_count
                    else:
                        spliced_counts_dict[level] += pattern_count
            
            elif skipped_count == 0: # no skipped introns
                level = spliced_count
                if skipped_count < total_introns:
                    if splice_status not in pattern_dict[level].keys():
                        pattern_dict[level][splice_status] = pattern_count
                    elif splice_status in pattern_dict[level].keys():
                        pattern_dict[level][splice_status] += pattern_count
                    if level not in spliced_counts_dict.keys():
                        spliced_counts_dict[level] = pattern_count
                    else:
                        spliced_counts_dict[level] += pattern_count           
            
    if len(pattern_dict[0]) == 0:
        unspliced_iso = "_".join(["NO" for i in range(len(introns_of_interest_list))])
        pattern_dict[0][unspliced_iso] = 0
    
    print(pattern_dict)
    
    # Initiate path_dict
    path_dict = {}
    for n in range(len(introns_of_interest_list)):
        path_dict[n] = {}
        
    # For each combination of isoforms between levels (e.g. number of splicing events),
    # retrieve the frequencies and calculate scores
    # Below, the following scores are calculated:
    # freq: number of reads for that isoform / total number of reads for that level
    # sub_path_score : freq for that isoform * freq for the isoform from which it is derived
    # path_score: freq for that isoform * path_score from the previous isoform (so multiplicative frequencies)
    # full_path_score: the path_score for the last level of that path
    for levels in itertools.product(pattern_dict.keys(), pattern_dict.keys()):
        level1 = levels[0]
        level2 = levels[1]
        if level1 != level2: # we don't want to compare isoforms from the same level
            
            # Iterate through each pair of isoforms that are from different levels
            for pair in itertools.product(pattern_dict[level1].keys(), pattern_dict[level2].keys()):
                pattern1 = pair[0].split("_")
                pattern2 = pair[1].split("_")
        
                # Get the splicing level of the isoform
                unspliced_introns_pattern1 = len([i for i, x in enumerate(pattern1) if x == "NO"])
                unspliced_introns_pattern2 = len([i for i, x in enumerate(pattern2) if x == "NO"])
                spliced_introns_pattern2 = Counter(pattern2)['YES'] + Counter(pattern2)['SKP']
                
                # Define the level that will be used below
                level = level1
        
                # retrieve the positions of the difference(s) between the two patterns
                diff_index_list = [i for i, x in enumerate(pattern2) if pattern1[i]!=x]
        
                # If pattern2 has one more spliced intron than pattern1:
                if len(diff_index_list) == 1:
                    diff_index = diff_index_list[0]
                    if pattern1[diff_index] == "NO" and pattern2[diff_index] == "YES":
                        count_pattern1 = pattern_dict[level1][pair[0]]
                        count_pattern2 = pattern_dict[level2][pair[1]]
                
                        if level == 0: # this means we are starting from the unspliced isoform
                            # define a new splicing path
                            new_intron_spliced = str(introns_of_interest_list[diff_index])
                            path_name = new_intron_spliced + "->"
                            freq = count_pattern2 / spliced_counts_dict[level+1]
                            path_score = freq
                            sub_path_score = freq
                            path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "spliced", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]
                    
                        elif level > 0 and level < len(introns_of_interest_list)-1: # this means we are at an intermediate isoform
                            for k in path_dict[level-1].keys():
                                if path_dict[level-1][k][1] == pattern1:
                                    new_intron_spliced = str(introns_of_interest_list[diff_index])
                                    if unspliced_introns_pattern2 == 0:
                                        path_name = str(k) + new_intron_spliced
                                    elif unspliced_introns_pattern2 > 0:
                                        path_name = str(k) + new_intron_spliced + "->"
                                    freq = count_pattern2 / spliced_counts_dict[level+1]
                                    path_score = path_dict[level-1][k][-1] * freq
                                    sub_path_score = path_dict[level-1][k][-3] * freq # only the frequency from the previous level and this level
                                    path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "spliced", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]
                            
                        elif level == len(introns_of_interest_list)-1: # this means we are at the fully spliced isoform        
                            for k in path_dict[level-1].keys():
                                if path_dict[level-1][k][1] == pattern1:
                                    new_intron_spliced = str(introns_of_interest_list[diff_index])
                                    path_name = str(k) + new_intron_spliced
                                    freq = 1
                                    path_score = path_dict[level-1][k][-1] * freq
                                    sub_path_score = path_dict[level-1][k][-3] * freq
                                    path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "spliced", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]
                        
                # If pattern2 has one more skipped event that includes more than one intron:
                # With the intron-centric annotations, this is necessarily the case since it requires one intron
                # to be SKP3SS and one intron to be SKP5SS
                elif len(diff_index_list) > 1:
                    count_diff = 0
                    skipped_indices_list = []
                    for a in diff_index_list:
                        if pattern1[a] == "NO" and pattern2[a] == "SKP":
                            if a > 0 and a < len(pattern2)-1:
                                if pattern1[a-1] != "SKP" and pattern1[a+1] != "SKP": 
                                    count_diff +=1
                                    skipped_indices_list.append(a)
                            elif a == 0:
                                if pattern1[a+1] != "SKP":
                                    count_diff +=1
                                    skipped_indices_list.append(a)
                            elif a == len(pattern1)-1:
                                if pattern1[a-1] != "SKP":
                                    count_diff += 1
                                    skipped_indices_list.append(a)
                    
                    # test if the number of new introns is the correct length and if they are all consecutive intron numbers
                    if count_diff == len(diff_index_list) and sorted(skipped_indices_list) == list(range(min(skipped_indices_list), max(skipped_indices_list)+1)):
                        count_pattern1 = pattern_dict[level1][pair[0]]
                        count_pattern2 = pattern_dict[level2][pair[1]]
                        if level == 0: # this means we are starting from the unspliced isoform
                            # define a new splicing path
                            n_new_intron_spliced = len(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                            new_intron_spliced = "_".join(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                            path_name = "[" + new_intron_spliced + "]->"
                            freq = count_pattern2 / spliced_counts_dict[level+1]
                            path_score = freq
                            sub_path_score = freq
                            path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "skipped", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]
                    
                        elif level > 0 and level < len(introns_of_interest_list)-1: # this means we are at an intermediate isoform
                            for k in path_dict[level-1].keys():
                                if path_dict[level-1][k][1] == pattern1:                                
                                    if unspliced_introns_pattern2 > 0:
                                        n_new_intron_spliced = len(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                                        new_intron_spliced = "_".join(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                                        path_name = str(k) + "[" + new_intron_spliced + "]->"
                                    elif unspliced_introns_pattern2 == 0:
                                        n_new_intron_spliced = len(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                                        new_intron_spliced = "_".join(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                                        path_name = str(k) + "[" + new_intron_spliced + "]"
                                    freq = count_pattern2 / spliced_counts_dict[level+1]
                                    path_score = path_dict[level-1][k][-1] * freq
                                    sub_path_score = path_dict[level-1][k][-3] * freq # only the frequency from the previous level and this level
                                    path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "skipped", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]
                            
                        elif level == len(introns_of_interest_list)-1: # this means we are at the fully spliced isoform
                            for k in path_dict[level-1].keys():
                                if path_dict[level-1][k][1] == pattern1:
                                    n_new_intron_spliced = len(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                                    new_intron_spliced = "_".join(introns_of_interest_list[diff_index_list[0]:diff_index_list[-1]+1])
                                    path_name = str(k) + "[" + new_intron_spliced
                                    freq = 1
                                    path_score = path_dict[level-1][k][-1] * freq
                                    sub_path_score = path_dict[level-1][k][-3] * freq # only the frequency from the previous level and this level
                                    path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "skipped", level, new_intron_spliced, spliced_counts_dict[new_level+1], freq, sub_path_score, path_score]
                            
    
    print(path_dict)

    # Now loop through the dictionary to match each pair with the possible final paths
    # The final level contains only final paths, so first retrieve those
    final_level = list(path_dict.keys())[-1]
    final_level_df = pd.DataFrame.from_dict(path_dict[final_level], orient='index').reset_index()
    final_level_df.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']
    final_level_df = final_level_df.drop(['pattern1_list','pattern2_list'],axis=1)
    final_level_df['full_path'] = final_level_df['path_name']
    
    final_df = final_level_df.copy()
    
    # Iterate through each of the levels to match the partial paths with all the possible final paths
    for lev in list(reversed(list(path_dict.keys())))[:-1]:
        # For the two final levels, merge the second last with the last to retrieve the final path and add the final path score
        if lev == final_level:
            df1 = pd.DataFrame.from_dict(path_dict[lev], orient='index').reset_index()
            df2 = pd.DataFrame.from_dict(path_dict[lev-1], orient='index').reset_index()
        
            df1.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']
            df2.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']
        
            fields = ['pattern1','path_name']
        
            new_df = df2.merge(df1[fields], left_on='pattern2', right_on='pattern1', how='left')
            new_df = new_df.rename(columns={'path_name_x':'path_name','path_name_y':'full_path','pattern1_x':'pattern1'}).drop(['pattern1_y','pattern1_list','pattern2_list'], axis=1)
            
            new_df = new_df.fillna(0)
            
            # If full_path is null, that means that it wasn't present in the level above and therefore the full path
            # is likely to be the current path, so replace it that way
            new_df_sub1 = new_df[new_df['full_path']==0].reset_index(drop=True)
            new_df_sub2 = new_df[new_df['full_path']!=0].reset_index(drop=True)
            
            new_df_sub1['full_path'] = new_df_sub1['path_name']
            
            new_df = pd.concat([new_df_sub1, new_df_sub2]).reset_index(drop=True)
        
            final_df = pd.concat([final_df, new_df]).reset_index(drop=True)
            
        # For any previous levels, repeat the previous steps
        elif (lev - 1) >= 0:
            df1 = new_df.copy()
            df2 = pd.DataFrame.from_dict(path_dict[lev-1], orient='index').reset_index()
        
            df2.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']
        
            fields = ['pattern1','full_path']
        
            new_df = df2.merge(df1[fields], left_on='pattern2', right_on='pattern1', how='left')
            new_df = new_df.rename(columns={'pattern1_x':'pattern1'}).drop(['pattern1_y','pattern1_list','pattern2_list'], axis=1)
            
            new_df = new_df.fillna(0)
            
            new_df_sub1 = new_df[new_df['full_path']==0].reset_index(drop=True)
            new_df_sub2 = new_df[new_df['full_path']!=0].reset_index(drop=True)
            
            new_df_sub1['full_path'] = new_df_sub1['path_name']
            
            new_df = pd.concat([new_df_sub1, new_df_sub2]).reset_index(drop=True)
        
            final_df = pd.concat([final_df, new_df]).reset_index(drop=True)
            
            
    # Now make sure that the start of the full path is the same as the path name, since merging on patterns as above
    # will give rows where that is not the case
    final_final_df = final_df[final_df.apply(lambda row: row.full_path.startswith(row.path_name), axis=1)].drop_duplicates().sort_values(by='level').reset_index(drop=True)
    
    
    # Get the final score for the path and express it so that the total score of all final isoforms is 1
    last_isos1 = final_final_df[~final_final_df['path_name'].str.endswith("->")][['full_path','path_score']].drop_duplicates().reset_index(drop=True)
    last_isos2 = final_final_df[(final_final_df['path_name']==final_final_df['full_path']) & (final_final_df['full_path'].str.endswith("->"))][['full_path','path_score']].drop_duplicates().reset_index(drop=True)
    last_isos = pd.concat([last_isos1, last_isos2]).reset_index(drop=True)
    last_isos['full_path_score'] = last_isos['path_score'] / last_isos['path_score'].sum()
    last_isos = last_isos.drop('path_score',axis=1).sort_values(by='full_path_score', ascending=False).reset_index(drop=True)
    last_isos['rank'] = last_isos.index + 1
    
    # Modify the levels and rows where the unspliced isoform is pattern2
    final_final_df['level'] = final_final_df['level'] + 1
    
    NO_df = final_final_df[final_final_df['level']==1].reset_index(drop=True)
    NO_df['pattern2'] = NO_df['pattern1']
    NO_df['count_pattern2'] = NO_df['count_pattern1']
    NO_df['level'] = 0
    NO_df['new_intron_spliced'] = 0
    
    final_final_df = pd.concat([final_final_df,NO_df]).sort_values(by='level').reset_index(drop=True)
    
    final_final_df = final_final_df.merge(last_isos, on='full_path')
               
    return(final_final_df)




# Apply the function above to all genes of interest in a given sample
def get_read_count_path_per_gene(multi_introns_df, genes_regions_dict, sample_list, hg38_intron_df):
    
    for gene_name in genes_regions_dict.keys():
        for sample_name in sample_list:
            print(gene_name, sample_name)
            df = multi_introns_df[(multi_introns_df['gene_name']==gene_name) & (multi_introns_df['target']==sample_name)].reset_index(drop=True)
            introns_of_interest = genes_regions_dict[gene_name]
            total_introns = len(introns_of_interest.split("_"))
            n_introns_gene = int(hg38_intron_df.loc[hg38_intron_df['gene_name']==gene_name]['intron_total'].drop_duplicates())
            strand = "".join(hg38_intron_df.loc[hg38_intron_df['gene_name']==gene_name]['strand'].drop_duplicates().tolist())
        
            results_df = get_score_per_path_pair_v2(df, total_introns, introns_of_interest, n_introns_gene, strand)
            results_df['gene_name'] = gene_name
            results_df['target'] = sample_name
        
            try:
                final_df = pd.concat([final_df,results_df]).reset_index(drop=True)
            except NameError:
                final_df = results_df
        
    return(final_df)

#############################

# Apply functions
genes_regions_dict = {'IFRD2':'7_6_5_4_3_2'}
sample_list = all_multi_introns_df['target'].drop_duplicates().tolist()
HeLa_test_paths = get_read_count_path_per_gene(all_multi_introns_df, genes_regions_dict, sample_list, hg38_intron_info_sub)


HeLa_test_paths.to_csv(out_splicing_paths, sep="\t", header=True, index=False)
