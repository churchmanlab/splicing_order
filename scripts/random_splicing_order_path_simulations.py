"""

Author: Karine Choquet

Date: March 14, 2022

This script will randomly splice reads from the intron groups used in splicing order analyses and compute splicing order paths from these simulations

Usage: python random_splicing_order_path_simulations.py $intron_df $gene_names_df $splicing_paths $out

gene_names_df is an annotation file containing transcript ids (NM_) in one column and the gene name in the other
intron_df is a BED file with coordinates of each intron of interest
splicing_paths is the output from K562_splicing_order_paths.sh 


"""


import numpy as np
import pandas as pd

import itertools

import random
from numpy.random import choice

import sys

from collections import Counter


# Load intron features and gene_names df
hg38_intron_df = pd.read_table(sys.argv[1])
gene_names_df = pd.read_table(sys.argv[2], header=None)
gene_names_df.columns = ['gene_name','gene_id']

hg38_intron_df['gene_id'] = hg38_intron_df['gene'].str.split("\\.").str[0]
hg38_intron_df = hg38_intron_df.merge(gene_names_df, on='gene_id')
hg38_intron_coord = hg38_intron_df.copy()[['chrom','start','end','gene','intron_pos']]

# Load results obtained from identifying splicing order paths in both replicates
path_df_tmp = pd.read_table(sys.argv[3])
path_df = path_df_tmp[['sample_name','gene','gene_name','analyzed_introns','n_analyzed_introns','full_path','full_path_score','rank']].drop_duplicates().sort_values(by=['gene_name','analyzed_introns','rank']).reset_index(drop=True)

# Remove duplicates intron groups because they belong to different transcripts
intron_groups = path_df[['gene','gene_name','analyzed_introns','n_analyzed_introns']].drop_duplicates().reset_index(drop=True)

intron_groups_3 = intron_groups[intron_groups['n_analyzed_introns']==3].reset_index(drop=True)
intron_groups_3['int1'] = intron_groups_3['analyzed_introns'].str.split("_").str[0].astype(int)
intron_groups_3['int2'] = intron_groups_3['analyzed_introns'].str.split("_").str[1].astype(int)
intron_groups_3['int3'] = intron_groups_3['analyzed_introns'].str.split("_").str[2].astype(int)

intron_groups_3 = intron_groups_3.merge(hg38_intron_coord, left_on=['gene','int1'], right_on=['gene','intron_pos']).rename(columns={'chrom':'chrom_1', 'start':'start_1', 'end':'end_1'})
intron_groups_3 = intron_groups_3.merge(hg38_intron_coord, left_on=['gene','int2'], right_on=['gene','intron_pos']).rename(columns={'chrom':'chrom_2', 'start':'start_2', 'end':'end_2'})
intron_groups_3 = intron_groups_3.merge(hg38_intron_coord, left_on=['gene','int3'], right_on=['gene','intron_pos']).rename(columns={'chrom':'chrom_3', 'start':'start_3', 'end':'end_3'})

intron_groups_3 = intron_groups_3.sort_values(by=['gene','int1','int2','int3']).drop_duplicates(subset=['chrom_1','start_1','end_1','chrom_2','start_2','end_2','chrom_3','start_3','end_3']).reset_index(drop=True)

intron_groups_4 = intron_groups[intron_groups['n_analyzed_introns']==4].reset_index(drop=True)
intron_groups_4['int1'] = intron_groups_4['analyzed_introns'].str.split("_").str[0].astype(int)
intron_groups_4['int2'] = intron_groups_4['analyzed_introns'].str.split("_").str[1].astype(int)
intron_groups_4['int3'] = intron_groups_4['analyzed_introns'].str.split("_").str[2].astype(int)
intron_groups_4['int4'] = intron_groups_4['analyzed_introns'].str.split("_").str[3].astype(int)

intron_groups_4 = intron_groups_4.merge(hg38_intron_coord, left_on=['gene','int1'], right_on=['gene','intron_pos']).rename(columns={'chrom':'chrom_1', 'start':'start_1', 'end':'end_1'})
intron_groups_4 = intron_groups_4.merge(hg38_intron_coord, left_on=['gene','int2'], right_on=['gene','intron_pos']).rename(columns={'chrom':'chrom_2', 'start':'start_2', 'end':'end_2'})
intron_groups_4 = intron_groups_4.merge(hg38_intron_coord, left_on=['gene','int3'], right_on=['gene','intron_pos']).rename(columns={'chrom':'chrom_3', 'start':'start_3', 'end':'end_3'})
intron_groups_4 = intron_groups_4.merge(hg38_intron_coord, left_on=['gene','int4'], right_on=['gene','intron_pos']).rename(columns={'chrom':'chrom_4', 'start':'start_4', 'end':'end_4'})

intron_groups_4 = intron_groups_4.sort_values(by=['gene','int1','int2','int3']).drop_duplicates(subset=['chrom_1','start_1','end_1','chrom_2','start_2','end_2','chrom_3','start_3','end_3','chrom_4','start_4','end_4']).reset_index(drop=True)

fields = ['gene','gene_name','analyzed_introns','n_analyzed_introns']
intron_groups_nodup = pd.concat([intron_groups_3[fields],intron_groups_4[fields]]).sort_values(by=['gene_name','analyzed_introns']).reset_index(drop=True)

# Merge back with paths
path_df_tmp_nodup = path_df_tmp.merge(intron_groups_nodup, on=['gene','gene_name','analyzed_introns','n_analyzed_introns'])
path_df_nodup = path_df_tmp_nodup[['sample_name','gene','gene_name','analyzed_introns','n_analyzed_introns','full_path','full_path_score','rank']].drop_duplicates().sort_values(by=['gene_name','analyzed_introns','rank']).reset_index(drop=True)


# Reformat the table to have replicates side by side
path_df_nodup_rep1 = path_df_nodup[path_df_nodup['sample_name']=='chr_rep1'].reset_index(drop=True)
path_df_nodup_rep2 = path_df_nodup[path_df_nodup['sample_name']=='chr_rep2'].reset_index(drop=True)

path_df_piv = path_df_nodup_rep1.merge(path_df_nodup_rep2, on=['gene_name','gene','analyzed_introns','n_analyzed_introns','full_path'])


# Retrieve the paths that are reproducible (same rank in both samples)
path_df_reprod = path_df_piv[(path_df_piv['rank_x']==path_df_piv['rank_y'])].reset_index(drop=True)

# Convert back into a long format dataframe
path_df_reprod_rep1 = path_df_reprod[['sample_name_x','gene','gene_name','analyzed_introns','n_analyzed_introns','full_path','full_path_score_x','rank_x']].rename(columns={'sample_name_x':'sample_name','full_path_score_x':'full_path_score','rank_x':'rank'})
path_df_reprod_rep2 = path_df_reprod[['sample_name_y','gene','gene_name','analyzed_introns','n_analyzed_introns','full_path','full_path_score_y','rank_y']].rename(columns={'sample_name_y':'sample_name','full_path_score_y':'full_path_score','rank_y':'rank'})


path_df_reprod_m = pd.concat([path_df_reprod_rep1,path_df_reprod_rep2]).reset_index(drop=True)



def get_score_per_path_non_consec(multi_introns_df, total_introns_of_interest, introns_of_interest_list):

    # Make a dictionary with the patterns and the counts and another with the number of introns spliced and the counts
    pattern_dict = {}
    results_list = []
    spliced_counts_dict = {}
    
    introns_of_interest = "_".join(introns_of_interest_list)
    
    # Initiate pattern_dict
    pattern_dict = {}
    for n in range(len(introns_of_interest_list)+1):
        pattern_dict[n] = {}
    
    # Iterate to get counts for isoforms
    for row in range(len(multi_introns_df)):
        gene_name = multi_introns_df.loc[row]['gene_name']
        splice_status_temp = multi_introns_df.loc[row]['splice_status']
        intron_numbers_list = multi_introns_df.loc[row]['intron_numbers'].split("_")
            
        # Determine if all introns of interest are present in the splice status
        common_introns = [a for a in introns_of_interest_list if a in intron_numbers_list]
        
        if len(common_introns) == total_introns_of_interest:
            introns_of_interest_pos = [i for i, x in enumerate(intron_numbers_list) if x in introns_of_interest_list]
            splice_status_list_temp1 = splice_status_temp.split("_")
            splice_status_list_temp = [splice_status_list_temp1[a] for a in introns_of_interest_pos]
            splice_status_list = ["SKP" if "SKP" in a else a for a in splice_status_list_temp]
            splice_status = "_".join(splice_status_list)
            pattern_count = multi_introns_df.loc[row]['count']
            skipped_count = Counter(splice_status_list)['SKP']
            spliced_count = Counter(splice_status_list)['YES']
            unspliced_count = Counter(splice_status_list)['NO']
            undetermined_count = Counter(splice_status_list)['UND']
            
            if skipped_count == 0 and undetermined_count == 0: # no skipped or undetermined introns
                level = spliced_count
                if skipped_count < total_introns_of_interest:
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

    
    # Filter for a certain number of reads at each intermediate level
    level_threshold = 1
    good_levels = []
    for level in sorted(list(spliced_counts_dict.keys()))[0:-1]: # exclude the final level (all spliced)
        if spliced_counts_dict[level] >= level_threshold:
            good_levels.append(level)
            
    if len(good_levels) == total_introns_of_interest: # there is at least one read at each level
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
    
        
        # Now loop through the dictionary to match each pair with the possible final paths
        # The final level contains only final paths, so first retrieve those
        try:
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
            last_isos = final_final_df[~final_final_df['path_name'].str.endswith("->")][['full_path','path_score']].drop_duplicates().reset_index(drop=True)
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
    
            final_final_df['analyzed_introns'] = introns_of_interest
            final_final_df['gene_name'] = gene_name
               
            return(final_final_df)
        
        except ValueError:
            pass


def random_splicing(df):
    
    random_reads = []
    
    for i in range(len(df)):
        sample_name = df.loc[i]['sample_name']
        gene = df.loc[i]['gene']
        gene_name = df.loc[i]['gene_name']
        analyzed_introns = df.loc[i]['analyzed_introns']
        level = df.loc[i]['level']
        count = df.loc[i]['count_pattern2']
    
        intron_list = analyzed_introns.split("_")
        for i in range(count):
            spliced_introns = random.sample(intron_list, level)
            read_splice_status = []
            for a in intron_list:
                if a not in spliced_introns:
                    read_splice_status.append("NO")
                elif a in spliced_introns:
                    read_splice_status.append("YES")
            read_splice_status_join = "_".join(read_splice_status)
            random_reads.append([sample_name, gene, gene_name, analyzed_introns, level, read_splice_status_join])
        
    random_reads_df = pd.DataFrame(random_reads)
    random_reads_df.columns = ['sample_name','gene','gene_name','analyzed_introns', 'level', 'read_splice_status_join']
    
    return random_reads_df



def compare_splicing_paths_to_random(gene, gene_name, analyzed_introns, n_iterations):
    
    random_paths_list = []
    
    # Retrieve the patterns corresponding to that gene and analyzed introns
    pattern_counts_by_level_sub = pattern_counts_by_level_all[(pattern_counts_by_level['gene']==gene) & (pattern_counts_by_level['analyzed_introns']==analyzed_introns)].reset_index(drop=True)
    
    # Compute random splicing order paths many times
    for i in range(1,n_iterations+1):
        
        random_reads_df = random_splicing(pattern_counts_by_level_sub)
        random_reads_counts = pd.DataFrame(random_reads_df.reset_index().groupby(['sample_name','gene','gene_name','analyzed_introns','level','read_splice_status_join']).count()).reset_index()
        random_reads_counts.columns = ['sample_name','gene','gene_name','intron_numbers','level','splice_status','count']
    
        analyzed_introns_list = analyzed_introns.split("_")
        n_introns = len(analyzed_introns_list)
    
        random_paths = get_score_per_path_non_consec(random_reads_counts, n_introns, analyzed_introns_list)
        random_paths_clean = random_paths[['full_path','full_path_score','rank']].drop_duplicates().sort_values(by='rank')
        random_paths_clean['iteration'] = i
        random_paths_clean['gene'] = gene
        random_paths_clean['gene_name'] = gene_name
        random_paths_clean['analyzed_introns'] = analyzed_introns
    
        random_paths_list.append(random_paths_clean)
        random_paths_df = pd.concat(random_paths_list)
    
    return(random_paths_df)


def random_and_SS_splicing(df, hg38_intron_df, gene, analyzed_introns):
    
    # Retrieve the splice site scores for the introns of interest
    intron_list = [int(a) for a in analyzed_introns.split("_")]
    hg38_test_df = hg38_intron_df[hg38_intron_df['gene']==gene][['intron_pos','MaxEnt_score_5SS','MaxEnt_score_3SS']]
    hg38_test_df.index = hg38_test_df['intron_pos']
    hg38_test_df = hg38_test_df.drop(columns=['intron_pos'])
    hg38_test_dict = hg38_test_df.to_dict()
        
    SS5 = [hg38_test_dict['MaxEnt_score_5SS'][a] for a in intron_list]
    SS3 = [hg38_test_dict['MaxEnt_score_3SS'][a] for a in intron_list]
        
    # Calculate sum of splice site scores
    SS_sum = [2**((SS5[a] + SS3[a])/2) for a in range(len(SS5))]
    SS_sum_bis = [a if a > 0 else 0.001 for a in SS_sum]
    SS_prob = [a/np.sum(SS_sum_bis) for a in SS_sum_bis]
    
    
    random_reads = []
    
    for i in range(len(df)):
        sample_name = df.loc[i]['sample_name']
        gene = df.loc[i]['gene']
        gene_name = df.loc[i]['gene_name']
        analyzed_introns = df.loc[i]['analyzed_introns']
        level = df.loc[i]['level']
        count = df.loc[i]['count_pattern2']
        
        for i in range(count):
            spliced_introns = choice(intron_list, level, p=SS_prob, replace=False)
            read_splice_status = []
            for a in intron_list:
                if a not in spliced_introns:
                    read_splice_status.append("NO")
                elif a in spliced_introns:
                    read_splice_status.append("YES")
            read_splice_status_join = "_".join(read_splice_status)
            random_reads.append([sample_name, gene, gene_name, analyzed_introns, level, read_splice_status_join])
        
    random_reads_df = pd.DataFrame(random_reads)
    random_reads_df.columns = ['sample_name','gene','gene_name','analyzed_introns', 'level', 'read_splice_status_join']
    
    return random_reads_df


def compare_splicing_paths_to_random_SS(gene, gene_name, analyzed_introns, n_iterations):
    
    random_paths_list = []
    
    # Retrieve the patterns corresponding to that gene and analyzed introns
    pattern_counts_by_level_sub = pattern_counts_by_level_all[(pattern_counts_by_level['gene']==gene) & (pattern_counts_by_level['analyzed_introns']==analyzed_introns)].reset_index(drop=True)
    
    # Compute random splicing order paths many times
    for i in range(1,n_iterations+1):
        
        random_reads_df = random_and_SS_splicing(pattern_counts_by_level_sub, hg38_intron_df, gene, analyzed_introns)
        random_reads_counts = pd.DataFrame(random_reads_df.reset_index().groupby(['sample_name','gene','gene_name','analyzed_introns','level','read_splice_status_join']).count()).reset_index()
        random_reads_counts.columns = ['sample_name','gene','gene_name','intron_numbers','level','splice_status','count']
    
        analyzed_introns_list = analyzed_introns.split("_")
        n_introns = len(analyzed_introns_list)
    
        random_paths = get_score_per_path_non_consec(random_reads_counts, n_introns, analyzed_introns_list)
        random_paths_clean = random_paths[['full_path','full_path_score','rank']].drop_duplicates().sort_values(by='rank')
        random_paths_clean['iteration'] = i
        random_paths_clean['gene'] = gene
        random_paths_clean['gene_name'] = gene_name
        random_paths_clean['analyzed_introns'] = analyzed_introns
    
        random_paths_list.append(random_paths_clean)
        random_paths_df = pd.concat(random_paths_list)
    
    return(random_paths_df)


def random_and_seq_splicing(df, hg38_intron_df, gene, analyzed_introns):
    
    # Define a score based on the position of the intron
    intron_list = [int(a) for a in analyzed_introns.split("_")]
    #pos_score = [1/(a+1) for a,x in enumerate(intron_list)]
    #pos_score_norm = [a/sum(pos_score) for a in pos_score]
    
    if len(intron_list) == 4:
        if intron_list[0] < intron_list[-1]:
            pos_score = [0.75,0.15,0.08,0.02]
        elif intron_list[0] > intron_list[-1]:
            pos_score = [0.02,0.08,0.15,0.75]
    elif len(intron_list) == 3:
        if intron_list[0] < intron_list[-1]:
            pos_score = [0.75,0.23,0.02]
        elif intron_list[0] > intron_list[-1]:
            pos_score = [0.02,0.23,0.75]

    
    random_reads = []
    
    for i in range(len(df)):
        sample_name = df.loc[i]['sample_name']
        gene = df.loc[i]['gene']
        gene_name = df.loc[i]['gene_name']
        analyzed_introns = df.loc[i]['analyzed_introns']
        level = df.loc[i]['level']
        count = df.loc[i]['count_pattern2']
        
        for i in range(count):
            spliced_introns = choice(intron_list, level, p=pos_score, replace=False)
            read_splice_status = []
            for a in intron_list:
                if a not in spliced_introns:
                    read_splice_status.append("NO")
                elif a in spliced_introns:
                    read_splice_status.append("YES")
            read_splice_status_join = "_".join(read_splice_status)
            random_reads.append([sample_name, gene, gene_name, analyzed_introns, level, read_splice_status_join])
        
    random_reads_df = pd.DataFrame(random_reads)
    random_reads_df.columns = ['sample_name','gene','gene_name','analyzed_introns', 'level', 'read_splice_status_join']
    
    return random_reads_df


def compare_splicing_paths_to_random_sequential(gene, gene_name, analyzed_introns, n_iterations):
    
    random_paths_list = []
    
    # Retrieve the patterns corresponding to that gene and analyzed introns
    pattern_counts_by_level_sub = pattern_counts_by_level_all[(pattern_counts_by_level['gene']==gene) & (pattern_counts_by_level['analyzed_introns']==analyzed_introns)].reset_index(drop=True)
    
    # Compute random splicing order paths many times
    for i in range(1,n_iterations+1):
        
        random_reads_df = random_and_seq_splicing(pattern_counts_by_level_sub, hg38_intron_df, gene, analyzed_introns)
        random_reads_counts = pd.DataFrame(random_reads_df.reset_index().groupby(['sample_name','gene','gene_name','analyzed_introns','level','read_splice_status_join']).count()).reset_index()
        random_reads_counts.columns = ['sample_name','gene','gene_name','intron_numbers','level','splice_status','count']
    
        analyzed_introns_list = analyzed_introns.split("_")
        n_introns = len(analyzed_introns_list)
    
        random_paths = get_score_per_path_non_consec(random_reads_counts, n_introns, analyzed_introns_list)
        random_paths_clean = random_paths[['full_path','full_path_score','rank']].drop_duplicates().sort_values(by='rank')
        random_paths_clean['iteration'] = i
        random_paths_clean['gene'] = gene
        random_paths_clean['gene_name'] = gene_name
        random_paths_clean['analyzed_introns'] = analyzed_introns
    
        random_paths_list.append(random_paths_clean)
        random_paths_df = pd.concat(random_paths_list)
    
    return(random_paths_df)


# To simulate reads, retrieve the number of reads per intron group for each level
pattern_counts = path_df_tmp_nodup[['sample_name','gene','gene_name','analyzed_introns','pattern2','count_pattern2','level']].drop_duplicates()
pattern_counts_by_level = pd.DataFrame(pattern_counts.groupby(['sample_name','gene','gene_name','analyzed_introns','level'])['count_pattern2'].sum()).reset_index()

pattern_counts_by_level_all = pd.DataFrame(pattern_counts_by_level.groupby(['sample_name','gene','gene_name','analyzed_introns','level'])['count_pattern2'].sum()).reset_index()




# Analyze all genes
cand_gene_df = path_df_reprod[['gene','gene_name','analyzed_introns']].drop_duplicates().reset_index(drop=True)

number_iter = 100

results_list = []

for i in range(len(cand_gene_df)):
    gene = cand_gene_df.loc[i]['gene']
    gene_name = cand_gene_df.loc[i]['gene_name']
    analyzed_introns = cand_gene_df.loc[i]['analyzed_introns']
    
    true_paths = path_df[(path_df['gene']==gene) & (path_df['analyzed_introns']==analyzed_introns) & (path_df['sample_name']=='chr_rep1')].reset_index(drop=True)
    random_paths = compare_splicing_paths_to_random(gene, gene_name, analyzed_introns, number_iter)
    seq_paths = compare_splicing_paths_to_random_sequential(gene, gene_name, analyzed_introns, number_iter)
    SS_paths = compare_splicing_paths_to_random_SS(gene, gene_name, analyzed_introns, number_iter)
    
    rand_df = pd.DataFrame(random_paths.groupby(['full_path','gene','gene_name','analyzed_introns'])['full_path_score'].mean()).reset_index() #.merge(true_paths, on=['full_path','gene','analyzed_introns'], how='left').fillna(0)
    seq_df = pd.DataFrame(seq_paths.groupby(['full_path','gene','gene_name','analyzed_introns'])['full_path_score'].mean()).reset_index() #.merge(true_paths, on=['full_path','gene','analyzed_introns'], how='left').fillna(0)
    SS_df = pd.DataFrame(SS_paths.groupby(['full_path','gene','gene_name','analyzed_introns'])['full_path_score'].mean()).reset_index() #.merge(true_paths, on=['full_path','gene','analyzed_introns'], how='left').fillna(0)
    
    # Add rank
    rand_df = rand_df.sort_values(by='full_path_score',ascending=False).reset_index(drop=True)
    rand_df['rank'] = rand_df.index + 1
    
    seq_df = seq_df.sort_values(by='full_path_score',ascending=False).reset_index(drop=True)
    seq_df['rank'] = seq_df.index + 1

    SS_df = SS_df.sort_values(by='full_path_score',ascending=False).reset_index(drop=True)
    SS_df['rank'] = SS_df.index + 1
    
    # Add category
    rand_df['approach'] = 'random'
    seq_df['approach'] = 'sequential'
    SS_df['approach'] = 'SS_scores'
    
    order_df = pd.concat([rand_df, seq_df, SS_df]).reset_index(drop=True)
    
    results_list.append(order_df)
    
results_df = pd.concat(results_list)
    


results_df.to_csv(sys.argv[-1], sep="\t", header=True, index=False)
