"""

Author: Karine Choquet

Date: January 29, 2022

This script will identify possible splicing paths and assign them a score based on the isoform counts from direct chromatin RNA sequencing

Usage: python splicing_order_paths_direct_RNA.sh gene_names_df intron_df rep1_multi_intron_counts rep2_multi_intron_counts rep1_splice_info rep2_splice_info splicing_paths

gene_names_df is an annotation file containing transcript ids (NM_) in one column and the gene name in the other
intron_df is a BED file with coordinates of each intron of interest
*_multi_intron_counts and *_splice_info are output files from get_splice_status_introns_direct_RNA_nanopore_seq.sh containing the splicing statuses of each intron and each read 


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

from interruptingcow import timeout

#############

gene_names_df = pd.read_table(sys.argv[1])
gene_names_df.columns = ['gene_name','gene_id']
intron_df = pd.read_table(sys.argv[2], header=None, names=['chrom','start','end','intron_name','score','strand'], dtype={'chrom':str, 'start':int, 'end':int, 'intron_name':str, 'score':str, 'strand':str})
intron_df['gene_id'] = intron_df['intron_name'].str.split("\\.").str[0]
intron_total_df = pd.DataFrame(intron_df.groupby('gene_id')['intron_name'].count()).rename(columns={'intron_name':'intron_total'})
gene_names_df = gene_names_df.merge(intron_total_df, on='gene_id')

# Set thresholds
# Minimum reads to consider introns to be post-transcriptionally spliced
min_reads = 10
level_threshold = 10

# Load multi_intron_df from both replicates
multi_intron_counts1 = pd.read_table(sys.argv[3])
multi_intron_counts1['gene_id'] = multi_intron_counts1['gene'].str.split("\\.").str[0]
multi_intron_counts1 = multi_intron_counts1.merge(gene_names_df, on=['gene_id'])

multi_intron_counts2 = pd.read_table(sys.argv[4])
multi_intron_counts2['gene_id'] = multi_intron_counts2['gene'].str.split("\\.").str[0]
multi_intron_counts2 = multi_intron_counts2.merge(gene_names_df, on=['gene_id'])

# Load splicing_info from both replicates
chr1_splice_info = pd.read_table(sys.argv[5], dtype={'read_name':str,'chrom':str,'intron_start':int,'intron_end':int,'strand':str,'gene_name':str,'intron_count':int,'read_overlap':int,'splice_status':str})
chr2_splice_info = pd.read_table(sys.argv[6], dtype={'read_name':str,'chrom':str,'intron_start':int,'intron_end':int,'strand':str,'gene_name':str,'intron_count':int,'read_overlap':int,'splice_status':str})


# Identify unspliced introns
chr1_splice_info_NO = chr1_splice_info[chr1_splice_info['splice_status']=='NO'].reset_index(drop=True)
chr1_gr_per_intron_NO = pd.DataFrame(chr1_splice_info_NO.groupby(['chrom','intron_start','intron_end','strand','gene_name','intron_count'])['read_name'].count()).reset_index().rename(columns={'read_name':'read_count'})
chr1_gr_per_intron_NO = chr1_gr_per_intron_NO[chr1_gr_per_intron_NO['read_count']>=min_reads].reset_index(drop=True).rename(columns={'gene_name':'gene'})

chr2_splice_info_NO = chr2_splice_info[chr2_splice_info['splice_status']=='NO'].reset_index(drop=True)
chr2_gr_per_intron_NO = pd.DataFrame(chr2_splice_info_NO.groupby(['chrom','intron_start','intron_end','strand','gene_name','intron_count'])['read_name'].count()).reset_index().rename(columns={'read_name':'read_count'})
chr2_gr_per_intron_NO = chr2_gr_per_intron_NO[chr2_gr_per_intron_NO['read_count']>=min_reads].reset_index(drop=True).rename(columns={'gene_name':'gene'})


# Merge the introns from the two replicates
gr_per_intron_NO = chr1_gr_per_intron_NO.merge(chr2_gr_per_intron_NO, on=['chrom','intron_start','intron_end','strand','gene','intron_count'])


# Identify transcripts that have at least 3 introns that meet this minimum
chr_transcripts = pd.DataFrame(gr_per_intron_NO.groupby(['gene','strand'])['intron_count'].count()).reset_index().rename(columns={'intron_count':'n_introns'})
chr_transcripts = chr_transcripts[chr_transcripts['n_introns']>=3].reset_index(drop=True)

# Merge with gene names and total intron counts
chr_transcripts['gene_id'] = chr_transcripts['gene'].str.split("\\.").str[0]
chr_transcripts = chr_transcripts.merge(gene_names_df, on='gene_id')

# Merge transcripts and introns
gr_per_intron_NO = gr_per_intron_NO.merge(chr_transcripts, on=['gene'])

# Make a dictionary that contains the identity of the introns to analyze for each transcript
# Alternatively, a list of intron groups (e.g. those with AS events) can be provided

transcript_dict = {}

for i in range(len(gr_per_intron_NO)):
    gene = gr_per_intron_NO.loc[i]['gene']
    intron_count = gr_per_intron_NO.loc[i]['intron_count']
    
    if gene not in transcript_dict.keys():
        transcript_dict[gene] = [intron_count]
    elif gene in transcript_dict.keys():
        transcript_dict[gene].append(intron_count)


####################

def get_score_per_path_non_consec(multi_introns_df, total_introns_of_interest, introns_of_interest_list_tmp, n_introns_gene, strand):

    # Make a dictionary with the patterns and the counts and another with the number of introns spliced and the counts
    pattern_dict = {}
    results_list = []
    spliced_counts_dict = {}
    
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
    #level_threshold = 10
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
    
            final_final_df['analyzed_introns'] = introns_of_interest_fix
            final_final_df['gene_name'] = gene_name
               
            return(final_final_df)
        
        except ValueError:
            pass

################################

results_df_list = []

s1 = 'chr_rep1'
s2 = 'chr_rep2'

for transcript in list(transcript_dict.keys()):
    
    gene_results_list = []
    print(transcript)
    
    # Retrieve necessary information for that transcript
    multi_introns_df1 = multi_intron_counts1[multi_intron_counts1['gene']==transcript].reset_index(drop=True)
    multi_introns_df2 = multi_intron_counts2[multi_intron_counts2['gene']==transcript].reset_index(drop=True)
    n_introns_gene = int(chr_transcripts[chr_transcripts['gene']==transcript]['intron_total'])
    strand = chr_transcripts[chr_transcripts['gene']==transcript]['strand'].tolist()[0]
    introns_of_interest = transcript_dict[transcript]
    total_introns_of_interest = len(introns_of_interest)
    
    # Iterate through groups of 3-4 introns in both replicates to retrieve the longest possible paths
    for i in range(3,5):
        for x in range(total_introns_of_interest-i+1):
            introns_of_interest_sub = [introns_of_interest[a] for a in range(x,x+i)]
            # Compute splicing order paths for the introns of interest
            results_df1 = get_score_per_path_non_consec(multi_introns_df1, i, introns_of_interest_sub, n_introns_gene, strand)
            results_df2 = get_score_per_path_non_consec(multi_introns_df2, i, introns_of_interest_sub, n_introns_gene, strand)

    
            if results_df1 is not None and results_df2 is not None:
                results_df1['sample_name'] = s1
                results_df2['sample_name'] = s2
                results_df = pd.concat([results_df1,results_df2]).reset_index(drop=True)
                results_df['gene'] = transcript
                results_df['n_analyzed_introns'] = i
                gene_results_list.append(results_df)
                
    if len(gene_results_list)>0:
        gene_results_df = pd.concat(gene_results_list).reset_index(drop=True)
    
        # Retrieve the intron groups that were successfully analyzed
        analyzed_introns_list = gene_results_df['analyzed_introns'].drop_duplicates().tolist()
    
        # Return the longest possible paths
        good_groups = []
        subgroup_list = []
    
    
        # Sort intron list from longest to shortest
        sorted_intron_list = sorted(analyzed_introns_list, key=len, reverse=True)
    
        for intron_group in sorted_intron_list:
            if len(good_groups) == 0:
                good_groups.append(intron_group)
            else:
                common_count = 0
                for a in good_groups:
                    if intron_group in a:
                        common_count += 1
                if common_count == 0 and intron_group not in good_groups:
                    good_groups.append(intron_group)                   
 
        # Retrieve those paths:
        gene_results_df_sub = gene_results_df[gene_results_df['analyzed_introns'].isin(good_groups)].reset_index(drop=True)
        results_df_list.append(gene_results_df_sub)
    

final_results_df = pd.concat(results_df_list).reset_index(drop=True)

final_results_df.to_csv(sys.argv[7], sep="\t", header=True, index=False)

