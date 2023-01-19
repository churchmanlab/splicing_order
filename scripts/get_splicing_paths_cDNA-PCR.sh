#!/bin/bash

# Author: Karine Choquet

# This script will call python scripts that will compute splicing order paths from cDNA-PCR sequencing data



###################

# CONFIG
baseDir="path_to_project_directory"
outDir="$baseDir/analysis/combined"

# Input files
ctrl_1="baseDir/analysis/sample1/sample1_hg38_multi_introns_counts.txt"
ctrl_2="baseDir/analysis/sample2/sample2_hg38_multi_introns_counts.txt"
U2_1="baseDir/analysis/sample3/sample3_hg38_multi_introns_counts.txt"
U2_2="baseDir/analysis/sample4/sample4_hg38_multi_introns_counts.txt"

# Output files
splicing_paths_with_SE="$outDir/all_samples_splicing_paths.with_SE.txt"
splicing_paths_no_SE="$outDir/all_samples_splicing_paths.no_SE.txt"

# Annotation files
annDir="path_to_annotations_directory"
gene_names_df="$annDir/hg38_UCSC_refGene_names.txt"
intron_info_df="$annDir/hg38_all_intron_features.txt"

# Python scripts
script_splice_paths_with_SE="get_splicing_paths.cDNA-PCR.with_SE.py"
script_splice_paths_no_SE="get_splicing_paths.cDNA-PCR.no_SE.py"

######################

module load samtools/1.9
module load bedtools/2.27.1
module load python/3.7.4
source my_virtual_env

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi

python $script_splice_paths_with_SE $gene_names_df $intron_info_df $ctrl_1 $ctrl_2 $U2_1 $U2_2 $splicing_paths_with_SE
python $script_splice_paths_no_SE $gene_names_df $intron_info_df $ctrl_1 $ctrl_2 $U2_1 $U2_2 $splicing_paths_no_SE
