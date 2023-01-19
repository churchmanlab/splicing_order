#!/bin/bash

# Author: Karine Choquet

# This script will compute splicing order paths from direct RNA nanopore sequencing data, for biological duplicates


# Usage: ./splicing_order_paths_direct_RNA.sh sample_name_rep1 sample_name_rep2




###################

# CONFIG
sample_rep1=$1
sample_rep2=$2

baseDir="path_to_project_directory"
outDir="$baseDir/analysis"

# Input files, obtained from get_splice_status_introns_direct_RNA_nanopore_seq.sh
rep1_splicing_info="$outDir/${sample_rep1}_hg38_splicing_info_per_intron.txt"
rep1_multi_introns_counts="$outDir/${sample_rep1}_hg38_multi_introns_isoforms_counts.txt"
rep2_splicing_info="$outDir/${sample_rep2}_hg38_splicing_info_per_intron.txt"
rep2_multi_introns_counts="$outDir/${sample_rep2}_hg38_multi_introns_isoforms_counts.txt"

# Output files
splicing_paths="$outDir/K562_chromatin_splicing_paths_non_consecutive_introns.max4introns.txt"

# Annotation files
annDir="path_to_annotations"
gene_names_df="$annDir/hg38_UCSC_refGene_names.txt"
intron_df="$annDir/NCBI_RefSeq_hg38_introns_parsed.bed"

# Python scripts
script_splice_paths="splicing_order_paths_direct_RNA.py"

######################

module load samtools/1.9
module load bedtools/2.27.1
module load python/3.7.4
source my_virtual_environment


# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi

python $script_splice_paths $gene_names_df $intron_df $rep1_multi_intron_counts $rep2_multi_intron_counts $rep1_splice_info $rep2_splice_info $splicing_paths

