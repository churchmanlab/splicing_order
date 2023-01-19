#!/bin/bash

# Author: Karine Choquet

# This script will perform the following steps:
# 1) Convert the sample BAM file to a BED file, including the CIGAR string
# 2) Intersect each read with introns of interest
# 3) Call a python script that determines the splicing status of each intron within the read
# 4) Call a python script that concatenates the splicing status of each intron per read

# Usage: ./get_splice_status_introns_direct_RNA_nanopore_seq.sh sample_name


###################

# CONFIG
sample=$1

baseDir="path_to_project_directory"
outDir="$baseDir/analysis/$sample"

# Input files
bam="$baseDir/alignment/$sample/minimap2/${sample}_minimap2_uniq_sort.bam"
bed="$baseDir/alignment/$sample/minimap2/${sample}_minimap2_uniq_sort.bed"


# Output files
splicing_info="$outDir/${sample}_hg38_splicing_info_per_intron.txt"
splicing_dict="$outDir/${sample}_hg38_splicing_dictionary"
multi_introns_df="$outDir/${sample}_hg38_multi_introns_isoforms_df"
multi_introns_counts="$outDir/${sample}_hg38_multi_introns_isoforms_counts"

# Annotation files
intron_bed_file="NCBI_RefSeq_hg38_introns_parsed_bedtool_for_intron_pairs.bed" # OR GENCODE V38

# Python scripts
script_splice_dict="get_splice_status_per_intron.directRNA.py"
script_multi_introns="get_multi_introns_isoforms_df.directRNA.py"

######################

module load samtools/1.9
module load bedtools/2.27.1
module load python/3.7.4
source /home/kc248/myenv_python3/bin/activate

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi



# Get bed file from BAM file
echo "Converting BAM file to BED..."
bedtools bamtobed -cigar -tag NM -i $bam > $bed

# Intersect with introns
echo "Intersecting with annotated introns..."
bedtools intersect -s -a $bed -b $intron_bed_file -wo > $outDir/$sample.intersect_introns.bed

# Make splicing info and splicing dictionary datasets
echo "Getting splicing info..."
python $script_splice_dict $outDir/$sample.intersect_introns.bed $splicing_info $splicing_dict
python $script_multi_introns $splicing_dict.npy $multi_introns_df.txt $multi_introns_counts.txt
