#!/bin/bash

# Author: Karine Choquet

# This script will perform the following steps:
# 1) Intersect each read with introns of interest
# 2) Call a python script that determines the splicing status of each intron within the read, one chromosome at a time
# 3) Call a python script that concatenates the splicing status of each intron per read

# Usage: ./get_splice_status_exons_direct_RNA_nanopore_seq.sh sample_name



###################

# CONFIG
sample=$1
baseDir="path_to_project_directory"
outDir="$baseDir/analysis/$sample"

# Input files
bed="$baseDir/alignment/$sample/minimap2/${sample}_minimap2_uniq_sort.bed"


# Output files
splicing_info="$outDir/${sample}_hg38_splicing_info_per_exon.txt"
splicing_dict="$outDir/${sample}_hg38_splicing_dictionary_exons"
multi_exons_df="$outDir/${sample}_hg38_multi_exons_isoforms_df"
multi_exons_counts="$outDir/${sample}_hg38_multi_exons_isoforms_counts"

# Annotation files
exon_bed_file="path_to_annotation_files/NCBI_RefSeq_hg38_exons_parsed.bed"

# Python scripts
script_splice_dict="get_splice_status_per_exon.directRNA.py"
script_multi_exons="get_multi_exons_isoforms_df.directRNA.py"

######################

module load samtools/1.9
module load bedtools/2.27.1
module load python/3.7.4
source /home/kc248/myenv_python3/bin/activate

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Intersect with exons
echo "Intersecting with annotated exons..."
bedtools intersect -s -a $bed -b $exon_bed_file -wo > $outDir/$sample.intersect_exons.bed

# Make splicing info and splicing dictionary datasets
echo "Getting splicing info..."
python $script_splice_dict $outDir/$sample.intersect_exons.bed $splicing_info $splicing_dict

# Make splicing info and splicing dictionary datasets
cd $outDir
echo -e "read\tgene\tstrand\texon_numbers\tsplice_status" > ${multi_exons_df}_all.txt
echo -e "gene\tstrand\texon_numbers\tsplice_status\tcount" > ${multi_exons_counts}_all.txt
for i in *exons_chr*.npy ; do
  chrom=$(echo $i | cut -f11- -d_ | cut -f1 -d.)
  python $script_multi_exons $i ${multi_exons_df}_${chrom}.txt ${multi_exons_counts}_${chrom}.txt
  cat ${multi_exons_df}_${chrom}.txt >> ${multi_exons_df}_all.txt
  cat ${multi_exons_counts}_${chrom}.txt >> ${multi_exons_counts}_all.txt
  rm ${multi_exons_df}_${chrom}.txt
  rm ${multi_exons_counts}_${chrom}.txt
done
