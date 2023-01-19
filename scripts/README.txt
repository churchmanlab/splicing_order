The scripts in this folder were used to analyze the nanopore sequencing data in the manuscript "Pre-mRNA splicing order is predetermined and maintains splicing fidelity across multi-intronic transcripts".

To compute splicing order from chromatin-associated polyA+ direct RNA sequencing, use the following scripts after aligning reads to the reference genome:
1) get_splice_status_per_intron.directRNA.sh, which calls the scripts get_splice_status_per_intron.directRNA.py and get_multi_introns_isoforms_df.directRNA.py
2) splicing_order_paths_direct_RNA.sh, which calls the script splicing_order_paths_direct_RNA.py

To compute splicing status for intron/exon pairs from total polyA+ direct RNA sequencing, use the following scripts after aligning reads to the reference genome:
1) get_splice_status_per_intron.directRNA.sh, which calls the scripts get_splice_status_per_intron.directRNA.py and get_multi_introns_isoforms_df.directRNA.py
2) get_splice_status_exons_direct_RNA_nanopore_seq.sh, which calls the scripts get_splice_status_per_exon.directRNA.py and get_multi_exons_isoforms_df.directRNA.py
3) get_splice_status_per_intron_exon_pair_directRNA.py

To compute splicing order from chromatin-associated cDNA-PCR nanopore sequencing, use the following scripts after aligning reads to the reference genome:
1) get_splice_status_introns_cDNA_PCR_nanopore_seq.sh, which calls the script get_splice_status_per_intron.PCR.py
2) get_splicing_paths_cDNA-PCR.sh, which calls the scripts get_splicing_paths.cDNA-PCR.no_SE.py and get_splicing_paths.cDNA-PCR.with_SE.py

To reproduce the analyses and plots presented in the figures, use the Jupyter notebooks and R markdown files that start with the corresponding Figure names. For some figures, a Jupyter notebook and an R markdown file were used to generate different parts of the figures.
