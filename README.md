# sgRNApoolanalyzer
sgRNA Pooled Crispr Target Analyzer

R script used to quickly detect pooled crispr targets of set lengths (in our case, 20) and output density and distribution plots.  This script checks for exact matches but is being modified to account for slight variations within sequence.  R script will most likely be replaced with python with small calls to R script with subprocess module. Takes in single BAM and BAI file, reference file containing target name and target sequence, and user inputted barcodes to identify region to be aligned to sequence within reference file.  If no barcodes are available for the target, then the barcodes used during sequencing can be used as well.   

Density plot.R - produces density plot curve to display distribution of density of correctly aligned reads.  Tighter curves means a more even distribution of reads for all correctly aligned sequences to reference file.  

NGS Exact vs Misaligned CrisprID Counter.R - Take counts of miscalls and exact matches to sequence by Target name specified in reference file

NGS Exact vs Misaligned Seuqnece Counter.R - Take counts of miscalls and exact matches

NGS CRISPR and Sequence Exact and Misaligned comparison.R - takes counts from both scripts and aligns both by crispr ID and sequence to output excel that shows ratios of miscalls vs exact alignments.  

NGS Exact vs Misaligned pipeline.R - File that runs other scrits in order.
