#############################################################

This is the repository for the paper "Probing TCR specificity using artificial in vivo diversification of CDR3 regions"
by Orlando B. Giorgetti, Annette Haas-Assenbaum & Thomas Boehm

#############################################################

The R code will generate figures from repertoire data tables containing CDR3 counts.

If you require additional help, contact OBG at orlandogiorgetti@gmail.com.

#############################################################



Raw data
========

In all cases raw data consists of Illumina sequencing of mouse TCRa and TCRb cDNA amplicons with UMI barcoding (PRJNA1128587,see paper).
Sequencing data comes from either TCR alpha- or TCR beta-CDR3 edited mice. Note there is data from both genes for each mouse, but only
the data from the pertinent edited gene was analyzed. The table of the not-edited gene consists of the SMARTA receptor sequence for that chain.

I provided examples of how to obtain the graph statistics for network plots in the code comments. The plots should match the published data, 
with the exception of random placement of network nodes and scatterplot dots, which do not affect interpretation of the data.

CDR3 counts and metadata
========================

The metadata (mainly FACS sort population percentages) is also available:

"XP.th.xlsx" # Metadata - Table with data from sorted mice Excel table with data from sorted mice - not an excel file!  
"XP.th.xlsx.2" # Metadata - Table with data from sorted mice Excel table with data from sorted mice - not an excel file!  
"MM.xp.3a" # Counts tables with the sequences from sort experiments (TCR alpha edited mice)  
"MM.xp.3b" # Counts tables with the sequences from sort experiments (TCR beta edited mice)  
"TCR.a.spp.MM.CDR3.amino.acid.sequence" # WT TCR alpha sequences from B6 controls (From Giorgetti et al. 2023)  
"TCR.b.spp.MM.CDR3.amino.acid.sequence" # WT TCR beta sequences from B6 controls (From Giorgetti et al. 2023)  
