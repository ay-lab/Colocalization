# Colocalization
Co-localization analysis between reference GWAS and eQTL summary statistics

Sourya Bhattacharyya

La Jolla Institute for Immunology, La Jolla, CA 92037, USA

----------------------

Implements colocalization between the reference GWAS summary statistics and input eQTL file.

Performs colocalization using the R package coloc (https://chr1swallace.github.io/coloc/index.html)


Running the script
==================

First edit the configuration file (configfile.yaml). 

Check individual parameters related to both GWAS and eQTL summary statistics.

***Note: Parameters marked as mandatory need to be provided.

***Note: Allele frequency (AF) information need to be provided, either for the GWAS or for the eQTL data.

Once the configuration file is edited, run the script 

Colocalization_Analysis_GWAS_Script.sh

It will invoke the colocalization along with the corresponding configuration file.

Output
========

With respect to the specified output directory (parameter "OutDir"), check the following two files:

1. FINAL_Summary_Coloc_Gene_SNP_Pairs.bed

Lists colocalized variants along with the GWAS and eQTL statistics, one per GWAS loci.

(By GWAS loci, we mean +/- 500 Kb from the significant GWAS SNPs)

2. FINAL_Summary_Coloc_Gene_SNP_Pairs_95pct_credible_set.bed

For each GWAS loci, we list not only the top colocalized variant but also all variants within the 95% credible causal set.

Motivated by the latest release of coloc, and described here:
https://chr1swallace.github.io/coloc/articles/a03_enumeration.html


Contact
==========

For any queries, please create an issue and we'll respond. Otherwise, please e-mail

Sourya Bhattacharyya: sourya@lji.org

Ferhat Ay: ferhatay@lji.org






