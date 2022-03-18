#!/usr/bin/env Rscript

#========================
# this script is for colocalization between eQTL and GWAS summary statistics
#========================
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(coloc))
suppressMessages(library(snpStats))
suppressMessages(library(dplyr))
# library(ggplot2)
# library(ggplotify)
suppressMessages(library(yaml))

options(scipen = 10)
options(datatable.fread.datatable=FALSE) 

##======== pp4 threshold for colocalization
##======== not used here - but still worth to note it
THR_POST_PROB <- 0.75

##===== within chr6, excluded MHC region from consideration
##===== as suggested in https://www.nature.com/articles/nature22403
MHC_chr6_LB <- 28477897
MHC_chr6_UB <- 33448354

##===========
## 500 Kb on both side of a significant GWAS SNP
GWAS_WINDOW <- 500000

## GWAS SE value (small floating point value) for zero entries
DUMMY_GWAS_SE_VAL <- 0.00001

# #=======================
# # this function estimates standard error of SNP using beta, MAF values
# ## parameters:
# ## b: beta
# ## p: MAF (or AF)
# ## n: sample size
# # references: 1) https://www.biostars.org/p/276869/ (last paragraph), 
# # 2) https://www.nature.com/articles/ng.3538, 3) https://www.biostars.org/p/319584/ (last paragraph)
# # here b = z / sqrt(2p(1-p)(n+z^2)); where b = beta (regression coefficient), z = z score, p = MAF n = sample size
# # so, z can be estimated as  2p(1-p)(b^2)n / sqrt(1 - 2p(1-p)b^2)
# # finally standard error is estimated as beta / z
# #=======================
# # estimate_SE <- function(inpdf, beta_column, MAF_column, size_column) { 	
# estimate_SE <- function(b, p, n) { 	

#  	# b <- inpdf[, beta_column]
#  	# p <- inpdf[, MAF_column]
#  	# n <- inpdf[, size_column]

#  	numerator <- (2 * p * (1-p) * (b ^ 2) * n)
#  	squared_denominator <- (1 - 2 * p * (1-p) * (b ^ 2))

#  	# Estimated_Z <- rep(0, nrow(inpdf)) 	
#  	Estimated_Z <- rep(0, length(numerator))
#  	zero_denom_idx <- which(squared_denominator == 0)
#  	# nonzero_denom_idx <- setdiff(seq(1,nrow(inpdf)), zero_denom_idx)
#  	nonzero_denom_idx <- setdiff(seq(1,length(numerator)), zero_denom_idx)
#  	Estimated_Z[nonzero_denom_idx] <- numerator[nonzero_denom_idx] / sqrt(abs(squared_denominator[nonzero_denom_idx]))
#  	Estimated_Z[zero_denom_idx] <- numerator[zero_denom_idx] 	# avoid zero division 

#  	# now estimate the standard error (se) which is basically beta / z
#  	zero_Z_idx <- which(Estimated_Z == 0)
#  	nonzero_Z_idx <- setdiff(seq(1,length(Estimated_Z)), zero_Z_idx)
#  	Estimated_SE <- rep(0, length(Estimated_Z))
#  	Estimated_SE[nonzero_Z_idx] <- b[nonzero_Z_idx] / Estimated_Z[nonzero_Z_idx]
#  	Estimated_SE[zero_Z_idx] <- b[zero_Z_idx]

#  	# standard error
#  	# outdf <- data.frame(SE=Estimated_SE)
#  	# testing with variance 
#  	# check https://sciencing.com/calculate-variance-standard-error-6372721.html 	
#  	outdf <- data.frame(SE=(Estimated_SE * Estimated_SE * n))

#  	return(outdf)
# }

# *************************************************************
# ******** input parameters ***********
# *************************************************************
args <- commandArgs(TRUE)

configfile <- args[1]

## read the configuration file
config <- yaml.load_file(configfile)

## create the output directory
system(paste("mkdir -p", config$OutDir))

##==== a few global parameters
p12_val <- 0.00001 	# 1e-5 
bool_Coloc_Summary_DF <- FALSE

## output summary file to contain the colocalized gene - SNP pairs
ColocSNPInfoFile <- paste0(config$OutDir, '/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed')

## output summary file to contain the colocalized gene - SNP pairs - 95% credible sets
ColocSNPInfoFile_95pct_credible <- paste0(config$OutDir, '/FINAL_Summary_Coloc_Gene_SNP_Pairs_95pct_credible_set.bed')

## output log file
textfile <- paste0(config$OutDir, '/Out_Summary.log')
sink(textfile)
cat(sprintf("\n\n Parameters for colocalization (GWAS) analysis **** \n input GWAS summary statistics file : %s \n reference EQTL file : %s \n BaseOutDir : %s ", config$GWASFile, config$eQTLFile, config$OutDir))

##=============
##==== read input GWAS file
##=============
GWASData <- data.table::fread(config$GWASFile, header=T)

collist_gwas <- c(config$GWASchrCol, config$GWASposCol)
colnameslist_gwas <- c('chr', 'pos')
if (length(config$GWASbetaCol) == 0) {
	GWASData$beta <- log(GWASData[, config$GWASORCol])
	collist_gwas <- c(collist_gwas, ncol(GWASData))
} else {
	collist_gwas <- c(collist_gwas, config$GWASbetaCol)
}
colnameslist_gwas <- c(colnameslist_gwas, 'beta_GWAS')
if (length(config$GWASSECol) > 0) {
	collist_gwas <- c(collist_gwas, config$GWASSECol)
	colnameslist_gwas <- c(colnameslist_gwas, 'SE_gwas')
}
collist_gwas <- c(collist_gwas, config$GWASpValCol)
colnameslist_gwas <- c(colnameslist_gwas, 'pval_GWAS')
if (length(config$GWASAFCol) > 0) {
	collist_gwas <- c(collist_gwas, config$GWASAFCol)
	colnameslist_gwas <- c(colnameslist_gwas, 'AF_GWAS')
}
if (length(config$GWASSampleSizeCol) > 0) {
	collist_gwas <- c(collist_gwas, config$GWASSampleSizeCol)
	colnameslist_gwas <- c(colnameslist_gwas, 'N_GWAS')
}
GWASData <- GWASData[, c(collist_gwas)]
colnames(GWASData) <- colnameslist_gwas
cat(sprintf("\n *** Number of reference GWAS SNPs: %s ", nrow(GWASData)))

## if the standard error in GWAS data is not provided, we can compute them
## using p-value (2-tailed distribution) and beta
## check: https://www.biostars.org/p/431875/
if (length(config$GWASSECol) == 0) {
	GWASData$SE_gwas <- abs(GWASData$beta_GWAS / qnorm(GWASData$pval_GWAS / 2))
}

## check the "SE_gwas" column - there should not be any entry with 0 standard error
## use a small floating point value (say 0.00001)
idx <- which(GWASData$SE_gwas == 0)
if (length(idx) > 0) {
	GWASData[idx, c("SE_gwas")] <- DUMMY_GWAS_SE_VAL
}

##======= if the GWAS data has chromosomes with numbers
##======= then append the "chr"
if (grepl("chr", GWASData[1,1]) == FALSE) {
	GWASData[,1] <- paste0("chr", GWASData[,1])
}

##====== keep only significant GWAS entries
## currently commented - use the full GWAS summary statistics
if (0) {
	GWASData <- GWASData[which(GWASData$pval_GWAS < 5e-8), ]
	cat(sprintf("\n *** Number of reference GWAS SNPs - p-value < 5e-8 : %s ", nrow(GWASData)))
}

##=============
##==== read input eQTL file
##=============
eQTLData <- data.table::fread(config$eQTLFile, header=T)
collist_eqtl <- c(config$eQTLchrCol, config$eQTLposCol, config$eQTLgeneIDCol)
colnameslist_eqtl <- c('chr', 'pos', 'eGeneID')
if (length(config$eQTLgeneNameCol) > 0) {
	collist_eqtl <- c(collist_eqtl, config$eQTLgeneNameCol)
	colnameslist_eqtl <- c(colnameslist_eqtl, 'eGeneName')
}
collist_eqtl <- c(collist_eqtl, config$eQTLrsIDCol)
colnameslist_eqtl <- c(colnameslist_eqtl, 'rsID')
if (length(config$eQTLTSSDistCol) > 0) {
	collist_eqtl <- c(collist_eqtl, config$eQTLTSSDistCol)
	colnameslist_eqtl <- c(colnameslist_eqtl, 'TSSDist')
}
collist_eqtl <- c(collist_eqtl, config$eQTLbetaCol)
colnameslist_eqtl <- c(colnameslist_eqtl, 'beta_eQTL')
if (length(config$eQTLSECol) > 0) {
	collist_eqtl <- c(collist_eqtl, config$eQTLSECol)
	colnameslist_eqtl <- c(colnameslist_eqtl, 'SE_eQTL')
}
if (length(config$eQTLAFCol) > 0) {
	collist_eqtl <- c(collist_eqtl, config$eQTLAFCol)
	colnameslist_eqtl <- c(colnameslist_eqtl, 'AF_eQTL')
}
collist_eqtl <- c(collist_eqtl, config$eQTLpValCol)
colnameslist_eqtl <- c(colnameslist_eqtl, 'pval_eQTL')
if (length(config$eQTLFDRCol) > 0) {
	collist_eqtl <- c(collist_eqtl, config$eQTLFDRCol)
	colnameslist_eqtl <- c(colnameslist_eqtl, 'FDR_eQTL')
}

eQTLData <- eQTLData[, c(collist_eqtl)]
colnames(eQTLData) <- colnameslist_eqtl
cat(sprintf("\n *** Number of reference eQTLs: %s ", nrow(eQTLData)))

## if the standard error in eQTL data is not provided, we can compute them
## using p-value (2-tailed distribution) and beta
## check: https://www.biostars.org/p/431875/
if (length(config$eQTLSECol) == 0) {
	eQTLData$SE_eQTL <- abs(eQTLData$beta_eQTL / qnorm(eQTLData$pval_eQTL / 2))
}

## check for any zero entries in "SE_eQTL" column
## and replace them with dummy floating point values
idx <- which(eQTLData$SE_eQTL == 0)
if (length(idx) > 0) {
	eQTLData[idx, c("SE_eQTL")] <- DUMMY_GWAS_SE_VAL
}

##======= if the eQTL data has chromosomes with numbers
##======= then append the "chr"
if (grepl("chr", eQTLData[1,1]) == FALSE) {
	eQTLData[,1] <- paste0("chr", eQTLData[,1])
}

##======================
##=== divide the genome into a fixed set of loci
##=== lead variant: 500 Kb region on both side (1 Mb locus)
##======================

GWASChrList <- as.vector(sort(unique(GWASData$chr)))
for (chridx in 1:length(GWASChrList)) {
	currchr <- GWASChrList[chridx]	
	
	## GWAS data for current chromosome
	GWASData_currchr <- GWASData[which(GWASData$chr == currchr), ]
	cat(sprintf("\n\n\n *********** ==>> processing chromosome : %s \n number of GWAS SNP entries : %s ********* ", currchr, nrow(GWASData_currchr)))
	
	## eQTL data for current chromosome
	eQTLData_currchr <- eQTLData[which(eQTLData$chr == currchr), ]
	cat(sprintf("\n ******** Number of reference eQTLs : %s ", nrow(eQTLData_currchr)))

	## append the AF in the eQTL data if both GWAS and eQTL data do not have allele frequency values
	bool_use_snpinfo_data <- FALSE
	if ((length(config$GWASAFCol) == 0) & (length(config$eQTLAFCol) == 0)) {
		## custom SNP information file for the current chromosome
		SNPInfoFile <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2018_HiChIP_FiveImmuneCell_Vivek/Data/Genotype_Ref_May_2020/SNPInfo/snpinfo_', currchr, '.txt')
		if (file.exists(SNPInfoFile)) {
			cat(sprintf("\n ** Both eQTL and GWAS input files do not have allele frequency information \n *** however, we are using a custom SNP information file to incorporate the allele frequency information : %s ", SNPInfoFile))
			SNPInfoData <- data.table::fread(SNPInfoFile, header=T, sep=" ")
			## colnames(SNPInfoData) <- c('chr', 'pos', 'variant_id', 'rs_id', 'ref', 'alt', 'AC', 'AF', 'AN')
			SNPInfoData <- SNPInfoData[, c(1,2,(ncol(SNPInfoData)-1))]	#, ncol(SNPInfoData))]
			colnames(SNPInfoData) <- c('chr', 'pos', 'AF_eQTL')	#, 'AN_eQTL')
			## merge this SNP information with eQTLData_currchr
			eQTLData_currchr <- dplyr::inner_join(eQTLData_currchr, SNPInfoData)
			## set the boolean variable true
			bool_use_snpinfo_data <- TRUE
			# ## estimate SE using the beta, AF and AN values
			# est_SE_DF <- estimate_SE(eQTLData_currchr$beta_eQTL, eQTLData_currchr$AF_eQTL, eQTLData_currchr$AN_eQTL)
			# eQTLData_currchr$est_SE_eQTL <- est_SE_DF$SE
		}
	}

	##======= define the set of GWAS loci for this chromosome
	## first, sort the GWAS entries by p-values	
	gwasdata <- GWASData_currchr[order(GWASData_currchr$pval_GWAS), ]
	
	## after the while loop, "GWAS_Loci_DF" contains the set of GWAS loci for this chromosome
	bool_GWAS_Loci <- FALSE
	while(1) {
		## the top most GWAS entry should be significant
		## 5th field stores the p-value
		if (gwasdata[1, c("pval_GWAS")] > 5e-8) {
			break
		}
		## the top most entry defines the current GWAS loci
		## significant GWAS SNP and 500 Kb in both side
		cat(sprintf("\n using GWAS SNP - chr : %s pos : %s ", currchr, gwasdata[1, 2]))		
		startpos <- max(0, (gwasdata[1, c("pos")] - GWAS_WINDOW))
		endpos <- (gwasdata[1, c("pos")] + GWAS_WINDOW)
		currLociDF <- data.frame(chr=currchr, start=startpos, end=endpos)
		if (bool_GWAS_Loci == FALSE) {
			GWAS_Loci_DF <- currLociDF
			bool_GWAS_Loci <- TRUE
		} else {
			GWAS_Loci_DF <- rbind.data.frame(GWAS_Loci_DF, currLociDF)
		}
		## discard any GWAS entry belonging to this interval
		idx <- which((gwasdata[, c("pos")] >= startpos) & (gwasdata[, c("pos")] <= endpos))
		gwasdata <- gwasdata[-idx, ]
		if (nrow(gwasdata) == 0) {
			break
		}
	}
	cat(sprintf("\n ===>> number of rows in GWAS_Loci_DF (GWAS loci for the current chromosome) : %s ", nrow(GWAS_Loci_DF)))	

	##========== process individual GWAS loci and corresponding eQTL files for colocalization
	for (lociidx in 1:nrow(GWAS_Loci_DF)) {
		startpos <- GWAS_Loci_DF[lociidx, 2]
		endpos <- GWAS_Loci_DF[lociidx, 3]
		cat(sprintf("\n *** processing GWAS loci ---  chr : %s start : %s end : %s ", currchr, startpos, endpos))
		## gwas data for the current loci
		currloci_GWASdata <- GWASData_currchr[which((GWASData_currchr[, c("pos")] >= startpos) & (GWASData_currchr[, c("pos")] <= endpos)), ]
		## eQTLs for the current loci
		currloci_eqtldata <- eQTLData_currchr[which((eQTLData_currchr[, c("pos")] >= startpos) & (eQTLData_currchr[, c("pos")] <= endpos)), ]
		cat(sprintf("\n *** Number of GWAS SNPs for this loci : %s number of eQTLs for this loci : %s ", nrow(currloci_GWASdata), nrow(currloci_eqtldata)))

		##===== common set of SNPs between GWAS and eQTL
		merge_SNP_GWAS_Data <- dplyr::inner_join(currloci_eqtldata, currloci_GWASdata)
		if (nrow(merge_SNP_GWAS_Data) == 0) {
			cat(sprintf("\n !!!! merged eQTL and GWAS data has 0 entries -- proceed to the next GWAS loci !!! "))			
			next
		}

		##===== discard entries within MHC region
		mhc_idx <- which((merge_SNP_GWAS_Data$chr == "chr6") & (merge_SNP_GWAS_Data$pos >= MHC_chr6_LB) & (merge_SNP_GWAS_Data$pos <= MHC_chr6_UB))
		if (length(mhc_idx) > 0) {
			merge_SNP_GWAS_Data <- merge_SNP_GWAS_Data[-mhc_idx, ]
		}
		if (nrow(merge_SNP_GWAS_Data) == 0) {
			cat(sprintf("\n !!!! After discarding the MHC region SNPs - merged eQTL and GWAS data has 0 entries -- proceed to the next GWAS loci !!! "))
			next
		}

		##===== check for non-NA fields
		NA_idx <- which((is.na(merge_SNP_GWAS_Data$pval_GWAS)) | (is.na(merge_SNP_GWAS_Data$pval_eQTL)) | (is.na(merge_SNP_GWAS_Data$beta_GWAS)) | (is.na(merge_SNP_GWAS_Data$SE_gwas)) | (is.na(merge_SNP_GWAS_Data$beta_eQTL)) | (is.na(merge_SNP_GWAS_Data$SE_eQTL)))
		if (length(NA_idx) > 0) {
			merge_SNP_GWAS_Data <- merge_SNP_GWAS_Data[-c(NA_idx), ]
		}
		if (nrow(merge_SNP_GWAS_Data) == 0) {
			cat(sprintf("\n !!!! After discarding the NA entries - merged eQTL and GWAS data has 0 entries -- proceed to the next GWAS loci !!! "))
			next
		}

		# ## check the minimum GWAS p-value of SNPs 
		# ## at least one SNP should be significant in GWAS
		# min_p_val_gwas <- min(merge_SNP_GWAS_Data$pval_GWAS)
		# if (min_p_val_gwas > 5e-8) {
		# 	cat(sprintf("\n !!!! The minimum GWAS p-value for the SNPs in the merge_SNP_GWAS_Data : %s - not GWAS significant - discard this loci ", min_p_val_gwas))
		# 	next 
		# }

		cat(sprintf("\n *** Number of entries in merge_SNP_GWAS_Data for this loci : %s ", nrow(merge_SNP_GWAS_Data)))		

		##=== dump the input SNPs for colocalization analysis in the specified output directory 
		##=== for processing the current gene
		CurrLoci_OutDir <- paste0(config$OutDir, '/GWAS_Loci_', currchr, '_', startpos, '_', endpos, "_GWASSNPPos_", (endpos - GWAS_WINDOW))
		system(paste("mkdir -p", CurrLoci_OutDir))

		## write the input GWAS SNP (which defines the current loci)
		write.table(dplyr::inner_join(data.frame(chr=currchr, pos=(endpos - GWAS_WINDOW)), GWASData_currchr), paste0(CurrLoci_OutDir, '/input_topmost_GWAS_SNP_defining_current_loci.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

		## write the merged GWAS and eQTL data
		write.table(merge_SNP_GWAS_Data, paste0(CurrLoci_OutDir, '/merged_SNPs_GWAS_input_colocalization.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

		##=================
		## main colocalization routine
		##=================	
		##==== quant for dataset1
		##====== N: number of SNPs for this gene
		cat(sprintf("\n ==>> Creating coloc dataset for GWAS === "))
		dataset1=list(pvalues=merge_SNP_GWAS_Data$pval_GWAS, type="quant", N=nrow(currloci_GWASdata), beta=merge_SNP_GWAS_Data$beta_GWAS, varbeta=(merge_SNP_GWAS_Data$SE_gwas * merge_SNP_GWAS_Data$SE_gwas))
		## check the p-values which will be derived by coloc 
		## using the beta and varbeta values
		pvalues_est_coloc_gwas <- pnorm( -abs( dataset1$beta/sqrt(dataset1$varbeta) ) ) * 2
		cat(sprintf("\n Minimum p-value for this coloc dataset GWAS : %s ", min(pvalues_est_coloc_gwas)))

		## quant for dataset2
		## N: number of SNPs for this chromosome (used for colocalization)
		cat(sprintf("\n ==>> Creating coloc dataset for eQTL === "))
		dataset2 <- list(pvalues=merge_SNP_GWAS_Data$pval_eQTL, type="quant", N=nrow(currloci_eqtldata), beta=merge_SNP_GWAS_Data$beta_eQTL, varbeta=(merge_SNP_GWAS_Data$SE_eQTL * merge_SNP_GWAS_Data$SE_eQTL))	#, sdY=stdev_gene_expr)		
		## check the p-values which will be derived by coloc 
		## using the beta and varbeta values
		pvalues_est_coloc_eQTL <- pnorm( -abs( dataset2$beta/sqrt(dataset2$varbeta) ) ) * 2
		cat(sprintf("\n Minimum p-value for this coloc dataset eQTL : %s ", min(pvalues_est_coloc_eQTL)))

		# ## debug 
		# ## check the alternative dataset
		# dataset2_alt <- list(pvalues=merge_SNP_GWAS_Data$pval_eQTL, type="quant", N=nrow(currloci_eqtldata), beta=merge_SNP_GWAS_Data$beta_eQTL, varbeta=(merge_SNP_GWAS_Data$est_SE_eQTL * merge_SNP_GWAS_Data$est_SE_eQTL))
		# pvalues_est_coloc_eQTL_alt <- pnorm( -abs( dataset2_alt$beta/sqrt(dataset2_alt$varbeta) ) ) * 2
		# cat(sprintf("\n Minimum p-value for this coloc dataset eQTL (alternative SE using AF and AN) : %s ", min(pvalues_est_coloc_eQTL_alt)))

		## colocalization function according to prior probability settings
		## also check if the allele frequency is provided in the GWAS data
		if (length(config$GWASAFCol) > 0) {			
			cat(sprintf("\n ===>> coloc.abf using AF from GWAS data"))
			result <- coloc.abf(dataset1=dataset1, dataset2=dataset2, MAF=merge_SNP_GWAS_Data$AF_GWAS, p12=p12_val)
		} else if ((length(config$eQTLAFCol) > 0) | (bool_use_snpinfo_data == TRUE)) {
			cat(sprintf("\n ===>> coloc.abf using AF from eQTL data"))
			result <- coloc.abf(dataset1=dataset1, dataset2=dataset2, MAF=merge_SNP_GWAS_Data$AF_eQTL, p12=p12_val)
		} else {
			cat(sprintf("\n ===>> coloc.abf without using AF"))
			result <- coloc.abf(dataset1=dataset1, dataset2=dataset2, p12=p12_val)
		}

		##======== colocalization results - plot sensitivity of the posterior probability computation
		if (0) {
			plotfile <- paste0(CurrLoci_OutDir, '/coloc_abf_res_sensitivity_plot.pdf')
			pdf(plotfile, width=8, height=6)
			coloc::sensitivity(result, rule="H4 > 0.5")
			dev.off()
		}

		##======== colocalization results - two data frames
		##======== let N = nrow(merge_eQTL_GWAS_Data)
		##======== result$summary: summary statistic of posterior probability (considering all SNPs)
		##======== result$results: detailed statistic of posterior probability per SNP
		##======== Column 1: SNP.NUM where NUM = index number of SNP - varies from 1 to N
		coloc_post_prob_summary_file <- paste0(CurrLoci_OutDir, '/coloc_abf_summary_DF1.bed')
		coloc_post_prob_summary_file_2 <- paste0(CurrLoci_OutDir, '/coloc_abf_results_DF2.bed')
		coloc_summary_file <- paste0(CurrLoci_OutDir, '/coloc_abf_results_detailed.bed')

		write.table(result$summary, coloc_post_prob_summary_file, row.names=T, col.names=T, sep="\t", quote=F, append=F)
		write.table(result$results, coloc_post_prob_summary_file_2, row.names=F, col.names=T, sep="\t", quote=F, append=F)

		##======== modify the results (DF2) file to add SNP IDs in addition to the SNP.NUM format
		resDF <- result$results 	
		CN <- colnames(resDF)
		temp_SNP_vec <- as.vector(paste0("SNP.", seq(1, nrow(merge_SNP_GWAS_Data))))
		m <- match(resDF[,1], temp_SNP_vec)
		idx_r <- which(!is.na(m))
		idx_t <- m[!is.na(m)]		
		## when data frame merging was done by chr and pos fields
		resDF <- cbind.data.frame(resDF[idx_r, 1], merge_SNP_GWAS_Data[idx_t, c("rsID", "pos")], resDF[idx_r, 2:ncol(resDF)])
		colnames(resDF) <- c(CN[1], 'rsID', 'pos', CN[2:length(CN)])
		write.table(resDF, coloc_summary_file, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		cat(sprintf("\n creating detailed results -- entries in resDF : %s matching SNP entries : %s ", nrow(resDF), length(idx_r)))		

		##======== check if there is a colocalization, w.r.t the PP4 threshold
		##======== in such a case, document the summary statistics

		##===== posterior probability (of shared variants)
		pp_H0 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==3) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H1 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==4) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H2 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==5) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H3 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==6) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H4 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==7) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		cat(sprintf("\n pp_H0 : %s pp_H1 : %s pp_H2 : %s pp_H3 : %s pp_H4 : %s ", pp_H0, pp_H1, pp_H2, pp_H3, pp_H4))

		##====== colocalization - pp.H4 > THR_POST_PROB - potential causal variant
		## now list the topmost causal variant
		## and 95% credible set causal variants
		## in two separate files
		## follow the instruction: https://chr1swallace.github.io/coloc/articles/a03_enumeration.html
		if ((!is.na(pp_H4)) & (pp_H4 > THR_POST_PROB)) {			
			cat(sprintf("\n\n *** COLOCALIZATION FOUND -------- ( pp.H4 > %s ) ", THR_POST_PROB))

			## for colocalization, check the SNPs and their SNP.PP.H4 field
			## subset for those SNPs which have SNP.PP.H4 > 0.01
			resDF_sub <- subset(resDF, SNP.PP.H4 > 0.01)
			if (nrow(resDF_sub) == 0) {
				next
			}

			## now sort the SNPs based on their SNP.PP.H4
			## such that the most significant SNP is at the first entry
			o <- order(resDF_sub$SNP.PP.H4, decreasing=TRUE)
			## resDF_sub[o, ] prints the SNPs in their sorting order

			## extract 95% credible set information
			cs <- cumsum(resDF_sub$SNP.PP.H4[o])
			w <- which(cs > 0.95)[1]
			if (is.na(w)) {
				w <- nrow(resDF_sub)
			}
			## resDF_sub[o, ][1:w, ] prints the 95% credible SNPs 

			## insert a field "snp" of the format "SNP.NUM" in the merge_SNP_GWAS_Data
			merge_SNP_GWAS_Data$snp <- paste0("SNP.", seq(1, nrow(merge_SNP_GWAS_Data)))

			## get the topmost colocalized SNP in a data frame
			## also append the colocalization summary (H0 to H4) outputs
			topColocDF <- cbind.data.frame(resDF_sub[o, ][1, ], data.frame(coloc_H0=pp_H0, coloc_H1=pp_H1, coloc_H2=pp_H2, coloc_H3=pp_H3, coloc_H4=pp_H4))

			## get the 95% credible set colocalized SNPs in a data frame
			credibleColocDF <- cbind.data.frame(resDF_sub[o, ][1:w, ], data.frame(coloc_H0=rep(pp_H0, w), coloc_H1=rep(pp_H1, w), coloc_H2=rep(pp_H2, w), coloc_H3=rep(pp_H3, w), coloc_H4=rep(pp_H4, w)))

			## merge with the input SNP and gene information
			## from merge_SNP_GWAS_Data data
			topColocDF <- dplyr::inner_join(topColocDF, merge_SNP_GWAS_Data)
			credibleColocDF <- dplyr::inner_join(credibleColocDF, merge_SNP_GWAS_Data)

			## create a summary data frame for all such colocalization outputs 
			if (bool_Coloc_Summary_DF == FALSE) {						
				Final_topColocDF <- cbind.data.frame(data.frame(GWASLoci=rep(paste0(currchr, ":", startpos, "-", endpos), nrow(topColocDF))), topColocDF)
				Final_credibleColocDF <- cbind.data.frame(data.frame(GWASLoci=rep(paste0(currchr, ":", startpos, "-", endpos), nrow(credibleColocDF))), credibleColocDF)			
				bool_Coloc_Summary_DF <- TRUE
			} else {
				Final_topColocDF <- rbind.data.frame(Final_topColocDF, cbind.data.frame(data.frame(GWASLoci=rep(paste0(currchr, ":", startpos, "-", endpos), nrow(topColocDF))), topColocDF))
				Final_credibleColocDF <- rbind.data.frame(Final_credibleColocDF, cbind.data.frame(data.frame(GWASLoci=rep(paste0(currchr, ":", startpos, "-", endpos), nrow(credibleColocDF))), credibleColocDF))				
			}			

			## dump the topColocDF
			topColocDF_OutFile <- paste0(CurrLoci_OutDir, '/coloc_SNP_Set.txt')
			write.table(topColocDF, topColocDF_OutFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

			## dump the credibleColocDF
			credibleColocDF_OutFile <- paste0(CurrLoci_OutDir, '/coloc_95pct_credible_SNP_Set.txt')
			write.table(credibleColocDF, credibleColocDF_OutFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)		
		}

		##===== delete temporary objects
		if (exists("currloci_GWASdata")) {
			rm(currloci_GWASdata)
		}
		if (exists("currloci_eqtldata")) {
			rm(currloci_eqtldata)
		}
		if (exists("merge_SNP_GWAS_Data")) {
			rm(merge_SNP_GWAS_Data)
		}
		if (exists("dataset1")) {
			rm(dataset1)
		}
		if (exists("dataset2")) {
			rm(dataset2)
		}
		if (exists("result")) {
			rm(result)
		}
		if (exists("resDF")) {
			rm(resDF)
		}
		if (exists("topColocDF")) {
			rm(topColocDF)
		}
		if (exists("credibleColocDF")) {
			rm(credibleColocDF)
		}

	}	# end GWAS loci processing loop

	##==== delete temporary objects
	if (exists("GWASData_currchr")) {
		rm(GWASData_currchr)
	}
	if (exists("eQTLData_currchr")) {
		rm(eQTLData_currchr)
	}
	if (exists("SNPInfoData")) {
		rm(SNPInfoData)
	}

}	# end GWAS chromosome loop

if (exists("GWASData")) {
	rm(GWASData)
}
if (exists("eQTLData")) {
	rm(eQTLData)
}

##==================
## print the summary statistics
##==================

# dump the final summary file
if (exists("Final_topColocDF")) {
	if (nrow(Final_topColocDF) > 0) {
		write.table(Final_topColocDF, ColocSNPInfoFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}
}

if (exists("Final_credibleColocDF")) {
	if (nrow(Final_credibleColocDF) > 0) {
		write.table(Final_credibleColocDF, ColocSNPInfoFile_95pct_credible, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}
}

