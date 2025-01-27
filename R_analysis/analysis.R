## Install required package if not installed
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)

if (!require("kSamples")) install.packages("kSamples")  # Required for ad.test
library(kSamples)

if (!require("VGAM")) install.packages("VGAM")  # L-shaped distribution
library(VGAM)

if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

if (!require("hierfstat")) install.packages("hierfstat")
library(hierfstat)

library(tidyr)

source("./script/custom_function.R") # Import your R script

## Set seed for reproducibility
set.seed(123)
#input dir
input_dir_path <- "./data/neutral_01"

# Get Data from VCF - raw data and with MAF filtering
wo_MAF <- process_vcf(input_dir_path)
MAF <-  maf(wo_MAF)

# Get base allele frequency
wo_MAF_freq <- get_freq(wo_MAF)
maf_freq <- get_freq(maf)

# Normally distributed effect size phenotype calc and allele frequency 

## Without maf
result_norm_wo <- norm_phenotype(wo_MAF)
norm_phenotype_wom <- result_norm_wo$phenotype
norm_allele_fre <- result_norm_wo$allele_fre

## with maf
result_norm_maf <- norm_phenotype(MAF)
norm_phenotype_maf <- result_norm_maf$phenotype
maf_norm_allele_freq <- result_norm_maf$allele_fre

## Calc QST (base equation) normally distributed effect size
qst_norm_wom <- qst_on_list(norm_phenotype_wom)
qst_norm_maf <- qst_on_list(norm_phenotype_maf)


# Uniformly distributed effect size phenotype calc and allele freq
#without maf
result_uni_wo <- uni_phenotype(wo_MAF)
uni_phenotype_wo <- result_uni_wo$phenotype
uni_allele_freq <- result_uni_wo$allele_fre

#with maf
result_uni_maf <- uni_phenotype(MAF)
uni_phenotype_maf <- result_uni_maf$phenotype
uni_allele_freq_maf <- result_uni_maf$allele_fre

# Calc QST
qst_uni_wo <- qst_on_list(uni_phenotype_wo)
qst_uni_maf <- qst_on_list(uni_phenotype_maf)


# L-shaped distributed affect size phenotype calc and allele fre
#without maf
l_result_wo <- l_phenotype(wo_MAF)
l_phenotype_wo <- l_result_wo$phenotype
l_allele_freq_wo <- l_result_wo$allele_fre

# with maf
l_result_maf <- l_phenotype(MAF)
l_phenotype_maf <- l_result_maf$phenotype
l_allele_freq_maf <- l_result_maf$allele_fre

# Calc QST
l_qst_wo <- qst_on_list(l_phenotype_wo)
l_qst_maf <- qst_on_list(l_phenotype_maf)


# Effect size link to allele frequency

#fst_result <- wc(sampled_df)
#fst_value <- fst_result$FST
#fst_value
#str(dnum_t)

