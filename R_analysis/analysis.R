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
library

library(tidyr)

source("./script/custom_function.R") # Import your R script

set.seed(123)

input_dir_path <- "./data/neutral_01"


wo_MAF <- process_vcf(input_dir_path)
MAF <-  maf(wo_MAF)

wo_MAF_freq <- get_freq(wo_MAF)


normally_dist <- norm_phenotype(wo_MAF)
#fst_result <- wc(sampled_df)
#fst_value <- fst_result$FST
#fst_value
#str(dnum_t)

