# Master project analysis script#

#Package install if needed
if (!require("vcfR")) {
  install.packages("vcfR")
  library(vcfR)
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}
if (!require("tidyr")) {
  install.packages("tidyr")
  library(tidyr)
}

## Seed for reproducibility / 
set.seed(123) #Reproducibility
setwd("C:/Users/User/Desktop/Master project/R_analysis")

#### Custom methods
#Convert genotype str to numeric
convert_genotype <- function(genotype) {
  alleles <- strsplit(genotype, "\\|")
  alleles <- lapply(alleles, as.numeric)
  return(sum(unlist(alleles)))
}


### MAIN##
# Import files
files <- list.files(path = "./data", pattern = "\\.vcf", full.names = T)

# import and read all the data (convert it to long-table for faster calculation)
data_list <- list()

for (file in files){
  data <- vcfR2tidy(read.vcfR(file))
  data_list[[file]] <- data
}

# Common loci list
# Extract loci from each long-table format data
loci_list <- lapply(data_list, function(data) unique(data$gt["POS"]))
flat_loci_list <- unlist(loci_list, recursive = FALSE)
common_loci <- Reduce(intersect, flat_loci_list)



# Get all the genomic data with common loci for all the 8 population
list_data_filtered <- list()
for (i in seq(1, length(data_list))){
  data <- data_list[[i]]$gt
  d_filtered <- data[data$POS %in% common_loci,]
  list_data_filtered[[i]] <- d_filtered
  
}

## Convert genetic data into numeric
genotype_numeric <- lapply(data_list, function(data) {
  data$gt$gt_numeric <- sapply(data$gt$gt_GT, convert_genotype)
  return(data)
})


#head(genotype_numeric[[1]]$gt)



### TO CHANGE from
#extract 1000, 100, 10 ,1 

# Sample the same loci for each VCF


# Row are the loci and column are the individuals
# For the next step my thinking is: alleles at same loci have the same effect size so the cumulative effect should stay the same
# exemple individual 1 has two allele 1 with effect size 1.5. : (1 + 1) * 1.5 = 1 * 1.5 + 1 * 1.5 = 3
#Convert genotype into numeric sum them by individuals
#Bonus convert all the data into dataframe make the calculation way faster


# Sample 1000, 100, 10, 1
# To ensure we will have the same row for each population
sampled_indices_1000 <- sample(seq_len(nrow(genotype_numeric[[1]])), 1000, replace = FALSE)
sampled_indices_100 <- sample(seq_len(nrow(genotype_numeric[[1]])), 100, replace = FALSE)
sampled_indices_10 <- sample(seq_len(nrow(genotype_numeric[[1]])), 10, replace = FALSE)
sampled_indices_1 <- sample(seq_len(nrow(genotype_numeric[[1]])), 1, replace = FALSE)

# Sample for each population

sampled_phenotype_100 <- list()
sampled_phenotype_10 <- list()
sampled_phenotype_1 <- list()

for (i in seq(1, length(genotype_numeric))){
  sampled_phenotype_100[[i]] <- genotype_numeric[[i]][sampled_indices_100, ]
  sampled_phenotype_10[[i]] <- genotype_numeric[[i]][sampled_indices_10, ]
  sampled_phenotype_1[[i]] <- genotype_numeric[[i]][sampled_indices_1, ]
}

## reproducibility
set.seed(123)

# Normal distribution of effect size based on all the loci
effects <- rnorm(nrow(genotype_numeric[[1]]), mean = 0, sd = 1)
effects_100 <- rnorm(nrow(sampled_phenotype_100[[1]]), mean = 0, sd = 1)
effects_10 <- rnorm(nrow(sampled_phenotype_10[[1]]), mean = 0, sd = 1)
effects_1 <- rnorm(nrow(sampled_phenotype_10[[1]]), mean = 0, sd = 1)



#Calculate phenotype
phenotype_list <- list()
pheno_100 <- list()
pheno_10 <- list()
pheno_1 <- list()
for (i in seq(1:length(genotype_numeric))){
  phenotype_list[[i]] <- colSums(genotype_numeric[[i]] * effects)
  pheno_100[[i]] <- colSums(sampled_phenotype_100[[i]] * effects_100)
  pheno_10[[i]] <- colSums(sampled_phenotype_10[[i]] * effects_10)
  pheno_1[[i]] <- colSums(sampled_phenotype_1[[i]] * effects_1)
  
}


### Unsure verify
#Calculate Qst

## Calculate the mean phenotype for each population
mean_phenotypes <- sapply(phenotype_list, mean) # 1000
mean_ph_100 <- sapply(pheno_100, mean)
mean_ph_10 <- sapply(pheno_10, mean)
mean_ph_1 <- sapply(pheno_1, mean)

#Calculate the variance within populations
var_within <- mean(sapply(phenotype_list,var))
var_wt_100 <- mean(sapply(pheno_100, var))
var_wt_10 <- mean(sapply(pheno_10, var))
var_wt_1 <- mean(sapply(pheno_1, var))

#Calculate the variance between populations
var_between <- var(mean_phenotypes)
var_bt_100 <- var(mean_ph_100)
var_bt_10 <- var(mean_ph_10)
var_bt_1 <- var(mean_ph_1)


#Calculate Qst

Q_st <- var_between / (var_between + 2*var_within)
Q_st_100 <- var_bt_100 / (var_bt_100 + 2*var_wt_100)
Q_st_10 <- var_bt_10 / (var_bt_10 + 2*var_wt_10)
Q_st_1 <- var_bt_1 / (var_bt_1 + 2*var_wt_1)


qqplot(Q_st, Q_st_100)
## Now i need to replicate a lot this experiment to be able to have QQplot with more that 1 point :D
