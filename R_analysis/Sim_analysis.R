# Master project analysis script#

#Package install if needed
if (!require("vcfR")) {
  install.packages("vcfR")
  library(vcfR)
}
setwd("C:/Users/User/Desktop/Master project/R_analysis")
## reproducibility
set.seed(123)
#### Custom methods
convert_genotype <- function(genotype) {
  alleles <- strsplit(genotype,  "\\|")
  alleles <- lapply(alleles, as.numeric)
  return(unlist(lapply(alleles,sum)))
}


## test function
data <- data.frame(
  locus = c("1_1_1", "1_1_2", "1_1_3"),
  ind1 = c("0|0", "0|1", "1|0"),
  ind2 = c("0|0", "1|1", "0|1"),
  ind3 = c("0|1", "0|0", "1|1"),
  stringsAsFactors = FALSE
)
# Test function
#data_numeric <- as.data.frame(lapply(data[,-1], convert_genotype))
#data_numeric

### MAIN##
# Import files
files <- list.files(path = "./data", pattern = "\\.vcf", full.names = T)

# import and read all the data
data_list <- list()

for (file in files){
  data <- read.vcfR(file)
  data_list[[file]] <- data
}


### Extract data for all the population
## Sample 1000 loci for the first vcf file
vcf_pop1 <- extract.gt(data_list[[1]])

# Ensure the sampled loci are valid for all VCF files
common_loci <- Reduce(intersect, lapply(data_list, function(vcf) rownames(extract.gt(vcf))))

#extract 1000, 100, 10 ,1 
sampled_loci <- sample(common_loci, 1000, replace = FALSE)
# Sample the same loci for each VCF
vcf_list <- list()
for(i in seq(1, length(data_list), by=1)) {
  data <- data_list[[i]]
  vcf_pop <- extract.gt(data)
  
  ## get the same loci for all the population

  sampled_data <- vcf_pop[common_loci, ]
  vcf_list[[i]] <- sampled_data
}
#vcf_list[[1]]


# Row are the loci and column are the individuals
# For the next step my thinking is: alleles at same loci have the same effect size so the cumulative effect should stay the same
# exemple individual 1 has two allele 1 with effect size 1.5. : (1 + 1) * 1.5 = 1 * 1.5 + 1 * 1.5 = 3
#Convert genotype into numeric sum them by individuals
#Bonus convert all the data into dataframe make the calculation way faster
genotype_numeric <- list()
for(i in seq(1, length(vcf_list))){
  data <- as.data.frame(vcf_list[[i]])
  
  genotype_numeric[[i]] <- as.data.frame(lapply(data[, -1], convert_genotype))
  
}

# Sample 100, 10, 1
# To ensure we will have the same row for each population
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
