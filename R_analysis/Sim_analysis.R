# Master project analysis script#

#Package install if needed
if (!require("vcfR")) {
  install.packages("vcfR")
  library(vcfR)
}
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
data_numeric <- as.data.frame(lapply(data[,-1], convert_genotype))
data_numeric

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

#extract 1000 loci
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
vcf_list[[1]]


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



# Normal distribution of effect size based on all the loci
effects <- rnorm(nrow(genotype_numeric[[1]]), mean = 0, sd = 1)


#Calculate phenotype
phenotype_list <- list()
for (i in seq(1:length(genotype_numeric))){
  phenotype_list[[i]] <- rowSums(genotype_numeric[[i]] * effects)
}

### Unsure verify
#Calculate Qst

## Calculate the mean phenotype for each population
mean_phenotypes <- sapply(phenotype_list, mean)

#Calculate the variance within populations
var_within <- mean(sapply(phenotype_list,var))

#Calculate the variance between populations
var_between <- var(mean_phenotypes)

#Calculate Qst

Q_st <- var_between / (var_between + 2*var_within)
Q_st
