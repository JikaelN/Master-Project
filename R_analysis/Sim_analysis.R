# Master project analysis script#

#Package install if needed
if (!require("vcfR")) {
  install.packages("vcfR")
  library(vcfR)
}

## reproducibility
set.seed(123)

# import file and read data
path_to_file <-  "./data/pop1.vcf"
vcf_data <- read.vcfR(path_to_file)

#Sample loci and get genotype
# 1 locus
sample_one <- vcf_data[sample(nrow(vcf_data), 1),]
genotype_one <- extract.gt(sample_one)
#split genotype data into alleles
alleles <- strsplit(genotype_one, "\\|")
# convert allele into numeric
alleles <- lapply(alleles, as.numeric)
#calculate genotype 
genotype_one_num <- unlist(lapply(alleles,sum))



# 10 loci
sample_ten <- vcf_data[sample(nrow(vcf_data), 10),]
genotype_ten <- extract.gt(sample_ten)

# 100 loci
sample_h <- vcf_data[sample(nrow(vcf_data), 100),]
genotype_h <- extract.gt(sample_h)

# 1000 loci
sample_th <- vcf_data[sample(nrow(vcf_data), 1000),]
genotype_th <- extract.gt(sample_th)


# Normal distribution of effect size
effects <- rnorm(2, mean = 0, sd = 1)

#phenotype 
phenotype <- effects[genotype_one_num +1 ]
