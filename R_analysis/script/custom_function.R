## Custom function for master project r analysis

if (!require("vcfR")) {
  install.packages("vcfR")
  library(vcfR)
}

## import vcf gt data into df
import_vcf <- function(directory) {
  replicate_no <- 1
  data_list <- list()
  
  for (vcf in directory) {
    print(paste0("Processing: ", vcf))
    if (!file.exists(vcf)) {
      warning(paste("File does not exist:", vcf))
      next
    }
    
    tryCatch({
      # Read VCF file
      data <- vcfR::read.vcfR(vcf)
      print(paste0("File read successfully: ", vcf))
      
      # Extract genotype information (matrix)
      genotype_inf <- vcfR::extract.gt(data)
      print(paste0("Genotype data extracted for: ", vcf))
      
      # Convert matrix to data frame
      df <- as.data.frame(genotype_inf)
      
      # Add row names as a separate column
      df <- data.frame(RowID = rownames(df), df, row.names = NULL)
      
      # Rename columns for clarity
      colnames(df)[-1] <- paste0("i", seq_len(ncol(df) - 1))
      
      # Store the processed data frame in the list
      data_list[[replicate_no]] <- df
      names(data_list)[replicate_no] <- basename(vcf)
      replicate_no <- replicate_no + 1
    }, error = function(e) {
      message(paste("Error processing", vcf, ":", e$message))
    })
  }
  
  if (length(data_list) == 0) {
    warning("No data was processed successfully.")
    return(NULL)
  }
  
  return(data_list)
}


### Convert GT information to numeric (0|0 = 0, 0|1 or 1|0 = 1 and 1|1 = 2)
data_to_numeric <- function(df) {
  # If empty
  if (is.null(df) || nrow(df) == 0) {
    stop("Input data frame is empty or invalid!")
  }
  
  # Replace genotypes with numeric values
  df_numeric <- cbind(
    RowID = df[, 1],  # Keep the first column as is
    as.data.frame(apply(df[, -1], 2, function(col) {  # Apply conversion only to the remaining columns
      col <- gsub("0\\|0", "0", col)  # Replace "0|0" with "0"
      col <- gsub("1\\|0|0\\|1", "1", col)  # Replace "1|0" or "0|1" with "1"
      col <- gsub("1\\|1", "2", col)  # Replace "1|1" with "2"
      as.numeric(col)  # Convert to numeric, introduce NA for invalid values
    }))
  )
  
  # Ensure same number of row as initial
  if (nrow(df_numeric) == nrow(df)) {
    rownames(df_numeric) <- rownames(df)  # Assign row names
  } else {
    # Print diagnostic information for debugging
    cat("Mismatch detected!\n")
    cat("Original data frame dimensions:", nrow(df), ncol(df), "\n")
    cat("Numeric data frame dimensions:", nrow(df_numeric), ncol(df_numeric), "\n")
    stop("Mismatch between row names and number of rows in df_numeric!")
  }
  
  return(df_numeric)
}

### Minor allele sampling
maf <- function(df ,threshold = 0.05) {
  filtered_data_numeric <- list()
  
  # Calculate MAF for each locus
  meta_maf <- apply(df[, -1], 1, function(row) {  # Exclude the first column (e.g., RowID)
    n_individuals <- length(row)
    count_0 <- sum(row == 0, na.rm = TRUE)
    count_1 <- sum(row == 1, na.rm = TRUE)
    count_2 <- sum(row == 2, na.rm = TRUE)
    
    # Calculate the frequency of allele A and B
    freq_A <- (count_0 * 2 + count_1) / (2 * n_individuals)
    freq_B <- 1 - freq_A
    
    return(min(freq_A, freq_B))  # Return the minor allele frequency
  })
  
  # Debug: Print the calculated MAF values
  #print("Calculated MAF values:")
  #print(meta_maf)
  
  # Identify loci to keep based on the threshold
  loci_to_keep <- which(meta_maf > threshold)
  
  # Debug: Print loci to keep
  #print("Loci to keep (indices):")
  #print(loci_to_keep)
  
  # Subset the data frame to keep selected loci
  filtered_df <- df[loci_to_keep, , drop = FALSE]  # Ensure all dimensions remain intact
  
  # Debug: Print the filtered data
  #print("Filtered data:")
  #print(filtered_df)
  
  # Store the filtered data in the list
  filtered_data_numeric <- filtered_df
  
  return(filtered_data_numeric)
}

#split our df into the 8 populations
split_into_pop <- function(df) {
  # Population range
  population_ranges <- list(
    pop1 = 2:21,   # i1 to i20
    pop2 = 22:41,  # i21 to i40
    pop3 = 42:61,  # i41 to i60
    pop4 = 62:81,  # i61 to i80
    pop5 = 82:101, # i81 to i100
    pop6 = 102:121,# i101 to i120
    pop7 = 122:141,# i121 to i140
    pop8 = 142:161 # i141 to i160
  )
  
  # Split genotype info into the 8 populations
  df_list <- lapply(names(population_ranges), function(pop_name) {
    # Extract columns corresponding to the population range
    cols <- population_ranges[[pop_name]]
    pop_df <- df[, c(1, cols), drop = FALSE]  # Include the first column along with the population range
    return(pop_df)
  })
  
  # Assign names to the list for reference purposes
  names(df_list) <- names(population_ranges)
  
  return(df_list)
}


#calculate phenotype
phenotype_calc_norm <- function(df, loci_number) {
  
  # Define which loci will be sampled
  sampled_loci <- sample(nrow(df[[1]]), loci_number, replace = FALSE)
  
  # Generate a normal distribution for effect sizes
  effect_size <- rnorm(loci_number) # Normal distribution effect size
  
  
  # Calculate phenotype for each population
  phenotype <- lapply(df, function(pop_data) { 
    # Ensure pop_data is numeric, excluding the first column
    pop_data_numeric <- as.matrix(pop_data[, -1])  # Exclude the first column (e.g., RowID)
    
    # Subset sampled loci
    sampled_data <- pop_data_numeric[sampled_loci, , drop = FALSE]
    
    # Calculate phenotype (weighted sum)
    colSums(sampled_data * effect_size)
  })
  
  return(phenotype)
}
# uniform distribution
phenotype_calc_uni <- function(df, loci_number) {
  
  # Define which loci will be sampled
  sampled_loci <- sample(nrow(df[[1]]), loci_number, replace = FALSE)
  
  # Generate an uniform distribution for effect sizes
  effect_size <- runif(loci_number, min = -0.5, max = 0.5)
  
  
  # Calculate phenotype for each population
  phenotype <- lapply(df, function(pop_data) { 
    # Ensure pop_data is numeric, excluding the first column
    pop_data_numeric <- as.matrix(pop_data[, -1])  # Exclude the first column (e.g., RowID)
    
    # Subset sampled loci
    sampled_data <- pop_data_numeric[sampled_loci, , drop = FALSE]
    
    # Calculate phenotype (weighted sum)
    colSums(sampled_data * effect_size)
  })
  
  return(phenotype)
}

# L-shaped distribution
phenotype_calc_l <- function(df, loci_number) {
  
  # Define which loci will be sampled
  sampled_loci <- sample(nrow(df[[1]]), loci_number, replace = FALSE)
  
  # Generate an uniform distribution for effect sizes
  effect_size <-rpareto(loci_number, scale = 0.1, shape = 2)
  
  
  # Calculate phenotype for each population
  phenotype <- lapply(df, function(pop_data) { 
    # Ensure pop_data is numeric, excluding the first column
    pop_data_numeric <- as.matrix(pop_data[, -1])  # Exclude the first column (e.g., RowID)
    
    # Subset sampled loci
    sampled_data <- pop_data_numeric[sampled_loci, , drop = FALSE]
    
    # Calculate phenotype (weighted sum)
    colSums(sampled_data * effect_size)
  })
  
  return(phenotype)
}

# Calculate Qst
Qst_cal <- function(pop_phenotype_list) {
  #transform list of pop datafram to matrix 
  phenotype_matrix <- do.call(cbind, pop_phenotype_list)
  
  pop_mean <- colMeans(phenotype_matrix)  # each col is a pop : mean of each population : result a vector
  overall_mean <- mean(pop_mean)     # mean of the resulting vector
  V_B <- mean((pop_mean - overall_mean)^2) #variance between population
  
  V_W <- mean(apply(phenotype_matrix, 2, var)) #compute the variance for each population (apply on column (2) since each column is a population) (diploid) and take the average
  
  # Calculate Qst
  Qst <- V_B / (V_B + 2 * V_W)
  return(Qst)
}

# Flatten the nested structure into a single data frame
flatten_replication_data <- function(replications) {
  do.call(rbind, lapply(seq_along(replications), function(replication_index) {
    replication <- replications[[replication_index]]
    do.call(rbind, lapply(seq_along(replication), function(replicate_index) {
      replicate <- replication[[replicate_index]]
      do.call(rbind, lapply(names(replicate), function(loci) {
        data.frame(
          Replication = paste0("Replication_", replication_index),
          Replicate = paste0("Replicate_", replicate_index),
          Loci = loci,
          QstValue = replicate[[loci]]
        )
      }))
    }))
  }))
}
