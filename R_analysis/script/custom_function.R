## Custom function for master project r analysis

if (!require("vcfR")) {
  install.packages("vcfR")
  library(vcfR)
}

library(tidyr)
library(dplyr)
library(ggplot2)


# Process the vcf into df column = loci, row = indiv, with pop attribution
# and genetic info as 0, 1, 2 where 0 = homozygous allele A, 1 = heterozygous, 2 = homozygous allele B
process_vcf <- function(directory) {
  replicate_no <- 1
  processed_data_list <- list()
  
  # List all VCF files in the directory
  vcfs <- list.files(path = directory, full.names = TRUE, pattern = '\\.vcf$')
  population <- rep(1:8, each = 20)  # 8 populations, 20 individuals each
  
  for (vcf in vcfs) {
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
      
      # Convert genotype data to numeric values
      df_numeric <- cbind(
        RowID = df[, 1],  # Keep the first column as is
        as.data.frame(apply(df[, -1], 2, function(col) {
          col <- gsub("0\\|0", "0", col)  # Replace "0|0" with "0"
          col <- gsub("1\\|0|0\\|1", "1", col)  # Replace "1|0" or "0|1" with "1"
          col <- gsub("1\\|1", "2", col)  # Replace "1|1" with "2"
          as.numeric(col)  # Convert to numeric, introduce NA for invalid values
        }))
      )
      
      # Ensure the number of rows matches the original
      if (nrow(df_numeric) != nrow(df)) {
        cat("Mismatch detected!\n")
        cat("Original data frame dimensions:", nrow(df), ncol(df), "\n")
        cat("Numeric data frame dimensions:", nrow(df_numeric), ncol(df_numeric), "\n")
        stop("Mismatch between row names and number of rows in df_numeric!")
      }
      
      
      # Transform df into df usable with hierfstat
      df_numeric <- t(df_numeric) #trasnpose
      df_numeric <- as.data.frame(df_numeric)
      names(df_numeric) <- lapply(df_numeric[1,], as.character) # locus as column
      df_numeric <- df_numeric[-1,] # remove locus row
      df_numeric[sapply(df_numeric, is.character)] <- lapply(df_numeric[sapply(df_numeric, is.character)], as.numeric) #convert data chr to numeric
      df_numeric <- df_numeric %>% #add population column
        mutate(pop = rep(1:8, each = 20))
      df_numeric <- df_numeric[, c(ncol(df_numeric), 1:(ncol(df_numeric) - 1))]
      
      # Store the processed data frame in the list
      
      processed_data_list[[paste0("replicate ",replicate_no)]] <- df_numeric
      replicate_no <- replicate_no + 1
      
    }, error = function(e) {
      message(paste("Error processing", vcf, ":", e$message))
    })
  }
  
  if (length(processed_data_list) == 0) {
    warning("No data was processed successfully.")
    return(NULL)
  }
  
  gc()  # Trigger garbage collection
  return(processed_data_list)
}

### Minor allele sampling
maf <- function(list_df, threshold = 0.05) {
  # Initialize a list to store filtered data frames
  all_filtered_data <- list()
  
  # Loop through each data frame in the list
  for (i in seq_along(list_df)) {
    df <- list_df[[i]]  # Get the current data frame
    
    # Exclude the first column (assumed to be non-locus data, e.g., population)
    loci_data <- df[, -1]
    
    # Calculate MAF for each locus
    meta_maf <- apply(loci_data, 2, function(col) {
      n_individuals <- length(col)
      count_0 <- sum(col == 0, na.rm = TRUE)
      count_1 <- sum(col == 1, na.rm = TRUE)
      count_2 <- sum(col == 2, na.rm = TRUE)
      
      # Calculate the frequency of allele A and B
      freq_A <- (count_0 * 2 + count_1) / (2 * n_individuals)
      freq_B <- 1 - freq_A
      
      return(min(freq_A, freq_B))  # Return the minor allele frequency
    })
    
    # Filter loci based on MAF threshold
    filtered_data_numeric <- loci_data[, meta_maf > threshold, drop = FALSE]
    
    # Combine the filtered loci with the population column
    filtered_df <- cbind(df[, 1, drop = FALSE], filtered_data_numeric)
    
    # Store the filtered data frame
    all_filtered_data[[paste0("replicate ",i)]] <- filtered_df
  }
  gc()
  # Return the list of filtered data frames
  return(all_filtered_data)
}

# normally distributed phenotype

norm_phenotype <- function(list_df, output_dir = "plots") {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  result_list_df <- list()
  allele_frequencies <- list()
  
  for (i in seq_along(list_df)) {
    # Extract loci data (excluding the population column)
    loci_data <- list_df[[i]][, -1]
    
    # Generate random effect sizes for 1, 10, 100, 1000 loci
    effect_size <- c(
      rnorm(1),       # Effect sizes for 1 locus
      rnorm(10),      # Effect sizes for 10 loci
      rnorm(100),     # Effect sizes for 100 loci
      rnorm(1000)     # Effect sizes for 1000 loci
    )
    
    # Phenotype calculations and allele frequency for different sampling sizes
    phenotypes <- list()
    freq_list <- list()
    
    for (n_loci in c(1, 10, 100, 1000)) {
      # Sample loci (use min() to handle cases where there are fewer loci available)
      sampled_loci <- sample(colnames(loci_data), min(n_loci, ncol(loci_data)))
      
      # Subset the data frame to include the sampled loci
      sampled_data <- loci_data[, sampled_loci, drop = FALSE]
      
      # Ensure sampled_data is numeric
      sampled_data <- as.matrix(sampled_data)
      sampled_data <- apply(sampled_data, 2, as.numeric)
      
      # Subset effect sizes to match the number of sampled loci
      sampled_effect_sizes <- effect_size[1:ncol(sampled_data)]
      
      # Calculate the phenotype as a weighted sum of loci values and effect sizes
      phenotype <- as.numeric(sampled_data %*% sampled_effect_sizes)
      phenotypes[[paste0(n_loci, "_loci")]] <- phenotype
      
      # Calculate allele frequencies for the sampled loci
      freqs <- colMeans(sampled_data, na.rm = TRUE) / 2  # Divide by 2 for diploid data
      freq_list[[paste0(n_loci, "_loci")]] <- freqs
      
      # Create a data frame for plotting
      freq_df <- data.frame(
        locus = names(freqs),
        frequency = freqs
      )
      
      # Determine plot width dynamically
      plot_width <- max(10, n_loci / 100)  # Scale width based on the number of loci
      
      # Generate the plot with vertical locus labels
      p <- ggplot(freq_df, aes(x = locus, y = frequency)) +
        geom_point(color = "blue") +
        theme_minimal() +
        labs(
          title = paste("Allele Frequencies of non variant allele for", n_loci, "Loci"),
          x = "Loci",
          y = "Allele Frequency"
        ) +
        theme(
          axis.text.x = element_text(angle = 90)
        )
      
      # Save the plot
      plot_file <- file.path(output_dir, paste0("replicate_",i, " allele_frequencies_", n_loci, "_loci.png"))
      ggsave(plot_file, plot = p, width = 10, height = 6)
    }
    
    # Combine phenotypes into a single data frame with population
    result_df <- data.frame(
      pop = list_df[[i]][, 1],
      do.call(cbind, phenotypes)  # Combine all phenotype columns
    )
    
    # Store the result and allele frequencies
    result_list_df[[i]] <- result_df
    allele_frequencies[[i]] <- freq_list
  }
  
  return(result_list_df)
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
