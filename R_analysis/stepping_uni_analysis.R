# Required Libraries
library("jsonlite")
library("dplyr")
library("ggplot2")
library("tidyr")

set.seed(123)

# Directory Paths
directory_path <- "./data/stepping_sampling/1_sampling"
maf_output_dir <- file.path(directory_path, "uni_maf_qq_plots")
wo_output_dir <- file.path(directory_path, "uni_wo_qq_plots")

# Ensure Directories Exist
if (!dir.exists(maf_output_dir)) dir.create(maf_output_dir, recursive = TRUE)
if (!dir.exists(wo_output_dir)) dir.create(wo_output_dir, recursive = TRUE)

# Function to Load JSON Data into Lists (MAF or WO)
load_json_data <- function(directory_path) {
  subdirs <- list.dirs(directory_path, recursive = FALSE)
  maf_list <- list()
  wo_list <- list()
  
  for (subdir in subdirs) {
    identifier <- sub(".*_(\\d+).*", "\\1", basename(subdir))
    list_type <- ifelse(grepl("maf", subdir), "maf", "wo")
    
    json_files <- list.files(subdir, pattern = "\\.json$", full.names = TRUE)
    for (file in json_files) {
      data <- fromJSON(file)
      if (list_type == "maf") {
        maf_list[[paste0("stepping_", identifier, "_", basename(file))]] <- data
      } else {
        wo_list[[paste0("stepping_", identifier, "_", basename(file))]] <- data
      }
    }
  }
  return(list(maf_list = maf_list, wo_list = wo_list))
}

# Function to Generate Normally Distributed Effect Sizes
generate_effect_sizes <- function(n_loci) {
  return(runif(n_loci, min = -1, max = 1))
}

# Function to Apply Effect Sizes and Compute Phenotypes
apply_effect_sizes <- function(df, effect_sizes) {
  if (ncol(df) == 0) return(NULL)
  phenotype_df <- data.frame(individual = 1:nrow(df))
  
  for (loci_size in names(effect_sizes)) {
    num_loci <- length(effect_sizes[[loci_size]])
    num_selected_loci <- min(num_loci, ncol(df))
    selected_loci <- sample(colnames(df), num_selected_loci, replace = FALSE)
    
    phenotype <- as.numeric(as.matrix(df[, selected_loci, drop = FALSE]) %*% effect_sizes[[loci_size]][1:num_selected_loci])
    phenotype_df[[loci_size]] <- phenotype
  }
  
  return(phenotype_df)
}

# Function to Compute Qst
calculate_qst <- function(phenotype_df, num_pops = 8) {
  phenotype_df$population <- rep(1:num_pops, each = 20)
  qst_values <- list()
  
  for (loci_size in c("1_loci", "10_loci", "100_loci", "1000_loci")) {
    if (!loci_size %in% colnames(phenotype_df)) next
    
    V_W <- phenotype_df %>%
      group_by(population) %>%
      summarise(variance = var(!!sym(loci_size), na.rm = TRUE)) %>%
      summarise(V_W = mean(variance, na.rm = TRUE)) %>%
      pull(V_W)
    
    V_B <- var(phenotype_df %>%
                 group_by(population) %>%
                 summarise(mean_phenotype = mean(!!sym(loci_size), na.rm = TRUE)) %>%
                 pull(mean_phenotype), na.rm = TRUE)
    
    Qst <- V_B / (V_B + 2 * V_W)
    qst_values[[loci_size]] <- Qst
  }
  
  return(qst_values)
}

# Function to Process and Compute Qst for a Given Data List
process_data <- function(output_dir, lists) {
  results_list <- list()
  
  for (file_name in names(lists)) {
    json_data <- lists[[file_name]]
    stepping_id <- sub("stepping_([0-9]+).*", "\\1", file_name)
    
    for (replicate in names(json_data)) {
      loci_data <- json_data[[replicate]]
      replicate_df <- data.frame(individual = 1:160)
      
      for (loci_size in names(loci_data)) {
        loci_value_data <- loci_data[[loci_size]]$value
        df <- as.data.frame(do.call(cbind, loci_value_data))
        
        num_loci <- as.numeric(sub("_loci", "", loci_size))  # Extract number from "1_loci", "10_loci", etc.
        effect_sizes <- generate_effect_sizes(num_loci)
        
        phenotype <- as.numeric(as.matrix(df) %*% effect_sizes[1:ncol(df)])
        replicate_df[[loci_size]] <- phenotype
      }
      
      replicate_df$replicate <- replicate
      results_list[[paste0(stepping_id, "_", replicate)]] <- replicate_df
    }
  }
  
  # Organize Data by stepping
  unique_steppings <- unique(sub("_replicate.*", "", names(results_list)))
  stepping_replicates <- list()
  for (stepping in unique_steppings) {
    stepping_replicates[[stepping]] <- list()
    stepping_replicate_names <- names(results_list)[grepl(paste0("^", stepping, "_replicate"), names(results_list))]
    for (replicate_name in stepping_replicate_names) {
      replicate_id <- sub(paste0(stepping, "_"), "", replicate_name)
      stepping_replicates[[stepping]][[replicate_id]] <- results_list[[replicate_name]]
    }
  }
  
  # Calculate Qst for Each stepping and Replicate
  qst_results <- list()
  for (stepping in names(stepping_replicates)) {
    qst_results[[stepping]] <- list()
    for (replicate in names(stepping_replicates[[stepping]])) {
      phenotype_df <- stepping_replicates[[stepping]][[replicate]]
      qst_values <- calculate_qst(phenotype_df)
      qst_results[[stepping]][[replicate]] <- qst_values
    }
  }
  
  # Convert Results to DataFrame
  qst_df <- bind_rows(lapply(names(qst_results), function(stepping) {
    bind_rows(lapply(names(qst_results[[stepping]]), function(replicate) {
      data.frame(
        stepping = stepping,
        replicate = replicate,
        loci_size = names(qst_results[[stepping]][[replicate]]),
        Qst = unlist(qst_results[[stepping]][[replicate]])
      )
    }), .id = "replicate_id")
  }), .id = "stepping_id")
  
  qst_df$loci_size <- factor(qst_df$loci_size, levels = c("1_loci", "10_loci", "100_loci", "1000_loci"))
  qst_by_stepping <- split(qst_df, qst_df$stepping)
  
  # Generate QQ Plots
  loci_pairs <- list(
    c("1_loci", "10_loci"),
    c("1_loci", "100_loci"),
    c("1_loci", "1000_loci"),
    c("10_loci", "100_loci"),
    c("10_loci", "1000_loci"),
    c("100_loci", "1000_loci")
  )
  
  plot_qq_comparison <- function(qst_data, stepping_label) {
    stepping_output_dir <- file.path(output_dir, stepping_label)
    if (!dir.exists(stepping_output_dir)) {
      dir.create(stepping_output_dir, recursive = TRUE)
    }
    
    for (pair in loci_pairs) {
      loci_1 <- pair[1]
      loci_2 <- pair[2]
      
      qst_subset <- qst_data %>%
        filter(loci_size %in% c(loci_1, loci_2)) %>%
        select(loci_size, Qst)
      
      qst_wide <- qst_subset %>%
        pivot_wider(names_from = loci_size, values_from = Qst, values_fn = list) %>%
        unnest(cols = everything())
      
      if (nrow(qst_wide) < 10) next
      
      p <- ggplot(qst_wide, aes(sample = !!sym(loci_1))) +
        stat_qq(aes(sample = !!sym(loci_2))) +
        stat_qq_line(aes(sample = !!sym(loci_2))) +
        labs(title = paste("QQ Plot:", loci_1, "vs", loci_2, "- stepping", stepping_label),
             x = paste("Theoretical Quantiles of", loci_1),
             y = paste("Sample Quantiles of", loci_2)) +
        theme_minimal()
      
      plot_filename <- file.path(stepping_output_dir, paste0("QQ_", loci_1, "_vs_", loci_2, ".png"))
      ggsave(plot_filename, plot = p, width = 6, height = 6)
    }
  }
  
  plot_qq_comparison(qst_by_stepping[["005"]], "005")
  plot_qq_comparison(qst_by_stepping[["01"]], "01")
  plot_qq_comparison(qst_by_stepping[["015"]], "015")
  plot_qq_comparison(qst_by_stepping[["02"]], "02")
  
  
  # Perform Statistical Tests
  perform_stat_tests <- function(qst_data, stepping_label) {
    stepping_results <- list()
    
    for (pair in loci_pairs) {
      loci_1 <- pair[1]
      loci_2 <- pair[2]
      
      qst_1 <- qst_data %>% filter(loci_size == loci_1) %>% pull(Qst)
      qst_2 <- qst_data %>% filter(loci_size == loci_2) %>% pull(Qst)
      
      if (length(qst_1) < 10 | length(qst_2) < 10) {
        warning(paste("Not enough data for test:", loci_1, "vs", loci_2, "in stepping", stepping_label))
        next
      }
      
      ks_test <- ks.test(qst_1, qst_2)
      wilcox_test <- wilcox.test(qst_1, qst_2, alternative = "two.sided")
      
      stepping_results[[paste0(loci_1, "_vs_", loci_2)]] <- data.frame(
        stepping = stepping_label,
        loci_comparison = paste(loci_1, "vs", loci_2),
        ks_p_value = ks_test$p.value,
        wilcox_p_value = wilcox_test$p.value
      )
    }
    
    results_df <- bind_rows(stepping_results)
    stepping_output_dir <- file.path(output_dir, stepping_label)
    if (!dir.exists(stepping_output_dir)) {
      dir.create(stepping_output_dir, recursive = TRUE)
    }
    
    test_results_file <- file.path(stepping_output_dir, paste0("statistical_test_results_", stepping_label, ".csv"))
    write.csv(results_df, file = test_results_file, row.names = FALSE)
    
    return(results_df)
  }
  
  stat_results <- list()
  for (stepping in names(qst_by_stepping)) {
    stat_results[[stepping]] <- perform_stat_tests(qst_by_stepping[[stepping]], stepping)
  }
  
  final_stat_results <- bind_rows(stat_results)
  print(head(final_stat_results))
}

data_lists <- load_json_data(directory_path)
process_data(maf_output_dir, data_lists$maf_list)
process_data(wo_output_dir, data_lists$wo_list)
