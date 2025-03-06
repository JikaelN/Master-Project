library(ggplot2)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(kSamples)  # for Anderson–Darling test
library(conflicted)  # force conflicts to become errors

create_qqplot <- function(data1, data2, label1, label2, title_text) {
  # Sort the data (i.e. get empirical quantiles)
  qq_df <- data.frame(
    quantile1 = sort(data1),
    quantile2 = sort(data2)
  )
  p <- ggplot(qq_df, aes(x = quantile1, y = quantile2)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = title_text, x = label1, y = label2) +
    theme_minimal()
  return(p)
}

create_qqplot_migration <- function(data1, data2, label1, label2, title_text) {
  # Use the minimum length from both groups
  n <- min(length(data1), length(data2))
  qq_df <- data.frame(
    sample1 = sort(data1)[1:n],
    sample2 = sort(data2)[1:n]
  )
  
  p <- ggplot(qq_df, aes(x = sample1, y = sample2)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = title_text, x = label1, y = label2) +
    theme_minimal()
  
  return(p)
}
compare_qst_files <- function(input_dir, output_dir) {
  # Prefer dplyr's filter
  conflicts_prefer(dplyr::filter)
  
  list_df <- list()
  # List all CSV files in the input directory
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  for(file in files){
    # take base name and extract info from file_name
    base_name <- basename(file)
    parts <- unlist(strsplit(base_name, "_"))
    effect_size <-parts[1] #effect size distribution
    filtering_type <- parts[2] #maf or without maf
    type_model <- unlist(strsplit(parts[4], "\\."))[1] #island or stepping
    
    #save directory
    filter_dir <- file.path(output_dir, filtering_type)
    if (!dir.exists(filter_dir)) {
      dir.create(filter_dir, recursive = TRUE)
    }
    
    
    # Correction to have migration rate value as character
    temp <- read.csv(file, nrows = 1)
    num_cols <- ncol(temp)
    df <- read.csv(file,colClasses = c("character", rep(NA, num_cols - 1)))
    #change colname to be uniform across all type of file
    colnames(df) <- c("MR", "replicate", "loci_size", "QST")
    
    #split df per migration rate
    df_list <- split(df, df$MR)
    
    for(data in df_list){
      MR <- unique(data$MR)
      df_loci <- split(df, df$loci_size)
      pairs_to_compare <- list(c("1_loci","10_loci"), c("1_loci","100_loci"), c("1_loci","1000_loci"),
                               c("10_loci","100_loci"), c("10_loci","1000_loci"), c("100_loci","1000_loci"))
      
      # Loop over each pair and create a QQ plot.
      for(pair in pairs_to_compare) {
        l1 <- pair[1]
        l2 <- pair[2]
        
        # Subset QST values for each loci size.
        qst1 <- data %>% filter(loci_size == l1) %>% pull(QST)
        qst2 <- data %>% filter(loci_size == l2) %>% pull(QST)
        
        l1 <- unlist(strsplit(l1, "_"))[1]
        l2 <- unlist(strsplit(l2, "_"))[1]
        
        if(length(qst1) > 0 && length(qst2) > 0) {
          title_text <- paste(filtering_type, effect_size, l1, "vs", l2, "loci for", type_model, "model at", MR, "migration rate")
          p <- create_qqplot(qst1, qst2, paste(l1, "loci"), paste(l2, "loci"), title_text)
          
          # Save the QQ plot
          filename <- file.path(filter_dir, 
                                paste0(effect_size,"_QQplot_", l1, "_vs_", l2, "_", type_model, "_", MR, "_Migration rate.png"))
          ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 300)
          
          # ----- Perform ANOVA -----
          values <- c(qst1, qst2)
          group <- factor(c(rep("group1", length(qst1)), rep("group2", length(qst2))))
          anova_model <- aov(values ~ group)
          anova_summary <- summary(anova_model)
          anova_p <- anova_summary[[1]]$`Pr(>F)`[1]
          if(length(anova_p) == 0) anova_p <- NA
          
          # ----- Perform Anderson-Darling Test -----
          ad_result <- tryCatch({
            kSamples::ad.test(list(qst1, qst2))
            kSamples::ad.test(list(qst1, qst2))
          }, error = function(e) { list(p.value = NA) })
          ad_p <- ad_result$ad[1,3]
          ad_p
          if(length(ad_p) == 0) ad_p <- NA
          
          # ----- Perform Kolmogorov-Smirnov Test -----
          ks_result <- ks.test(qst1, qst2)
          ks_p <- ks_result$p.value
          if(length(ks_p) == 0) ks_p <- NA
          
          # Prepare a row for the CSV file for each test.
          loci_comp <- paste(l1, "vs", l2)
          
          anova_row <- data.frame(migration_rate = MR,effect_size = effect_size, loci_comparison = loci_comp, p_value = anova_p, stringsAsFactors = FALSE)
          ad_row    <- data.frame(migration_rate = MR,effect_size = effect_size, loci_comparison = loci_comp, p_value = ad_p, stringsAsFactors = FALSE)
          ks_row    <- data.frame(migration_rate = MR,effect_size = effect_size, loci_comparison = loci_comp, p_value = ks_p, stringsAsFactors = FALSE)
          
          # Define filenames for the CSV output (one file per test)
          anova_file <- file.path(filter_dir, "anova_results.csv")
          ad_file    <- file.path(filter_dir, "adtest_results.csv")
          ks_file    <- file.path(filter_dir, "kstest_results.csv")
          
          write_mode <- function(file, data) {
            if (!file.exists(file)) {
              write.table(data, file = file, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
            } else {
              write.table(data, file = file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
            }
          }
          
          write_mode(anova_file, anova_row)
          write_mode(ad_file, ad_row)
          write_mode(ks_file, ks_row)
          
        } else {
          message("Not enough data for loci sizes ", l1, " and ", l2)
        }
        
        
      }
    }
    
  }
}

compare_qst_files_migration <- function(input_dir, output_dir) {
  # Prefer dplyr's filter
  conflicts_prefer(dplyr::filter)
  
  # List all CSV files in the input directory
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for(file in files){
    # Take base name and extract info from file name
    base_name <- basename(file)
    parts <- unlist(strsplit(base_name, "_"))
    effect_size <- parts[1]                # effect size distribution (e.g., "inverse", "l", "linked", "uni", "norm")
    filtering_type <- parts[2]               # e.g., "maf" or "wo"
    type_model <- unlist(strsplit(parts[4], "\\."))[1]  # e.g., "island" or "stepping"
    
    # Define output directory: use output_dir/filtering_type
    filter_dir <- file.path(output_dir, filtering_type)
    if (!dir.exists(filter_dir)) {
      dir.create(filter_dir, recursive = TRUE)
    }
    
    # Read file forcing the first column as character
    temp <- read.csv(file, nrows = 1)
    num_cols <- ncol(temp)
    df <- read.csv(file, colClasses = c("character", rep(NA, num_cols - 1)))
    # Rename columns uniformly
    colnames(df) <- c("MR", "replicate", "loci_size", "QST")
    
    # ----- Compare the same loci between different migration rates -----
    # Split data by loci_size so that for each loci (e.g., "1_loci") we can compare across migration rates.
    loci_groups <- split(df, df$loci_size)
    
    for(loci in names(loci_groups)){
      data_loci <- loci_groups[[loci]]
      # Get unique migration rates present for this loci group
      mrs <- unique(data_loci$MR)
      if(length(mrs) >= 2){
        # Create all pairwise combinations of migration rates for this loci group
        migration_pairs <- combn(mrs, 2, simplify = FALSE)
        for(pair in migration_pairs){
          mr1 <- pair[1]
          mr2 <- pair[2]
          
          # Extract QST values for each migration rate in the same loci group.
          qst1 <- data_loci %>% filter(MR == mr1) %>% pull(QST)
          qst2 <- data_loci %>% filter(MR == mr2) %>% pull(QST)
          
          # For labeling, remove the "_loci" suffix (e.g., "1_loci" becomes "1")
          loci_label <- unlist(strsplit(loci, "_"))[1]
          
          if(length(qst1) > 0 && length(qst2) > 0){
            title_text <- paste(filtering_type, effect_size, loci_label, "loci for", type_model, "model comparing migration rates", mr1, "vs", mr2)
            p <- create_qqplot_migration(qst1, qst2, paste(mr1, "migration"), paste(mr2, "migration"), title_text)
            
            # Save the QQ plot
            filename <- file.path(filter_dir, 
                                  paste0(effect_size, "_QQplot_", loci_label, "_loci_", mr1, "_vs_", mr2, "_MigrationComparison_", type_model, ".png"))
            ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 300)
            
            # ----- Perform Statistical Tests -----
            # ANOVA: compare the two groups
            values <- c(qst1, qst2)
            group <- factor(c(rep("group1", length(qst1)), rep("group2", length(qst2))))
            anova_model <- aov(values ~ group)
            anova_summary <- summary(anova_model)
            anova_p <- anova_summary[[1]]$`Pr(>F)`[1]
            if(length(anova_p) == 0) anova_p <- NA
            
            # Anderson-Darling Test: extract the asymptotic p-value from the "ad" matrix.
            ad_result <- tryCatch({
              kSamples::ad.test(list(qst1, qst2))
            }, error = function(e) { list(ad = matrix(NA, nrow = 2, ncol = 3)) })
            ad_p <- ad_result$ad[1, 3]
            if(length(ad_p) == 0) ad_p <- NA
            
            # Kolmogorov-Smirnov Test:
            ks_result <- ks.test(qst1, qst2)
            ks_p <- ks_result$p.value
            if(length(ks_p) == 0) ks_p <- NA
            
            # Prepare rows for the CSV output. The columns are:
            # "Migration rate X v Y", "effect_size", "loci comparison", "p_value"
            migration_comp <- paste(mr1, "v", mr2)
            loci_comp <- paste(loci_label, "v", loci_label)
            
            anova_row <- data.frame(`Migration rate X v Y` = migration_comp,
                                    effect_size = effect_size,
                                    `loci comparison` = loci_comp,
                                    p_value = anova_p,
                                    stringsAsFactors = FALSE)
            ad_row <- data.frame(`Migration rate X v Y` = migration_comp,
                                 effect_size = effect_size,
                                 `loci comparison` = loci_comp,
                                 p_value = ad_p,
                                 stringsAsFactors = FALSE)
            ks_row <- data.frame(`Migration rate X v Y` = migration_comp,
                                 effect_size = effect_size,
                                 `loci comparison` = loci_comp,
                                 p_value = ks_p,
                                 stringsAsFactors = FALSE)
            
            # Define output filenames for the CSV results within filter_dir.
            anova_file <- file.path(filter_dir, "anova_results_migration_comparison.csv")
            ad_file <- file.path(filter_dir, "adtest_results_migration_comparison.csv")
            ks_file <- file.path(filter_dir, "kstest_results_migration_comparison.csv")
            
            write_mode <- function(file, data) {
              if (!file.exists(file)) {
                write.table(data, file = file, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
              } else {
                write.table(data, file = file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
              }
            }
            
            write_mode(anova_file, anova_row)
            write_mode(ad_file, ad_row)
            write_mode(ks_file, ks_row)
            
          } else {
            message("Not enough data for loci ", loci, " comparing migration rates ", mr1, " vs ", mr2)
          }
        }
      } else {
        message("Not enough migration rate groups for loci ", loci)
      }
    }
    
    # ----- Summary Statistics -----
    # Compute summary statistics (mean, standard deviation, and variance of QST)
    # for each loci_size for each migration rate.
    summary_stats <- df %>%
      group_by(MR, loci_size) %>%
      summarise(
        mean_QST = mean(as.numeric(QST), na.rm = TRUE),
        sd_QST = sd(as.numeric(QST), na.rm = TRUE),
        var_QST = var(as.numeric(QST), na.rm = TRUE),
        n = n()
      ) %>%
      ungroup() %>%
      mutate(effect_size = effect_size,
             filtering_type = filtering_type,
             type_model = type_model)
    
    summary_file <- file.path(filter_dir, paste0(effect_size, "_", type_model, "_", filtering_type, "_summary.csv"))
    write.csv(summary_stats, file = summary_file, row.names = FALSE)
  }
}

compare_qst_effect_size <- function(input_dir, output_dir){
  # Prefer dplyr's filter
  conflicts_prefer(dplyr::filter)
  
  #####TEST####
  temp_df <- list()
  
  # List all CSV files in the input directory
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  for(file in files){
    # Extract information from file name
    base_name <- basename(file)
    parts <- unlist(strsplit(base_name, "_"))
    effect_size <- parts[1]                
    filtering_type <- parts[2]               
    type_model <- unlist(strsplit(parts[4], "\\."))[1]  
    
    # Define output directory: use output_dir/filtering_type
    filter_dir <- file.path(output_dir, filtering_type)
    if (!dir.exists(filter_dir)) {
      dir.create(filter_dir, recursive = TRUE)
    }
    
    # Read file
    temp <- read.csv(file, nrows = 1)
    num_cols <- ncol(temp)
    df <- read.csv(file, colClasses = c("character", rep(NA, num_cols - 1)))
    colnames(df) <- c("MR", "replicate", "loci_size", "QST")
    temp_df[[filtering_type]][[effect_size]] <- df
  }
  
  # Get unique values of MR and loci_size
  unique_MR <- unique(temp_df$maf$uni$MR)
  unique_loci_size <- unique(temp_df$maf$uni$loci_size)
  
  # Loop over all unique MR and loci_size values
  for (MR_value in unique_MR) {
    for (loci_value in unique_loci_size) {
      
      # Define output paths
      maf_plot_path <- file.path(output_dir, paste0("maf_plot_MR_", MR_value, "_loci_", loci_value, ".png"))
      wo_plot_path <- file.path(output_dir, paste0("wo_plot_MR_", MR_value, "_loci_", loci_value, ".png"))
      stats_path <- file.path(output_dir, paste0("stats_MR_", MR_value, "_loci_", loci_value, ".txt"))
      
      # Open text file for writing statistical results
      sink(stats_path)
      
      # Extract data for "maf"
      maf_inverse <- temp_df$maf$inverse[temp_df$maf$inverse$MR == MR_value & temp_df$maf$inverse$loci_size == loci_value, ]
      maf_l <- temp_df$maf$l[temp_df$maf$l$MR == MR_value & temp_df$maf$l$loci_size == loci_value, ]
      maf_linked <- temp_df$maf$linked[temp_df$maf$linked$MR == MR_value & temp_df$maf$linked$loci_size == loci_value, ]
      maf_norm <- temp_df$maf$norm[temp_df$maf$norm$MR == MR_value & temp_df$maf$norm$loci_size == loci_value, ]
      maf_uni <- temp_df$maf$uni[temp_df$maf$uni$MR == MR_value & temp_df$maf$uni$loci_size == loci_value, ]
      
      # Extract data for "wo"
      wo_inverse <- temp_df$wo$inverse[temp_df$wo$inverse$MR == MR_value & temp_df$wo$inverse$loci_size == loci_value, ]
      wo_l <- temp_df$wo$l[temp_df$wo$l$MR == MR_value & temp_df$wo$l$loci_size == loci_value, ]
      wo_linked <- temp_df$wo$linked[temp_df$wo$linked$MR == MR_value & temp_df$wo$linked$loci_size == loci_value, ]
      wo_norm <- temp_df$wo$norm[temp_df$wo$norm$MR == MR_value & temp_df$wo$norm$loci_size == loci_value, ]
      wo_uni <- temp_df$wo$uni[temp_df$wo$uni$MR == MR_value & temp_df$wo$uni$loci_size == loci_value, ]
      
      # 📌 Compare maf only
      qst_values_maf <- c(maf_inverse$QST, maf_l$QST, maf_linked$QST, maf_norm$QST, maf_uni$QST)
      groups_maf <- rep(c("maf_inverse", "maf_l", "maf_linked", "maf_norm", "maf_uni"),
                        c(nrow(maf_inverse), nrow(maf_l), nrow(maf_linked), nrow(maf_norm), nrow(maf_uni)))
      
      if (length(qst_values_maf) > 0) {
        # Save maf boxplot
        png(maf_plot_path, width = 800, height = 600)
        boxplot(qst_values_maf ~ groups_maf, col = rainbow(5), 
                main = paste("maf - MR =", MR_value, ", Loci Size =", loci_value),
                ylab = "QST values", las = 1) 
        dev.off()
        
        # Perform statistical tests for maf
        anova_maf <- aov(qst_values_maf ~ groups_maf)
        print(summary(anova_maf))
        
        # Tukey post-hoc test (if ANOVA is significant)
        if (summary(anova_maf)[[1]][["Pr(>F)"]][1] < 0.05) {
          print(TukeyHSD(anova_maf))
        }
        
        # Kruskal-Wallis for maf
        kw_maf <- kruskal.test(qst_values_maf ~ groups_maf)
        print(kw_maf)
        
        # Pairwise Wilcoxon if Kruskal-Wallis is significant
        if (kw_maf$p.value < 0.05) {
          print(pairwise.wilcox.test(qst_values_maf, groups_maf, p.adjust.method = "bonferroni"))
        }
      }
      
      # 📌 Compare wo only
      qst_values_wo <- c(wo_inverse$QST, wo_l$QST, wo_linked$QST, wo_norm$QST, wo_uni$QST)
      groups_wo <- rep(c("wo_inverse", "wo_l", "wo_linked", "wo_norm", "wo_uni"),
                       c(nrow(wo_inverse), nrow(wo_l), nrow(wo_linked), nrow(wo_norm), nrow(wo_uni)))
      
      if (length(qst_values_wo) > 0) {
        # Save wo boxplot
        png(wo_plot_path, width = 800, height = 600)
        boxplot(qst_values_wo ~ groups_wo, col = rainbow(5), 
                main = paste("wo - MR =", MR_value, ", Loci Size =", loci_value),
                ylab = "QST values", las = 1)
        dev.off()
        
        # Perform statistical tests for wo
        anova_wo <- aov(qst_values_wo ~ groups_wo)
        print(summary(anova_wo))
        
        # Tukey post-hoc test (if ANOVA is significant)
        if (summary(anova_wo)[[1]][["Pr(>F)"]][1] < 0.05) {
          print(TukeyHSD(anova_wo))
        }
        
        # Kruskal-Wallis for wo
        kw_wo <- kruskal.test(qst_values_wo ~ groups_wo)
        print(kw_wo)
        
        # Pairwise Wilcoxon if Kruskal-Wallis is significant
        if (kw_wo$p.value < 0.05) {
          print(pairwise.wilcox.test(qst_values_wo, groups_wo, p.adjust.method = "bonferroni"))
        }
      }
      
      # Close text file output
      sink()
      
    } # End of loci_size loop
  } # End of MR loop
}









compare_qst_files("./data/qst_result/1_island", "stat_results/loci_distribution/1_island_between_loci_stat_result")
compare_qst_files("./data/qst_result/1_stepping", "stat_results/loci_distribution/1_stepping_between_loci_stat_result")
compare_qst_files("./data/qst_result/island", "stat_results/loci_distribution/island_between_loci_stat_result")
compare_qst_files("./data/qst_result/stepping", "stat_results/loci_distribution/stepping_between_loci_stat_result")

compare_qst_files_migration("./data/qst_result/1_island", "stat_results/migration/1_island_between_loci_stat_result")
compare_qst_files_migration("./data/qst_result/island", "stat_results/migration/island_between_loci_stat_result")
compare_qst_files_migration("./data/qst_result/1_stepping", "stat_results/migration/1_stepping_between_loci_stat_result")
compare_qst_files_migration("./data/qst_result/stepping", "stat_results/migration/stepping_between_loci_stat_result")

compare_qst_effect_size("./data/qst_result/1_island", "stat_results/effect_size/island")
compare_qst_effect_size("./data/qst_result/1_stepping", "stat_results/effect_size/stepping")

