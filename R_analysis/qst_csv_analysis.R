library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

stepping_vizualization <- function(directory){
  # List CSV files matching the pattern *_Qst_stepping.csv in the given directory
  files <- list.files(directory, pattern = "_Qst_stepping\\.csv$", full.names = TRUE)
  
  # Create an output directory for the plots (stepping_plot)
  output_dir <- file.path(directory, "stepping_plot")
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define migration labels mapping for facets (adjust keys if needed)
  migration_labels <- c("5" = "Migration rate: 5%",
                        "1" = "Migration rate: 10%",
                        "15" = "Migration rate: 15%",
                        "2" = "Migration rate: 20%")
  
  # Loop over each CSV file in the directory
  for(file in files){
    # Extract the prefix from the filename.
    # For example: "XX_maf_Qst_stepping.csv" -> "XX_maf"
    prefix <- sub("_Qst_stepping\\.csv$", "", basename(file))
    
    # Read the CSV file
    qst_df <- read.csv(file, stringsAsFactors = FALSE)
    
    # Convert columns to proper types:
    # Ensure 'loci_size' is a factor with the desired ordering
    qst_df$loci_size <- factor(qst_df$loci_size, levels = c("1_loci", "10_loci", "100_loci", "1000_loci"))
    # Convert stepping to a factor if it isn't already
    qst_df$stepping <- as.factor(qst_df$stepping)
    
    # Replace negative Qst values with 0
    qst_df_modified <- qst_df %>%
      mutate(Qst = ifelse(Qst < 0, 0, Qst))
    
    # Create the box plot faceted by stepping with custom facet labels
    box_plot <- ggplot(qst_df_modified, aes(x = loci_size, y = Qst, fill = loci_size)) +
      geom_boxplot() +
      facet_wrap(~ stepping, labeller = as_labeller(migration_labels)) +
      labs(title = "Box Plot of Qst by Loci Size for each Migration Rate",
           x = "Loci Size",
           y = "Qst") +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Save the box plot with a filename that includes the prefix
    box_filename <- file.path(output_dir, paste0(prefix, "_box_plot_stepping.png"))
    ggsave(filename = box_filename, plot = box_plot, width = 8, height = 6)
    
    # Create the density plot faceted by stepping
    density_plot <- ggplot(qst_df_modified, aes(x = Qst, color = loci_size)) +
      geom_density(size = 1) +
      facet_wrap(~ stepping, labeller = as_labeller(migration_labels)) +
      labs(title = "Density Plot of Qst by Loci Size for each Migration Rate",
           x = "Qst",
           y = "Density") +
      theme_minimal()
    
    # Save the density plot with a filename that includes the prefix
    density_filename <- file.path(output_dir, paste0(prefix, "_density_plot_stepping.png"))
    ggsave(filename = density_filename, plot = density_plot, width = 8, height = 6)
    
    cat("Plots saved for", prefix, "\n")
  }
}

stepping_vizualization_fst <- function(directory){
  # List CSV files matching *_Fst_stepping.csv in the specified directory
  files <- list.files(directory, pattern = "_Fst_stepping\\.csv$", full.names = TRUE)
  
  # Create an output directory for the plots (e.g., "stepping_plot")
  output_dir <- file.path(directory, "stepping_plot")
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define a mapping for the stepping values to migration rate labels.
  # Adjust the keys to match your stepping values.
  migration_labels <- c("5"  = "Migration rate: 5%",
                        "1"  = "Migration rate: 10%",
                        "15" = "Migration rate: 15%",
                        "2"  = "Migration rate: 20%")
  
  # Iterate over each CSV file
  for(file in files){
    # Extract the prefix from the filename.
    # E.g. "XX_maf_Fst_stepping.csv" -> "XX_maf"
    prefix <- sub("_Fst_stepping\\.csv$", "", basename(file))
    
    # Read the CSV file. 
    # If the first column (e.g., stepping) should be forced to character, adjust colClasses as needed.
    fst_df <- read.csv(file, stringsAsFactors = FALSE)
    
    # Convert stepping to factor (preserving any leading zeros if you have them).
    # Adjust the factor levels if you want a specific order.
    fst_df$stepping <- factor(fst_df$stepping, levels = c("5", "1", "15", "2"))
    
    # Set loci_size as an ordered factor
    fst_df$loci_size <- factor(fst_df$loci_size, levels = c("1_loci", "10_loci", "100_loci", "1000_loci"))
    
    # Replace negative Fst values with 0
    fst_df_modified <- fst_df %>%
      mutate(Fst = ifelse(Fst < 0, 0, Fst))
    
    ## Faceted Plots ##
    # 1) Box plot faceted by stepping
    box_plot <- ggplot(fst_df_modified, aes(x = loci_size, y = Fst, fill = loci_size)) +
      geom_boxplot() +
      facet_wrap(~ stepping, labeller = as_labeller(migration_labels)) +
      labs(title = "Box Plot of Fst by Loci Size for each Migration Rate",
           x = "Loci Size",
           y = "Fst") +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Save the box plot
    box_filename <- file.path(output_dir, paste0(prefix, "_box_plot_stepping.png"))
    ggsave(filename = box_filename, plot = box_plot, width = 8, height = 6)
    
    # 2) Density plot faceted by stepping
    density_plot <- ggplot(fst_df_modified, aes(x = Fst, color = loci_size)) +
      geom_density(size = 1) +
      facet_wrap(~ stepping, labeller = as_labeller(migration_labels)) +
      labs(title = "Density Plot of Fst by Loci Size for each Migration Rate\n(Negatives set to 0)",
           x = "Fst",
           y = "Density") +
      theme_minimal()
    
    density_filename <- file.path(output_dir, paste0(prefix, "_density_plot_stepping.png"))
    ggsave(filename = density_filename, plot = density_plot, width = 8, height = 6)
    
    # 3) Histogram faceted by stepping
    histogram_plot <- ggplot(fst_df_modified, aes(x = Fst, fill = loci_size)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "dodge") +
      facet_wrap(~ stepping, labeller = as_labeller(migration_labels)) +
      labs(title = "Histogram of Fst by Loci Size for each Migration Rate\n(Negatives set to 0)",
           x = "Fst",
           y = "Count") +
      theme_minimal()
    
    histogram_filename <- file.path(output_dir, paste0(prefix, "_histogram_plot_stepping.png"))
    ggsave(filename = histogram_filename, plot = histogram_plot, width = 8, height = 6)
    
    ## Overall Comparison Plot by Loci Size ##
    # Define comparisons for stepping groups (within each loci_size)
    # Adjust pairs to match your factor levels.
    comparisons <- list(
      c("5", "1"),
      c("5", "15"),
      c("5", "2"),
      c("1", "15"),
      c("1", "2"),
      c("15", "2")
    )
    
    # Compute max Fst to position annotations
    max_fst <- max(fst_df_modified$Fst, na.rm = TRUE)
    
    overall_box_plot_by_loci <- ggplot(fst_df_modified, aes(x = stepping, y = Fst, fill = stepping)) +
      geom_boxplot() +
      facet_wrap(~ loci_size) +
      # Kruskal-Wallis test
      stat_compare_means(
        method = "kruskal.test",
        label.x.npc = "right",
        label.y = max_fst * 1.1,  # place above the highest box
        size = 3
      ) +
      # Pairwise Wilcoxon comparisons
      stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        step.increase = 0.1,
        size = 3
      ) +
      labs(
        title = "Comparison of Fst by Stepping (Migration Rate) for each Loci Size",
        subtitle = "Pairwise comparisons (Wilcoxon tests) shown per loci size",
        x = "Stepping (Migration Rate)",
        y = "Fst (Negatives set to 0)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        # Increase overall text sizes for better readability
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
    
    overall_box_filename_by_loci <- file.path(output_dir, paste0(prefix, "_overall_box_plot_by_loci.png"))
    ggsave(filename = overall_box_filename_by_loci, plot = overall_box_plot_by_loci,
           width = 16, height = 8, dpi = 300)
    
    cat("Plots saved for", prefix, "\n")
  }
}


island_vizualization <- function(directory){
  # List CSV files matching the pattern *_Qst_island.csv in the given directory
  files <- list.files(directory, pattern = "_Qst_island\\.csv$", full.names = TRUE)
  
  # Create an output directory for the plots (island_plot)
  output_dir <- file.path(directory, "island_plot")
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define migration labels mapping for facets (adjust keys if needed)
  migration_labels <- c("001" = "Migration rate: 1%",
                        "005" = "Migration rate: 5%",
                        "008" = "Migration rate: 8%",
                        "01" = "Migration rate: 10%")
  
  # Loop over each CSV file in the directory
  for(file in files){
    # Extract the prefix from the filename.
    # For example: "XX_maf_Qst_island.csv" -> "XX_maf"
    prefix <- sub("_Qst_island\\.csv$", "", basename(file))
    
    # Read the CSV file
    qst_df <- read.csv(file, stringsAsFactors = FALSE ,colClasses = c("character", NA, NA, NA))
    
    
    # Convert columns to proper types:
    # Ensure 'loci_size' is a factor with the desired ordering
    qst_df$loci_size <- factor(qst_df$loci_size, levels = c("1_loci", "10_loci", "100_loci", "1000_loci"))
    # Convert island to a factor if it isn't already
    qst_df$island <- as.factor(qst_df$island)
    
    # Replace negative Qst values with 0
    qst_df_modified <- qst_df %>%
      mutate(Qst = ifelse(Qst < 0, 0, Qst))
    
    # Create the box plot faceted by island with custom facet labels
    box_plot <- ggplot(qst_df_modified, aes(x = loci_size, y = Qst, fill = loci_size)) +
      geom_boxplot() +
      facet_wrap(~ island, labeller = as_labeller(migration_labels)) +
      labs(title = "Box Plot of Qst by Loci Size for each Migration Rate",
           x = "Loci Size",
           y = "Qst") +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Save the box plot with a filename that includes the prefix
    box_filename <- file.path(output_dir, paste0(prefix, "_box_plot_island.png"))
    ggsave(filename = box_filename, plot = box_plot, width = 8, height = 6)
    
    # Create the density plot faceted by island
    density_plot <- ggplot(qst_df_modified, aes(x = Qst, color = loci_size)) +
      geom_density(size = 1) +
      facet_wrap(~ island, labeller = as_labeller(migration_labels)) +
      labs(title = "Density Plot of Qst by Loci Size for each Migration Rate",
           x = "Qst",
           y = "Density") +
      theme_minimal()
    
    # Save the density plot with a filename that includes the prefix
    density_filename <- file.path(output_dir, paste0(prefix, "_density_plot_island.png"))
    ggsave(filename = density_filename, plot = density_plot, width = 8, height = 6)
    
    cat("Plots saved for", prefix, "\n")
  }
}

island_vizualization_fst <- function(directory){
  # List CSV files matching *_Fst_island.csv in the specified directory
  files <- list.files(directory, pattern = "_Fst_island\\.csv$", full.names = TRUE)
  
  # Create an output directory for the plots (e.g., "island_plot")
  output_dir <- file.path(directory, "island_plot")
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define a mapping for the island values to custom labels.
  # Adjust keys and labels to your desired presentation.
  island_labels <- c("001" = "Migration rate: 1%",
                     "005" = "Migration rate: 5%",
                     "008" = "Migration rate: 8%",
                     "01"  = "Migration rate: 10%")
  
  # Iterate over each CSV file
  for(file in files){
    # Extract the prefix from the filename.
    # E.g. "XX_maf_Fst_island.csv" -> "XX_maf"
    prefix <- sub("_Fst_island\\.csv$", "", basename(file))
    
    # Read the CSV file. Force the first column (island) to character.
    fst_df <- read.csv(file, stringsAsFactors = FALSE,
                       colClasses = c("character", rep(NA, ncol(read.csv(file, nrows = 1))-1)))
    
    # Convert the island column to a factor with a specific order.
    fst_df$island <- factor(fst_df$island, levels = c("005", "001", "008", "01"))
    
    # Set loci_size as an ordered factor
    fst_df$loci_size <- factor(fst_df$loci_size, levels = c("1_loci", "10_loci", "100_loci", "1000_loci"))
    
    # Replace negative Fst values with 0
    fst_df_modified <- fst_df %>%
      mutate(Fst = ifelse(Fst < 0, 0, Fst))
    
    ## Faceted Plots ##
    # Box plot faceted by island
    box_plot <- ggplot(fst_df_modified, aes(x = loci_size, y = Fst, fill = loci_size)) +
      geom_boxplot() +
      facet_wrap(~ island, labeller = as_labeller(island_labels)) +
      labs(title = "Box Plot of Fst by Loci Size for each Migration rate",
           x = "Loci Size",
           y = "Fst") +
      theme_minimal() +
      theme(legend.position = "none")
    
    box_filename <- file.path(output_dir, paste0(prefix, "_box_plot_island.png"))
    ggsave(filename = box_filename, plot = box_plot, width = 8, height = 6)
    
    # Density plot faceted by island
    density_plot <- ggplot(fst_df_modified, aes(x = Fst, color = loci_size)) +
      geom_density(size = 1) +
      facet_wrap(~ island, labeller = as_labeller(island_labels)) +
      labs(title = "Density Plot of Fst by Loci Size for each Migration rate\n(Negatives set to 0)",
           x = "Fst",
           y = "Density") +
      theme_minimal()
    
    density_filename <- file.path(output_dir, paste0(prefix, "_density_plot_island.png"))
    ggsave(filename = density_filename, plot = density_plot, width = 8, height = 6)
    
    # Histogram faceted by island
    histogram_plot <- ggplot(fst_df_modified, aes(x = Fst, fill = loci_size)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "dodge") +
      facet_wrap(~ island, labeller = as_labeller(island_labels)) +
      labs(title = "Histogram of Fst by Loci Size for each Migration rate\n(Negatives set to 0)",
           x = "Fst",
           y = "Count") +
      theme_minimal()
    
    histogram_filename <- file.path(output_dir, paste0(prefix, "_histogram_plot_island.png"))
    ggsave(filename = histogram_filename, plot = histogram_plot, width = 8, height = 6)
    
    ## Overall Comparison Plot by Loci Size ##
    # Define comparisons for island groups (within each loci_size)
    comparisons <- list(
      c("005", "001"),
      c("005", "008"),
      c("005", "01"),
      c("001", "008"),
      c("001", "01"),
      c("008", "01")
    )
    
    overall_box_plot_by_loci <- ggplot(fst_df_modified, aes(x = island, y = Fst, fill = island)) +
      geom_boxplot() +
      facet_wrap(~ loci_size) +
      # Kruskal-Wallis test annotation at the top left, smaller text size
      stat_compare_means(
        method = "kruskal.test",
        label.x.npc = "right",
        label.y = max(fst_df_modified$Fst, na.rm = TRUE) * 1.1,
        size = 3  # Reduce text size for the KW test annotation
      ) +
      # Pairwise Wilcoxon comparisons
      stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        step.increase = 0.1,  # or manually specify label.y values
        size = 3
      ) +
      labs(
        title = "Comparison of Fst by migration rate for each Loci Size",
        subtitle = "Pairwise comparisons (Wilcoxon tests) shown per loci size",
        x = "Island (Migration Rate)",
        y = "Fst (Negatives set to 0)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        # Increase overall text sizes for better readability
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
    
    
    
    overall_box_filename_by_loci <- file.path(output_dir, paste0(prefix, "_overall_box_plot_by_loci.png"))
    ggsave(filename = overall_box_filename_by_loci, plot = overall_box_plot_by_loci,
           width = 17, height = 8, dpi = 300)
    
    cat("Plots saved for", prefix, "\n")
  }
}

# Call the function using the directory containing your CSV files
stepping_vizualization("./data/qst_result/stepping")
island_vizualization("./data/qst_result/island")

stepping_vizualization_fst("./data/fst_result/stepping")
island_vizualization_fst("./data/fst_result/island")
