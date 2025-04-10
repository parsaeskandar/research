library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)

# -----------------------------
# Chromosome mapping
# -----------------------------
chromosome_map <- c(
  "0" = "chr1",  "1" = "chr10", "2" = "chr11", "3" = "chr12", "4" = "chr13",
  "5" = "chr14", "6" = "chr15", "7" = "chr16", "8" = "chr17", "9" = "chr18",
  "10" = "chr19", "11" = "chr2", "12" = "chr20", "13" = "chr21", "14" = "chr22",
  "15" = "chr3", "16" = "chr4", "17" = "chr5", "18" = "chr6", "19" = "chr7",
  "20" = "chr8", "21" = "chr9", "22" = "chrM", "23" = "chrX", "24" = "chrY"
)

# -----------------------------
# Parse a single log file
# -----------------------------
parse_log_file <- function(file_path, version) {
  content <- read_file(file_path)
  chrom_match <- str_match(content, "(chr[\\w\\dMXYEBV]+|[0-9]+)\\.gbz")[,2]
  total_match <- as.numeric(str_match(content, "The size of the BWT is: (\\d+)")[,2])
  unique_match <- as.numeric(str_match(content, "covered by unique kmers is: (\\d+) / (\\d+) = ([0-9.]+)")[,2])
  extended_match <- as.numeric(str_match(content, "covered after extending the kmers is: (\\d+) / (\\d+) = ([0-9.]+)")[,2])

  if (!is.na(chrom_match) && !is.na(total_match) && !is.na(unique_match) && !is.na(extended_match)) {
    chrom_readable <- ifelse(chrom_match %in% names(chromosome_map), chromosome_map[[chrom_match]], chrom_match)
    return(data.frame(
      Chromosome = chrom_readable,
      Total = total_match,
      Unique = unique_match,
      Extended = extended_match,
      Version = version,
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}

# -----------------------------
# Load all logs from directory
# -----------------------------
read_all_logs <- function(base_dir, version) {
  files <- list.files(base_dir, pattern = "\\.log$", full.names = TRUE)
  log_data <- lapply(files, parse_log_file, version = version)
  do.call(rbind, log_data)
}

# -----------------------------
# Input dirs and parse
# -----------------------------
v1_logs <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/V1.1", "V1.1")
v2_logs <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/V2", "V2")
df <- bind_rows(v1_logs, v2_logs)

# -----------------------------
# Compute coverage breakdown
# -----------------------------
df <- df %>%
  mutate(
    Unique_Coverage = Unique / Total,
    Extension_Only = (Extended - Unique) / Total,
    Traversal_Only = 1 - (Extended / Total)
  )

# -----------------------------
# Long format
# -----------------------------
df_long <- df %>%
  select(Chromosome, Version, Unique_Coverage, Extension_Only, Traversal_Only) %>%
  pivot_longer(cols = starts_with("U"):starts_with("T"), names_to = "Stage", values_to = "Coverage") %>%
  mutate(Stage = factor(Stage, levels = c("Traversal_Only", "Extension_Only", "Unique_Coverage")))

# Natural order for chromosomes
chrom_order <- paste0("chr", c(1:22, "X", "Y", "M", "EBV"))
df_long$Chromosome <- factor(df_long$Chromosome, levels = chrom_order)

# -----------------------------
# Plot function for one version
# -----------------------------
plot_stacked_version <- function(df_long, version_label) {
  data <- df_long %>% filter(Version == version_label)

  stage_labels <- c(
    "Unique_Coverage" = "Unique k-mers",
    "Extension_Only" = "Extension",
    "Traversal_Only" = "Traversal"
  )

  stage_colors <- c(
    "Unique_Coverage" = "#4C78A8",
    "Extension_Only" = "#F58518",
    "Traversal_Only" = "#54A24B"
  )

  ggplot(data, aes(x = Chromosome, y = Coverage, fill = Stage)) +
    geom_bar(stat = "identity", width = 0.7, colour = "black", alpha = 0.95) +
    scale_fill_manual(values = stage_colors, labels = stage_labels, name = NULL) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 1, by= 0.2)) +
    labs(
      title = paste("Tag Array Coverage per Chromosome —", version_label),
      x = "Chromosome",
      y = "Fraction of BWT Covered"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 13),
      axis.title = element_text(face = "bold")
    )
}

# -----------------------------
# Render both plots (separately)
# -----------------------------
plot_stacked_version(df_long, "V1.1")
plot_stacked_version(df_long, "V2")
