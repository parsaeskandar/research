library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(ggpattern)
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
# Time parser (h:mm:ss or mm:ss)
# -----------------------------
parse_time <- function(line) {
  if (str_detect(line, ":\\d+:\\d+")) {
    parts <- as.numeric(unlist(str_split(str_extract(line, "\\d+:\\d+:\\d+"), ":")))
    return(parts[1] + parts[2]/60 + parts[3]/3600)
  } else if (str_detect(line, "\\d+:\\d+")) {
    parts <- as.numeric(unlist(str_split(str_extract(line, "\\d+:\\d+"), ":")))
    return(parts[1]/60 + parts[2]/3600)
  } else {
    return(NA)
  }
}

# -----------------------------
# Parse a single log file
# -----------------------------
parse_log <- function(file_path, graph) {
  content <- read_lines(file_path)

  # Chromosome name
  chrom_match <- str_subset(content, "\\.gbz")
  chrom <- if (length(chrom_match)) str_match(chrom_match[1], "(chr[\\w\\dMXYEBV]+|[0-9]+)\\.gbz")[,2] else NA
  if (!is.na(chrom) && chrom %in% names(chromosome_map)) chrom <- chromosome_map[[chrom]]

  # Tag array time
  time_line <- str_subset(content, "Elapsed \\(wall clock\\) time")
  tag_time_hr <- if (length(time_line)) parse_time(time_line[1]) else NA

  # Max RSS
  mem_line <- str_subset(content, "Maximum resident set size")
  mem_gb <- if (length(mem_line)) as.numeric(str_match(mem_line[1], "(\\d+)")[,2]) / 1e6 else NA

  # grlBWT time
  grlbwt_time_hr <- parse_grlbwt_time(content)

  if (!is.na(chrom) && !is.na(tag_time_hr) && !is.na(mem_gb)) {
    return(data.frame(
      Chromosome = chrom,
      Tag_Time_Hours = tag_time_hr + grlbwt_time_hr/2,
      grlBWT_Time_Hours = grlbwt_time_hr / 2,
      Max_Memory_GB = mem_gb,
      graph = graph
    ))
  } else {
    return(NULL)
  }
}


# -----------------------------
# Parse all grlBWT elapsed times
# -----------------------------
parse_grlbwt_time <- function(content) {
  time_lines <- str_subset(content, "Elapsed time")

  time_seconds <- sapply(time_lines, function(line) {
    if (str_detect(line, "\\(mm:ss\\.ms\\):")) {
      time_str <- str_match(line, "(\\d+):(\\d+\\.\\d+)")[,2:3]
      if (!any(is.na(time_str))) {
        return(as.numeric(time_str[1]) * 60 + as.numeric(time_str[2]))
      }
    } else if (str_detect(line, "\\(ss\\.ms\\):")) {
      time_str <- str_match(line, "(\\d+\\.\\d+)")[,2]
      if (!is.na(time_str)) {
        return(as.numeric(time_str))
      }
    } else if (str_detect(line, "\\(ms\\):")) {
      time_str <- str_match(line, "(\\d+)")[,2]
      if (!is.na(time_str)) {
        return(as.numeric(time_str) / 1000)
      }
    }
    return(0)  # Default fallback
  })

  return(sum(unlist(time_seconds), na.rm = TRUE) / 3600)  # Convert seconds to hours
}


# -----------------------------
# Batch parse logs from a directory
# -----------------------------
read_all_logs <- function(base_dir, graph) {
  files <- list.files(base_dir, pattern = "\\.log$", full.names = TRUE)
  logs <- lapply(files, parse_log, graph = graph)
  do.call(rbind, logs)
}

# -----------------------------
# Load your logs here (change these paths)
# -----------------------------
v1_logs <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/kmer_extension_08-25/v1.1", "V1.1")
v2_logs <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/kmer_extension_08-25/v2", "V2")
df_res <- bind_rows(v1_logs, v2_logs)

# Sort chromosomes naturally
chrom_order <- paste0("chr", c(1:22, "X", "Y", "M", "EBV"))
# df_res$Chromosome <- factor(df_res$Chromosome, levels = chrom_order)

# Clean color palette (color-blind friendly)
graph_colors <- c(
  "V1.1" = "#1f78b4",  # A bold blue
  "V2"   = "#e31a1c"   # A strong red
)

# Reorder chromosomes

## Memory usage comparison (side-by-side, not stacked)
# Ensure ordering and factor levels
df_res$Chromosome <- factor(df_res$Chromosome, levels = chrom_order)
df_res$graph <- factor(df_res$graph, levels = c("V1.1", "V2"))

ggplot(df_res, aes(x = Chromosome, y = Max_Memory_GB, fill = graph)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "black") +
  scale_fill_manual(
    name = "Graph",
    values = c("V1.1" = "#1f77b4",  # blue
               "V2"   = "#ff7f0e")  # orange
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Max RSS Memory Usage per Chromosome",
    x = "Chromosome",
    y = "Max Memory (GB)"
  ) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16)
  )

df_res %>%
  group_by(graph) %>%
  summarise(
    Max_grlBWT_Time_Hours = max(grlBWT_Time_Hours, na.rm = TRUE),
    Max_Tag_Time_Hours = max(Tag_Time_Hours, na.rm = TRUE),
    total_time = sum(Tag_Time_Hours, na.rm = TRUE),
  ) %>%
  print()     


# -----------------------------
# Compare old vs new memory usage (paper_output vs kmer_extension_08-25)
# -----------------------------

# Load OLD (paper_output) logs
old_v1 <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/paper_output/v1.1", "V1.1")
old_v2 <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/paper_output/v2", "V2")
df_old <- bind_rows(old_v1, old_v2) %>% mutate(dataset = "Old")

# Load NEW (kmer_extension_08-25) logs
new_v1 <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/kmer_extension_08-25/v1.1", "V1.1")
new_v2 <- read_all_logs("/Users/seeskand/PycharmProjects/python_plots/data/kmer_extension_08-25/v2", "V2")
df_new <- bind_rows(new_v1, new_v2) %>% mutate(dataset = "New")

# Combine and prepare ordering
df_all <- bind_rows(df_old, df_new) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = chrom_order),
    graph = factor(graph, levels = c("V1.1", "V2")),
    dataset = factor(dataset, levels = c("Old", "New"))
  )

## Overlay comparison: Old (faded) vs New (solid) per Chromosome and Graph
# Average duplicates to one value per Chromosome/graph/dataset
df_all_mem_avg <- df_all %>%
  select(Chromosome, graph, dataset, Max_Memory_GB) %>%
  group_by(Chromosome, graph, dataset) %>%
  summarise(Max_Memory_GB = mean(Max_Memory_GB, na.rm = TRUE), .groups = "drop")

# Shared position so bars for V1.1 and V2 align; Old shown wide and faded behind New
pos <- position_dodge(width = 0.8)

ggplot(df_all_mem_avg, aes(x = Chromosome, y = Max_Memory_GB, fill = graph)) +
  geom_col(
    data = subset(df_all_mem_avg, dataset == "Old"),
    aes(alpha = dataset),
    position = pos,
    width = 0.7,
    colour = "black"
  ) +
  geom_col(
    data = subset(df_all_mem_avg, dataset == "New"),
    aes(alpha = dataset),
    position = pos,
    width = 0.45,
    colour = "black"
  ) +
  scale_fill_manual(values = c("V1.1" = "#1f77b4", "V2" = "#ff7f0e"), name = "Graph") +
  scale_alpha_manual(values = c("Old" = 0.35, "New" = 1.0), name = "Dataset",
                     labels = c("Old" = "Old (paper_output)", "New" = "New (kmer_extension_08-25)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Max RSS Memory per Chromosome: Old (faded) vs New (solid)",
    x = "Chromosome",
    y = "Max Memory (GB)"
  ) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16)
  )

# Print average absolute reduction (GB) per graph and overall
df_reduction <- df_all_mem_avg %>%
  select(Chromosome, graph, dataset, Max_Memory_GB) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = Max_Memory_GB) %>%
  filter(!is.na(Old) & !is.na(New)) %>%
  mutate(Abs_Reduction_GB = Old - New)

df_reduction %>%
  group_by(graph) %>%
  summarise(avg_abs_reduction_gb = round(mean(Abs_Reduction_GB, na.rm = TRUE), 3)) %>%
  print()

df_reduction %>%
  summarise(overall_avg_abs_reduction_gb = round(mean(Abs_Reduction_GB, na.rm = TRUE), 3)) %>%
  print()

# Print average percent reduction per graph and overall
df_reduction_pct <- df_reduction %>%
  mutate(Pct_Reduction = 100 * (Old - New) / Old)

df_reduction_pct %>%
  group_by(graph) %>%
  summarise(avg_pct_reduction = round(mean(Pct_Reduction, na.rm = TRUE), 2)) %>%
  print()

df_reduction_pct %>%
  summarise(overall_avg_pct_reduction = round(mean(Pct_Reduction, na.rm = TRUE), 2)) %>%
  print()