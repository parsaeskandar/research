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
parse_log_compression <- function(file_path, graph) {
  content <- read_lines(file_path)

  # Chromosome name: try several patterns found in logs
  # Examples:
  #  - "running for chrM in .../chrM"
  #  - "running for 24"
  chrom <- NA
  run_line <- str_subset(content, "^running for ")
  if (length(run_line)) {
    # Prefer explicit 'chr' name if present
    chr_name <- str_match(run_line[1], "running for (chr[\\w\\dMXYEBV]+)")
    if (!is.na(chr_name[1])) chrom <- chr_name[2]
    if (is.na(chrom)) {
      # Fallback to numeric code (e.g., 24) and map to canonical chromosome name
      chr_num <- str_match(run_line[1], "running for ([0-9]+)")
      if (!is.na(chr_num[1])) {
        key <- chr_num[2]
        if (key %in% names(chromosome_map)) chrom <- chromosome_map[[key]]
      }
    }
  }
  # Final fallback: try to extract from a path containing '/chrX' if present
  if (is.na(chrom)) {
    path_chr <- str_match(paste(content, collapse = "\n"), "/(chr[\\w\\dMXYEBV]+)/?")
    if (!is.na(path_chr[1])) chrom <- path_chr[2]
  }

  # Number of bases in the tag array: use the denominator from the most refined coverage line available
  # Prefer "after extending the kmers"; fallback to "by unique kmers"
  base_line <- str_subset(content, "after extending the kmers is: ")
  if (!length(base_line)) base_line <- str_subset(content, "by unique kmers is: ")
  num_bases <- if (length(base_line)) as.numeric(str_match(base_line[1], "/ (\\d+) =")[,2]) else NA

  # Number of items in the BPlusTree (Tag Runs): handle both 'final' and generic messages
  bpt_lines <- str_subset(content, "The (final )?number of items in the BPlusTree is")
  num_tag_runs <- if (length(bpt_lines)) as.numeric(str_match(bpt_lines[1], "(\\d+)$")[,1]) else NA

  # Number of logical runs in the grlBWT (BWT Runs)
  grlbwt_line <- str_subset(content, "FastLocate::FastLocate\\(\\): [0-9]+ logical runs in the grlBWT")
  num_bwt_runs <- if (length(grlbwt_line)) as.numeric(str_match(grlbwt_line[1], "([0-9]+) logical runs")[,2]) else NA

  if (!is.na(chrom) && !is.na(num_bases) && !is.na(num_tag_runs) && !is.na(num_bwt_runs)) {
    return(data.frame(
      Chromosome = chrom,
      graph = graph,
      Bases = num_bases,
      Tag_Runs = num_tag_runs,
      BWT_Runs = num_bwt_runs,
      Compression_Ratio = num_bases / num_tag_runs 
    ))
  } else {
    return(NULL)
  }
}

# -----------------------------
# Parse all logs in a directory
# -----------------------------
read_compression_logs <- function(base_dir, graph) {
  files <- list.files(base_dir, pattern = "\\.log$", full.names = TRUE)
  logs <- lapply(files, parse_log_compression, graph = graph)
  do.call(rbind, logs)
}

# -----------------------------
# Load logs (set your paths here)
# -----------------------------
v1_logs <- read_compression_logs("/Users/seeskand/PycharmProjects/python_plots/data/kmer_extension_08-25/v1.1", "V1.1")
v2_logs <- read_compression_logs("/Users/seeskand/PycharmProjects/python_plots/data/kmer_extension_08-25/v2", "V2")
df_compression <- bind_rows(v1_logs, v2_logs)

# -----------------------------
# Override/augment Tag Runs with provided values
# -----------------------------

# V1.1 provided Tag Runs keyed by numeric IDs, map to chromosome names
manual_v11_raw <- c(
  "4" = 344191611,
  "18" = 499530293,
  "12" = 290116066,
  "20" = 560633493,
  "7" = 499415185,
  "19" = 699044282,
  "21" = 567407001,
  "15" = 612385650,
  "13" = 203269708,
  "24" = 212854926,
  "10" = 335331466,
  "0" = 1028847222,
  "23" = 632128932,
  "9" = 254419843,
  "5" = 430287496,
  "22" = 40608,
  "17" = 723972655,
  "1" = 468343884,
  "8" = 427434282,
  "2" = 415301595,
  "11" = 888834855,
  "6" = 510525049,
  "3" = 376710868,
  "16" = 641273885,
  "14" = 311209232
)

manual_v11_chr <- unname(vapply(names(manual_v11_raw), function(k) {
  if (k %in% names(chromosome_map)) chromosome_map[[k]] else NA_character_
}, character(1)))
manual_v11_df <- data.frame(
  Chromosome = manual_v11_chr,
  graph = "V1.1",
  Tag_Runs = as.numeric(unname(manual_v11_raw)),
  stringsAsFactors = FALSE
)
manual_v11_df <- manual_v11_df[!is.na(manual_v11_df$Chromosome), , drop = FALSE]

# V2 provided Tag Runs keyed by chromosome names
manual_v2_df <- data.frame(
  Chromosome = c(
    "chr8","chr18","chr7","chr6","chr5","chr17","chrM","chr11","chr9","chr22",
    "chr16","chr1","chr20","chrY","chr21","chr13","chr14","chr4","chr10","chrX",
    "chr2","chr3","chr15","chr12","chr19"
  ),
  graph = "V2",
  Tag_Runs = c(
    1078055627,390493821,1605454409,770820485,1452800237,945447774,61722,693326083,1024825353,971253650,
    1201908691,2398096328,549133236,833501557,366057594,624072183,971997604,955119062,722900495,1290320302,
    1888553087,933649765,1242373705,558265008,697696026
  ),
  stringsAsFactors = FALSE
)

manual_tag_runs_df <- bind_rows(manual_v11_df, manual_v2_df)

# Merge manual Tag_Runs, overriding parsed Tag_Runs when provided
if (!is.null(df_compression) && nrow(df_compression) > 0) {
  df_compression <- df_compression %>%
    select(Chromosome, graph, Bases, Tag_Runs, BWT_Runs) %>%
    full_join(manual_tag_runs_df, by = c("Chromosome", "graph"), suffix = c("", "_manual")) %>%
    mutate(Tag_Runs = coalesce(Tag_Runs_manual, Tag_Runs)) %>%
    select(Chromosome, graph, Bases, Tag_Runs, BWT_Runs)
} else {
  # Build a minimal df from manual Tag_Runs if parsing yielded nothing
  df_compression <- manual_tag_runs_df %>%
    mutate(Bases = NA_real_, BWT_Runs = NA_real_) %>%
    select(Chromosome, graph, Bases, Tag_Runs, BWT_Runs)
}

# -----------------------------
# Plotting
# -----------------------------

# Natural chromosome order
chrom_order <- paste0("chr", c(1:22, "X", "Y", "M"))
if (!is.null(df_compression) && nrow(df_compression) > 0) {
  df_compression$Chromosome <- factor(df_compression$Chromosome, levels = chrom_order)
} else {
  stop("No parsed data found from logs. Please check parsing patterns.")
}

print(df_compression)

# Color palette
graph_colors <- c("V1.1" = "#4C78A8", "V2" = "#F58518")

# --- Plot 1: Compression Ratio (Bases per Tag Run)
# ggplot(df_compression, aes(x = Chromosome, y = Compression_Ratio, fill = graph)) +
#   geom_point(
#     aes(fill = graph, shape = graph),
#     size = 3.8,
#     color = "black",
#     position = position_identity()
#   ) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 250), breaks = seq(0, 250, by= 25)) +
#   scale_shape_manual(name = "Graph", values = c("V1.1" = 21, "V2" = 21)) +
#   scale_fill_manual(name = "Graph", values = graph_colors) +
#   labs(
#     title = "Compression Efficiency of Tag Arrays",
#     subtitle = "Bases per tag run across chromosomes in HPRC v1.1 and v2 graphs",
#     x = "Chromosome",
#     y = "Compression Ratio (Bases per Tag Run)"
#   ) +
#   theme_minimal(base_size = 15) +
#   theme(
#     plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
#     plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#     axis.text.y = element_text(size = 12),
#     axis.title.y = element_text(face = "bold", size = 14),
#     axis.title.x = element_text(face = "bold", size = 14),
#     panel.grid.major.x = element_blank(),
#     panel.grid.major.y = element_line(color = "gray85", linetype = "dashed"),
#     legend.position = "top",
#     legend.title = element_text(size = 13, face = "bold"),
#     legend.text = element_text(size = 12)
#   )

# --- Plot 2: V2 / V1.1 Ratios for Bases, Tag Runs, and BWT Runs

df_wide <- df_compression %>%
  select(Chromosome, graph, Bases, Tag_Runs, BWT_Runs) %>%
  pivot_wider(names_from = graph, values_from = c(Bases, Tag_Runs, BWT_Runs)) %>%
  mutate(
    base_ratio = Bases_V2 / Bases_V1.1,
    tag_ratio  = Tag_Runs_V2 / Tag_Runs_V1.1,
    bwt_ratio  = BWT_Runs_V2 / BWT_Runs_V1.1
  )

df_ratio_long <- df_wide %>%
  select(Chromosome, base_ratio, tag_ratio, bwt_ratio) %>%
  pivot_longer(cols = c(base_ratio, tag_ratio, bwt_ratio), names_to = "Metric", values_to = "Ratio") %>%
  mutate(
    Metric = recode(Metric,
                    "base_ratio" = "Bases",
                    "tag_ratio" = "Tag Runs",
                    "bwt_ratio"  = "BWT Runs")
  )


df_ratio_long$Chromosome <- factor(df_ratio_long$Chromosome, levels = chrom_order)

print(df_ratio_long)
ggplot(df_ratio_long, aes(x = Chromosome, y = Ratio, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  scale_fill_manual(values = c("Bases" = "#4C78A8", "Tag Runs" = "#F58518", "BWT Runs" = "#54A24B")) +
  labs(
    title = "Ratio of Bases and Run Counts between V1.1 and V2 graphs",
    x = "Chromosome",
    y = "Ratio (V2 / V1.1)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.title.x = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linetype = "dashed"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )
