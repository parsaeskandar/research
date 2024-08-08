.libPaths(c("/opt/homebrew/lib/R/4.3/site-library", .libPaths()))

# Install and load required packages
# install.packages("ggplot2")
# install.packages("tidyr")
# install.packages("dplyr")
library(ggplot2)
library(tidyr)
library(dplyr)

# Load the data
data <- read.csv("output_incomplete.csv")
names(data) <- gsub("^\\s+|\\s+$", "", names(data))

# Function to convert chromosome labels to numeric-like values
convert_chr <- function(chr) {
  if (chr == "chrX") {
    return(23)
  } else if (chr == "chrY") {
    return(24)
  } else {
    return(as.numeric(gsub("chr", "", chr)))
  }
}

# Apply the function and create a new factor for Chromosome
data$Chromosome_numeric <- sapply(data$Chromosome, convert_chr)
data$Chromosome <- factor(data$Chromosome, levels = data$Chromosome[order(data$Chromosome_numeric)])
# gsub("chr", "", data$Chromosome)

# Reshape the data for ggplot2
long_data <- data %>%
  select(Chromosome, Tag_arrays_unique_kmers, Tag_arrays_after_extension) %>%
  gather(key = "Stage", value = "Coverage", -Chromosome)

# Mean coverage of stage T2 for all chromosomes
mean(long_data$Coverage[long_data$Stage == "Tag_arrays_after_extension"])

# Plot the data
ggplot(long_data, aes(y = Chromosome, x = Coverage, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  scale_fill_brewer(palette = "Paired") +
  coord_cartesian(xlim = c(0.6, 1)) +
  theme_minimal() +
  labs(x = "Coverage", y = "Chromosome", fill = "Stage") +
  ggtitle("Tag Arrays Coverage for Each Chromosome")
