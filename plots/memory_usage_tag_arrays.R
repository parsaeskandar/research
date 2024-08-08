.libPaths(c("/opt/homebrew/lib/R/4.3/site-library", .libPaths()))




library(ggplot2)
library(tidyr)
library(dplyr)

# Load the data
data <- read.csv("output_incomplete.csv")
names(data) <- gsub("^\\s+|\\s+$", "", names(data))
data$Mem <- as.numeric(sub("MiB", "", data$Mem)) / 1024 # Convert GB to MB
# data <- data %>%
#   mutate(
#     Max_Mem = Max_Mem / BWT
#   )
data
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

# Plot for Max Memory
ggplot(data, aes(y = Chromosome, x = Mem, fill = Chromosome)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = round(Mem, 1)), hjust = 1.1, color = "black") +
  scale_x_continuous(labels = scales::comma) + # Use comma for large numbers
  theme_minimal() +
  labs(x = "Memory (GiB)", y = "Chromosome") +
  ggtitle("Max Memory Usage for Each Chromosome")
