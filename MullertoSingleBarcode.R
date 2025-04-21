#Muller to Single Barcode
#Goal: Maintain color code across plot types

experiment <- "04.02.2025_06.15_AM"
fish <- 6

# Import packages
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(ggpubr)
library(rstudioapi)

#1. IMPORT DATA
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load metadata
meta_directory <- paste0("./experiments/", experiment, '/passed_samples_info.csv')
meta <- read.table(meta_directory, header=TRUE, sep=",")
meta <- meta %>% arrange(sample_id)

# Load barcode data
barcode_directory <- paste0("./experiments/", experiment, "/matrix/fish", fish, "_barcodeMatrix_3filters.tsv")
barcode.mat <- fread(barcode_directory, header=TRUE, sep="\t")
barcode.mat[is.na(barcode.mat)] <- 0

# Convert to long format
barcodes_long <- melt(barcode.mat, id.vars = "barcode", variable.name = "sample_id", value.name = "count")

# Convert barcode column to factor
barcodes_long$barcode <- as.factor(barcodes_long$barcode)

# Merge metadata
barcodes_long <- merge(barcodes_long, meta, by = "sample_id")

# Compute total reads per sample
total_reads <- barcodes_long %>% group_by(sample_id) %>% summarize(total = sum(count))
barcodes_long <- merge(barcodes_long, total_reads, by = "sample_id")
barcodes_long <- barcodes_long %>% mutate(proportion = count / total)

# Determine the first timepoint for each fish
first_timepoint <- barcodes_long %>% group_by(fish) %>% summarize(first_age = min(age))
barcodes_long <- merge(barcodes_long, first_timepoint, by = "fish")

# Get barcode proportions at the first timepoint
barcode_order <- barcodes_long %>%
  filter(age == first_age) %>%
  group_by(barcode) %>%
  summarize(first_proportion = mean(proportion), .groups = 'drop') %>%
  arrange(desc(first_proportion)) %>%
  pull(barcode)

# Convert barcode to ordered factor
barcodes_long$barcode <- factor(barcodes_long$barcode, levels = barcode_order)

# Generate color palette for barcodes
library(Polychrome)
set.seed(25)
nb.cols <- length(unique(barcodes_long$barcode))
palette <- createPalette(nb.cols - 1, c("#FF4D6FFF", "#579EA4FF", "#DF7713FF", "#F9C000FF", "#86AD34FF", "#5D7298FF", "#81B28DFF", "#7E1A2FFF", "#2D2651FF", "#C8350DFF", "#BD777AFD"), M = 1000)
color_map <- c("NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE" = "grey", setNames(palette, setdiff(unique(barcodes_long$barcode), "NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE_NONE")))
colors <- color_map[barcodes_long$barcode]

#MULLER PLOTS
#subset by replicates
rep1 <- barcodes_long[barcodes_long$replicate %in% "1", ]
rep1$age <- as.factor(rep1$age)
p4 <- ggplot(rep1, aes(x=age, y=proportion, fill=barcode, group=barcode)) + 
  geom_area() +
  scale_fill_manual(values=colors) +
  scale_x_discrete(breaks = unique(rep1$age), expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x= "age") + 
  #ggtitle(paste0("Fish ", fish, ", replicate 1")) + 
  ggtitle(paste0("Fish ", fish)) + 
  theme(text=element_text(size=15)) + 
  theme(legend.position="none")
p4

rep2 <- barcodes_long[barcodes_long$replicate %in% "2", ]
rep2$age <- as.factor(rep2$age)
p5 <- ggplot(rep2, aes(x=age, y=proportion, fill=barcode, group=barcode)) + 
  geom_area() +
  scale_fill_manual(values=colors) +
  scale_x_discrete(breaks = unique(barcodes_long$age), expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x= "age") + 
  #ggtitle(paste0("Fish ", fish, ", replicate 2")) + 
  ggtitle(paste0("Fish ", fish)) + 
  theme(text=element_text(size=15)) + 
  theme(legend.position="none")
p5

rep3 <- barcodes_long[barcodes_long$replicate %in% "3", ]
rep3$age <- as.factor(rep3$age)
p6 <- ggplot(rep3, aes(x=age, y=proportion, fill=barcode, group=barcode)) + 
  geom_area() +
  scale_fill_manual(values=colors) +
  scale_x_discrete(breaks = unique(barcodes_long$age), expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x= "age") + 
  #ggtitle(paste0("Fish ", fish, ", replicate 3")) + 
  ggtitle(paste0("Fish ", fish)) + 
  theme(text=element_text(size=15)) + 
  theme(legend.position="none")
p6

ggarrange(p4,p5,p6, nrow=1, legend=NULL)
#with legend
ggarrange(p4,p5,p6, ncol=1, common.legend=TRUE, legend ="right", legend.text=element_text(size=5))

#SINGLE BARCODE PLOTS
# Compute mean proportion per barcode at each age
barcode_summary <- barcodes_long %>%
  group_by(fish, barcode, age) %>%
  summarize(mean_proportion = mean(proportion), .groups = 'drop')

# Plot faceted by barcode with shared y-axis and all x-axis ticks labeled
p1 <- ggplot(barcodes_long, aes(x = age, y = proportion)) +
  geom_point(color = "black", alpha = 0.5) +  # Individual replicate data points in black
  geom_line(data = barcode_summary, aes(x = age, y = mean_proportion, color = barcode), linewidth = 1) +
  scale_color_manual(values = colors) +  # Ensure only the "NONE" barcode is grey
  facet_wrap(~barcode, scales = "free_y") +  # Shared y-axis = fixed
  scale_x_continuous(breaks = unique(barcodes_long$age), expand = c(0,0)) +  # Ensure all x values are shown
  ggh4x::facet_wrap2(~ barcode, scales = "free_y", strip.position = "top", axes = "all") +  # Show x-axis labels in all facets
  theme_minimal() +
  labs(x = "Age", y = "Proportion") +
  ggtitle(paste0("Fish ", fish)) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 7),  # Keep facet labels readable
    axis.text.x = element_text(angle = 45, vjust = 1, size = 6),  # Rotate x-axis labels for readability
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),  # Larger, bold x-axis label
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),  # Larger, bold y-axis label
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Larger, bold, centered title
    panel.spacing = unit(1, "lines"),  # Ensure enough space for labels
    strip.placement = "outside",  # Keep facet labels visible
    panel.grid.minor.x = element_blank(),  # Remove excess gridlines
    )

print(p1)