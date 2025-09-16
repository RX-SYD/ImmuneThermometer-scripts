############################################
# t-SNE analysis of tumor immune signatures
# Author: RX2025
# Description:
#   Perform t-SNE on immune- and stromal-related signature scores 
#   across tumor samples, and visualize by cancer type.
############################################

# Set working directory (modify as needed)
setwd("YOUR FOLDER")

# Install tsne package if not already installed
if (!require("tsne")) {
  install.packages("tsne")
  library(tsne)
} else {
  library(tsne)
}

# Load signature matrix (rows: samples, columns: signatures)
# Ensure 'cancer_type' column exists in the input file
tumor_immune_signatures <- read.csv("YOUR MATRIX.csv", row.names = 1)

# Preview first few entries of selected columns (example: columns 21–22)
head(tumor_immune_signatures[, 21:22])

# Define colors for each cancer type
colors <- rainbow(length(unique(tumor_immune_signatures$cancer_type)))
names(colors) <- unique(tumor_immune_signatures$cancer_type)

# Run t-SNE (2D embedding with perplexity = 50)
tsne_result <- tsne(tumor_immune_signatures[, 1:21], k = 2, perplexity = 50)

# Plot t-SNE results
plot(tsne_result,
     col = colors[tumor_immune_signatures$cancer_type],
     pch = 16,
     axes = FALSE,
     cex = 0.7)

# Add legend
legend("topright",
       title = "Cancer Type",
       inset = 0.01,
       legend = unique(tumor_immune_signatures$cancer_type),
       pch = 16,
       cex = 0.6,
       col = unique(colors[tumor_immune_signatures$cancer_type]))

# Save t-SNE results
saveRDS(tsne_result, "tsne_tumor_immune_signatures.rds")
