# Load required libraries

library(phangorn)      # For phylogenetic analysis
library(readxl)        # For reading Excel files
library(ape)           # For tree manipulations and plotting
library(vegan)         # For computing distances
library(itol.toolkit)  # For working with iTOL
library(cluster)       # For clustering methods
library(data.table)    # For data manipulation

# read gene presence-absence gene data from excel file
panaroo_data <- as.data.frame(read_excel("all_panaroo.xlsx"))

# Extract the gene names and transpose the data
genes <- panaroo_data$Gene
panaroo_data <- t(as.data.frame(panaroo_data[, -1]))

# Assign the gene names as column names
colnames(panaroo_data) <- genes

# Store the gene names as a separate column
genes <- row.names(panaroo_data)
panaroo_data <- data.frame(genes = genes, panaroo_data)
row.names(panaroo_data) <- NULL

# Clean up the gene names to remove the numbering
# panaroo_data$genes <- gsub("\\.\\d+\\.1", "", panaroo_data$genes)
# panaroo_data$genes <- gsub("\\.\\d+\\.1.fasta", "", panaroo_data$genes)
# panaroo_data$genes <- gsub("\\.\\d+\\.fasta", "", (panaroo_data$genes))
# panaroo_data$genes <- gsub("^NZ_", "", panaroo_data$genes)

# Replace missing values with 0
panaroo_data[is.na(panaroo_data)] <- 0

# Set the gene names as row names
row.names(panaroo_data) <- panaroo_data$genes

# Compute binary distance matrix
d = dist(panaroo_data, method = "binary")

# Perform hierarchical clustering
hc_2 = hclust(d, method = "ward.D2")

# Plot the hierarchical clustering dendrogram
plot(hc_2)

# Convert hierarchical clustering result to phylogenetic tree format
newwick_tree <- as.phylo(hc_2)

# Save the phylogenetic tree in Newick format
write.tree(newwick_tree, file = "final_tree.newick")




