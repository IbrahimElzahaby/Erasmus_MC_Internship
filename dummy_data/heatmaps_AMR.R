
# load the required package
library(readxl)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(reshape2)
# Read the Excel data into a matrix
excel_data <- read_excel("AMR_heatmap.xlsx")

# Extract the row and column names
row_names <- excel_data[[1]]
col_names <- colnames(excel_data)[-1]  # Exclude the first column which is for row names

# Convert the data to a numeric matrix
data_matrix <- as.matrix(excel_data[-1])  # Exclude the first column

# Create a heatmap with limited legend levels
AMR <- pheatmap(data_matrix,
                cluster_rows = F, cluster_cols = TRUE, 
                col = colorRampPalette(c("#A03C78", "#FDFF00", "#38E54D"))(100),
                labels_row = row_names, labels_col = col_names,
                fontsize_row = 5, fontsize_col = 5,
                cellwidth = 5, cellheight = 5,
                border_color = "white",
                legend_breaks = c(0, 1, 2))  # Set the value scale to 0, 1, 2

# Print the heatmap
print(AMR)


ggsave("AMR_heatmap.png", AMR, width = 10, height = 6, dpi = 400)


