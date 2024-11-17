#Author: Chrats Melkonian
# =============================================
# Part A
# =============================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(factoextra)
library(ggrepel)
library(uwot)
library(patchwork)

# Set the working directory
setwd(file.path(Sys.getenv("HOME"), "Desktop/Git/Genome_mining_PPB_Halo"))

# Load data
HALO_KO_matrix <- read.csv("Data/HALO_KO_matrix.csv")
data_all_KOs <- read.csv("Data/data_functions_KOs_markers.csv")
metadata_final_selection <- read.csv("Data/metadata_phenotype_annotation_final.csv")

# Update genome names in the HALO_KO_matrix
HALO_KO_matrix <- HALO_KO_matrix %>%
  mutate(X = recode(X,
                    "Rhodocista_1709002.5.PATRIC" = "NCBI_6",
                    "Rhabdochromatium_48729.6.PATRIC" = "NCBI_5",
                    "Ectothiorhodospira.variabilis" = "NCBI_4",
                    "Ectothiorhodospira.shaposhnikovii" = "NCBI_3",
                    "Ectothiorhodospira.lacustris" = "NCBI_2",
                    "Ectothiorhodospira.haloalkaliphila" = "NCBI_1"))

# Extract genus from genome_name and remove duplicate entries
metadata_final_selection$genus <- unlist(lapply(strsplit(metadata_final_selection$genome_name, " "), function(x) x[1]))
metadata_final_selection <- distinct(metadata_final_selection, genome_id, .keep_all = TRUE)

# Filter HALO_KO_matrix columns for relevant orthologies
filtered_matrix <- HALO_KO_matrix[, colnames(HALO_KO_matrix) %in% paste("ko.", data_all_KOs$Orthology, sep = "")]

# Function to remove columns with all zeros or ones
remove_all_zeros_ones_columns <- function(matrix) {
  non_constant_columns1 <- apply(matrix, 2, function(col) all(col == 1))
  return(matrix[, !non_constant_columns1, drop = FALSE])
}

# Update row names and filter rows and columns based on metadata
rownames(HALO_KO_matrix) <- HALO_KO_matrix$X
HALO_KO_matrix <- HALO_KO_matrix[, -1]
rows_to_keep <- metadata_final_selection$genome_id %in% rownames(HALO_KO_matrix)
metadata_final_selection_filtered <- metadata_final_selection[rows_to_keep, ]
metadata_final_selection_unique <- distinct(metadata_final_selection_filtered, genome_name, .keep_all = TRUE)
rows_to_keep1 <- rownames(HALO_KO_matrix) %in% metadata_final_selection_unique$genome_id
HALO_KO_matrix_filtered <- HALO_KO_matrix[rows_to_keep1, , drop = FALSE]
filtered_matrix <- remove_all_zeros_ones_columns(HALO_KO_matrix_filtered)

# Define marker groups and assign classifications
groups_markers <- c("CO2 fixation", "CO fixation anaerobic", "CO2/CO fixation", "CO fixation aerobic")
list_markers <- list()
j <- 1
for (i in groups_markers) {
  temp <- rep("Other", nrow(metadata_final_selection))
  KOs <- data_all_KOs$Orthology[grep(i, data_all_KOs$Marker_info)]
  filtered_matrix_temp <- filtered_matrix[, colnames(filtered_matrix) %in% paste("ko.", KOs, sep = "")]
  ids <- names(which(rowSums(filtered_matrix_temp) > 0)) # at least 1 marker gene present 
  temp[metadata_final_selection$genome_id %in% ids] <- i
  list_markers[[j]] <- temp
  j <- j + 1
}
names(list_markers) <- groups_markers

# Add markers to metadata
metadata_final_selection_f <- metadata_final_selection
metadata_final_selection_f$Markers_CO2_f <- list_markers$`CO2 fixation`
metadata_final_selection_f$Markers_CO_f_anaerobic <- list_markers$`CO fixation anaerobic`
metadata_final_selection_f$Markers_CO2_CO_f <- list_markers$`CO2/CO fixation`
metadata_final_selection_f$Markers_CO_f_aerobic <- list_markers$`CO fixation aerobic`

# Perform UMAP for dimensionality reduction
umap_result <- umap(filtered_matrix, n_neighbors = 15, spread = 20, n_components = 2, metric = "hamming")
umap_df <- as.data.frame(umap_result)
umap_df$id <- rownames(umap_df)

# Merge UMAP results with metadata
umap_df_all <- merge(umap_df, metadata_final_selection_f, by.x = "id", by.y = "genome_id", all.x = TRUE)
umap_df_all <- umap_df_all[!grepl("strain not applicable", umap_df_all$genome_name), ]

# Define colors and plot settings
colors <- c('#1f78b4', '#a6cee3', "#7f0000", "#d7301f", "#33a02c")
alpha_values <- c(0.8, 0.4)
size_values <- c(3, 1)
shape_values <- c(18, 1)
# Create UMAP visualizations for each marker
g1<-ggplot(umap_df_all, aes(x = V1, y = V2, color= Phenotype_final,shape=Markers_CO2_f , alpha = Markers_CO2_f, size = Markers_CO2_f)) +
  geom_point() +
  theme_minimal() +
  labs(title = "",x="",y="UMAP2")+
  scale_color_manual(values = colors)+
  scale_alpha_manual(values = alpha_values) + # Set custom alpha values
  scale_size_manual(values = size_values)+
  scale_shape_manual(values = shape_values)+theme(
    legend.title = element_text(size = 14), # Increase legend title size
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 12),    # Increase axis text size
    plot.title = element_text(size = 16)    # Increase plot title size
  ) 

g2<-ggplot(umap_df_all, aes(x = V1, y = V2, color= Phenotype_final,shape=Markers_CO_f_anaerobic , alpha = Markers_CO_f_anaerobic, size = Markers_CO_f_anaerobic)) +
  geom_point() +
  theme_minimal() +
  labs(title = "",x="",y="")+
  scale_color_manual(values = colors)+
  scale_alpha_manual(values = alpha_values) + # Set custom alpha values
  scale_size_manual(values = size_values)+
  scale_shape_manual(values = shape_values)+theme(
    legend.title = element_text(size = 14), # Increase legend title size
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 12),    # Increase axis text size
    plot.title = element_text(size = 16)    # Increase plot title size
  ) 
shape_values <- c(16, 1)
g3<-ggplot(umap_df_all, aes(x = V1, y = V2, color= Phenotype_final,shape=Markers_CO2_CO_f , alpha = Markers_CO2_CO_f, size = Markers_CO2_CO_f)) +
  geom_point() +
  theme_minimal() +
  labs(title = "",x="UMAP1",y="UMAP2")+
  scale_color_manual(values = colors)+
  scale_alpha_manual(values = alpha_values) + # Set custom alpha values
  scale_size_manual(values = size_values)+
  scale_shape_manual(values = shape_values)+theme(
    legend.title = element_text(size = 14), # Increase legend title size
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 12),    # Increase axis text size
    plot.title = element_text(size = 16)    # Increase plot title size
  ) 
shape_values <- c(15, 1)
g4<-ggplot(umap_df_all, aes(x = V1, y = V2, color= Phenotype_final,shape=Markers_CO_f_aerobic , alpha = Markers_CO_f_aerobic, size = Markers_CO_f_aerobic)) +
  geom_point() +
  theme_minimal() +
  labs(title = "",x="UMAP1",y="")+
  scale_color_manual(values = colors)+
  scale_alpha_manual(values = alpha_values) + # Set custom alpha values
  scale_size_manual(values = size_values)+
  scale_shape_manual(values = shape_values)+theme(
    legend.title = element_text(size = 14), # Increase legend title size
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 12),    # Increase axis text size
    plot.title = element_text(size = 16)    # Increase plot title size
  ) 

library(patchwork)

combined_plot <- (g1 + g2) / (g3 + g4) + 
  plot_annotation(title = "UMAP of microbial functions",
                  subtitle = "2x2 Panel Layout with Different Highlighted Markers")+
  plot_layout(guides = "collect") 

# Print the combined plot


# Set the output to a PDF file with specified dimensions
pdf("Figures/UMAP.pdf", width = 12, height = 10)

print(combined_plot)

# Close the PDF device to save the file
dev.off()
# =============================================
# Part B
# =============================================

# Load necessary libraries

library(tidyr)
library(ggupset)
library(pheatmap)


# =============================================
# Part B1: UMAP Plot Visualization
# =============================================

# Create a UMAP plot with custom aesthetics
f1 <- ggplot(umap_df_all, aes(x = V1, y = V2, color = Phenotype_final)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(
    title = "",
    x = "UMAP1",
    y = "UMAP2",
    color = "Phenotype"
  ) +
  scale_color_manual(values = colors) +
  scale_alpha_manual(values = alpha_values) +  # Set custom alpha values
  scale_size_manual(values = size_values) +
  scale_shape_manual(values = shape_values) +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16)
  )

# =============================================
# Part B2: Data Transformation and UpSet Plot
# =============================================

# Extract specific columns from metadata
mydata <- metadata_final_selection_f[, c(92:95)]

# Replace "Other" with 0 and convert other values to 1
mydata[mydata == "Other"] <- 0
mydata[mydata != 0] <- 1

# Convert to a logical data frame
my_df <- as.data.frame(apply(mydata, 2, function(x) as.logical(as.numeric(x))))
rownames(my_df) <- metadata_final_selection_f[, 1]

# Rename columns to meaningful function names
colnames(my_df) <- c(
  "CO2 fixation", 
  "CO fixation anaerobic", 
  "CO2/CO fixation", 
  "CO fixation aerobic"
)

# Transform data to a tidy format for visualization
tidy_pathway_member <- my_df %>%
  as_tibble(rownames = "Isolates") %>%
  gather(Function, Member, -Isolates) %>%
  filter(Member) %>%
  select(-Member)

# Summarize and create an UpSet plot
f2 <- tidy_pathway_member %>%
  group_by(Isolates) %>%
  summarize(Function = list(Function)) %>%
  ggplot(aes(x = Function)) +
  geom_bar() +
  scale_x_upset() +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16)
  )

# =============================================
# Part B3: Combine Plots
# =============================================

# Combine UMAP and UpSet plots with annotations
combined_plot1 <- f1 + f2 + 
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(face = 'bold'),
    text = element_text(size = 20)
  )

# Print the combined plot
print(combined_plot1)

# Set the output to a PDF file with specified dimensions
pdf("Figures/UMAP_upset.pdf", width = 12, height = 6)

print(combined_plot1)

# Close the PDF device to save the file
dev.off()
# =============================================
# Part B4: Data Preparation for Heatmaps
# =============================================

# Clean column names in the filtered matrix
colnames(filtered_matrix) <- gsub("ko.", "", colnames(filtered_matrix))

# Match KOs with column names in the filtered matrix
matching_indices <- match(data_all_KOs$Orthology, colnames(filtered_matrix))
matching_indices[is.na(matching_indices)] <- 0
filtered_matrix_sel <- filtered_matrix[, matching_indices]

# Create annotation data for pathways and functions
annotation_data <- data.frame(
  Pathway = data_all_KOs$Pathway[match(
    unlist(lapply(strsplit(colnames(filtered_matrix_sel), "[.]"), function(x) x[1])),
    data_all_KOs$Orthology
  )],
  Function = data_all_KOs$Function[match(
    unlist(lapply(strsplit(colnames(filtered_matrix_sel), "[.]"), function(x) x[1])),
    data_all_KOs$Orthology
  )]
)
rownames(annotation_data) <- colnames(filtered_matrix_sel)

# Create row annotations for phenotypes
annotation_data_row <- data.frame(
  Phenotype = metadata_final_selection_f$Phenotype_final[
    match(rownames(filtered_matrix_sel), metadata_final_selection_f$genome_id)
  ]
)
rownames(annotation_data_row) <- rownames(filtered_matrix_sel)

# =============================================
# Part B5: Pathway Completeness Analysis
# =============================================

# Initialize a completeness matrix
filtered_matrix_completness <- c()

# Calculate completeness percentage for each pathway
for (i in unique(data_all_KOs$Pathway)) {
  temp <- data_all_KOs[data_all_KOs$Pathway == i, ]
  length_kos <- length(temp$Function)
  
  matching_indices <- match(temp$Orthology, colnames(filtered_matrix))
  matching_indices[is.na(matching_indices)] <- 0
  filtered_matrix_temp <- filtered_matrix[, matching_indices]
  
  if (is.null(ncol(filtered_matrix_temp))) {
    filtered_matrix_temp_perc <- filtered_matrix_temp * 100
  } else {
    filtered_matrix_temp_perc <- apply(filtered_matrix_temp, 1, function(x) (sum(x) / length_kos) * 100)
  }
  filtered_matrix_completness <- cbind(filtered_matrix_completness, filtered_matrix_temp_perc)
}

# Rename columns in the completeness matrix
colnames(filtered_matrix_completness) <- unique(data_all_KOs$Pathway)

# Select a subset for visualization
final <- filtered_matrix_completness[, c(8:14, 16:17)]

ann_colors = list(
  # Function = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33'),
  Phenotype = c(Halophile = '#1f78b4',Halophile_extra = '#a6cee3', Halotolerant = "#7f0000", Halotolerant_extra = "#d7301f", Unknown = "#33a02c")
)

# Generate a heatmap for pathway completeness

# Set the output to a PDF file with specified dimensions
pdf("Figures/heatmap_all.pdf", width = 12, height = 10)

print(pheatmap::pheatmap(final, cluster_cols = TRUE, annotation_row = annotation_data_row, annotation_colors = ann_colors)
)

# Close the PDF device to save the file
dev.off()

# =============================================
# Part C1: Logical Conditions for Fixation Types
# =============================================

# Define logical conditions for different fixation types
CO2fixation <- apply(my_df, 1, function(x) x[1] == TRUE && x[2] == FALSE && x[3] == FALSE && x[4] == FALSE)
CO2fixationCOfixaerob <- apply(my_df, 1, function(x) x[1] == TRUE && x[2] == FALSE && x[3] == FALSE && x[4] == TRUE)
COfixaerob <- apply(my_df, 1, function(x) x[1] == FALSE && x[2] == FALSE && x[3] == FALSE && x[4] == TRUE)

# =============================================
# Part C2: Generate Heatmaps for Fixation Types
# =============================================

# 2.1: CO2 Fixation Heatmap
final_CO2fixation <- final[names(which(CO2fixation == TRUE)), ]
pheatmap::pheatmap(
  final_CO2fixation,
  cluster_cols = FALSE,
  annotation_names_row = FALSE,
  main = "CO2 fixation",
  show_rownames = FALSE,
  cutree_cols = 9,
  annotation_row = annotation_data_row,
  annotation_colors = ann_colors
)

# 2.2: CO2 Fixation & Aerobic CO Fixation Heatmap
final_CO2fixationCOfixaerob <- final[names(which(CO2fixationCOfixaerob == TRUE)), ]
pheatmap::pheatmap(
  final_CO2fixationCOfixaerob,
  cluster_cols = FALSE,
  annotation_names_row = FALSE,
  main = "CO2 fixation & CO fixation aerobic",
  show_rownames = FALSE,
  cutree_cols = 9,
  annotation_row = annotation_data_row,
  annotation_colors = ann_colors
)

# 2.3: Aerobic CO Fixation Heatmap
final_COfixaerob <- final[names(which(COfixaerob == TRUE)), ]
pheatmap::pheatmap(
  final_COfixaerob,
  cluster_cols = FALSE,
  annotation_names_row = FALSE,
  main = "CO fixation aerobic",
  show_rownames = FALSE,
  cutree_cols = 9,
  annotation_row = annotation_data_row,
  annotation_colors = ann_colors
)

# =============================================
# Part C3: Merge Metadata and Completeness Data
# =============================================

merged_df <- merge(
  filtered_matrix_completness, 
  metadata_final_selection_f, 
  by.x = "row.names", 
  by.y = "genome_id", 
  all.x = TRUE
)

# =============================================
# Part C4: Aerobic CO Fixation Analysis
# =============================================

# Extract and subset data for aerobic CO fixation
final_COfixaerob <- merged_df[match(names(which(COfixaerob == TRUE)), merged_df$Row.names), ]
final_COfixaerob_sub <- final_COfixaerob[, c(9:15, 17:18, 33)]

# Calculate average presence and counts
avg_df <- final_COfixaerob_sub %>%
  group_by(genus) %>%
  summarise(count = n(), across(everything(), mean, na.rm = TRUE))

avg_df_long <- avg_df %>%
  pivot_longer(cols = -c(genus, count), names_to = "variable", values_to = "value")

# Generate heatmap for aerobic CO fixation
heatmap_plot <- ggplot(avg_df_long, aes(x = variable, y = genus, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#e34a33") +
  theme_minimal() +
  labs(title = "", x = "Functions", y = "Genus", fill = "Avg. presence (%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

# Generate count bar plot
avg_df_long_count <- avg_df_long %>%
  group_by(genus) %>%
  summarise(count = unique(count))

count_plot <- ggplot(as.data.frame(avg_df_long_count), aes(x = genus, y = count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "", y = "Genome count") +
  theme(legend.position = "none") +
  coord_flip() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Combine heatmap and count plot
h1<-heatmap_plot + count_plot + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', title = 'CO fixation aerobic') &
  theme(plot.tag = element_text(face = 'bold'))

# Set the output to a PDF file with specified dimensions
pdf("Figures/CO_fixation_aerobic.pdf", width = 12, height = 10)

print(h1)


# Close the PDF device to save the file
dev.off()
# =============================================
# Part C5: CO2 Fixation Analysis
# =============================================

# Subset data for CO2 fixation based on logical condition
final_CO2fixation <- merged_df[match(names(which(CO2fixation == TRUE)), merged_df$Row.names), ]
final_CO2fixation_sub <- final_CO2fixation[, c(9:15, 17:18, 33)]

# Calculate average values and genome count by genus
avg_df_CO2 <- final_CO2fixation_sub %>%
  group_by(genus) %>%
  summarise(count = n(), across(everything(), mean, na.rm = TRUE))

# Reshape data for heatmap
avg_df_long_CO2 <- avg_df_CO2 %>%
  pivot_longer(cols = -c(genus, count), names_to = "variable", values_to = "value")

# Create heatmap plot
heatmap_plot_CO2 <- ggplot(avg_df_long_CO2, aes(x = variable, y = genus, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#2b8cbe") +
  theme_minimal() +
  labs(x = "Functions", y = "Genus", fill = "Avg. presence (%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

# Calculate genome count by genus
avg_df_long_count_CO2 <- avg_df_long_CO2 %>%
  group_by(genus) %>%
  summarise(count = unique(count))

# Create genome count bar plot
count_plot_CO2 <- ggplot(as.data.frame(avg_df_long_count_CO2), aes(x = genus, y = count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y = "Genome count") +
  theme(legend.position = "none") + 
  coord_flip() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Combine heatmap and count plot
h2 <- heatmap_plot_CO2 + count_plot_CO2 + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', title = 'CO2 fixation') &
  theme(plot.tag = element_text(face = 'bold'))

# Set the output to a PDF file with specified dimensions
pdf("Figures/CO2_fixation.pdf", width = 12, height = 10)

print(h2)


# Close the PDF device to save the file
dev.off()
# =============================================
# Part C6: CO2 Fixation & Aerobic CO Fixation Analysis
# =============================================


# Subset data for CO2 fixation and aerobic CO fixation
final_CO2fixationCOfixaerob <- merged_df[match(names(which(CO2fixationCOfixaerob == TRUE)), merged_df$Row.names), ]
final_CO2fixationCOfixaerob_sub <- final_CO2fixationCOfixaerob[, c(9:15, 17:18, 33)]

# Calculate average values and genome count by genus
avg_df_COfixaerob <- final_CO2fixationCOfixaerob_sub %>%
  group_by(genus) %>%
  summarise(count = n(), across(everything(), mean, na.rm = TRUE))

# Reshape data for heatmap
avg_df_long_COfixaerob <- avg_df_COfixaerob %>%
  pivot_longer(cols = -c(genus, count), names_to = "variable", values_to = "value")

# Create heatmap plot
heatmap_plot_COfixaerob <- ggplot(avg_df_long_COfixaerob, aes(x = variable, y = genus, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#31a354") +
  theme_minimal() +
  labs(x = "Functions", y = "Genus", fill = "Avg. presence (%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

# Calculate genome count by genus
avg_df_long_count_COfixaerob <- avg_df_long_COfixaerob %>%
  group_by(genus) %>%
  summarise(count = unique(count))

# Create genome count bar plot
count_plot_COfixaerob <- ggplot(as.data.frame(avg_df_long_count_COfixaerob), aes(x = genus, y = count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y = "Genome count") +
  theme(legend.position = "none") + 
  coord_flip() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Combine heatmap and count plot
h3 <- heatmap_plot_COfixaerob + count_plot_COfixaerob + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', title = 'CO2 fixation & CO fixation aerobic') &
  theme(plot.tag = element_text(face = 'bold'))

# Set the output to a PDF file with specified dimensions
pdf("Figures/CO2_fixation_&_CO_fixation_aerobic.pdf", width = 12, height = 10)

print(h3)

# Close the PDF device to save the file
dev.off()