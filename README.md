# spatial_transcriptomics_analysis.R
#This repository provides a comprehensive R script for analyzing single-cell and spatial transcriptomics data. It includes steps for quality control, normalization, dimensionality reduction, clustering, and visualization. The script also integrates cell type annotation using SingleR and celldex and identifies spatially variable features, 
# Spatial Transcriptomics Analysis

This repository contains an R script for comprehensive analysis of single-cell and spatial transcriptomics data using Seurat and related packages.

## Features
- Quality control, normalization, and scaling of spatial data.
- Dimensionality reduction, clustering, and visualization.
- Cell type annotation using SingleR and celldex.
- Identification of spatially variable features.

## Prerequisites
- R (version 4.0 or higher)
- Required R packages: `Seurat`, `SingleR`, `celldex`, `ggplot2`, `patchwork`, `dplyr`

## Installation
Install the necessary R packages:
```R
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
BiocManager::install(c("SingleR", "celldex"))

library(Seurat)
library(dplyr) # Data manipulation
library(ggplot2)
library(patchwork)
library(limma) # Optional

# Set file paths
data_dir <- "C:\\Users\\rushi\\Downloads\\V1_Mouse_Brain_Sagittal_Posterior_Section_2_spatial (3)"

h5_mat_name <- file.path(data_dir, "filtered_feature_bc_matrix.h5")
spatial_folder_path <- file.path(data_dir, "spatial")

# Load a 10x Genomics Visium Spatial Experiment into a Seurat object
brain_data <- Seurat::Load10X_Spatial(
  data.dir = data_dir, 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  slice = "slice1", 
  filter.matrix = TRUE, 
  to.upper = FALSE
)

# Output the Seurat object dimensions and other details
print(dim(brain_data)) # Dimensions of the Seurat object
print(nrow(brain_data)) # Number of features (genes)
print(ncol(brain_data)) # Number of samples (spots)
head(rownames(brain_data), n = 5)
tail(colnames(brain_data), n = 5)

# Examine metadata and counts
print(class(brain_data[[]]))
print(colnames(brain_data[[]]))
head(brain_data@meta.data)
brain_data$nCount_Spatial[1:3]
brain_data[['Spatial']] 
brain_data[['slice1']] 
brain_data@assays$Spatial@counts[5:10, 1:3]

brain_data[['Spatial']]@meta.features
head(brain_data[['Spatial']][[]])

# Calculate mitochondrial percentage
brain_data[["percent.mt"]] <- PercentageFeatureSet(brain_data, pattern = "^mt-")

# Plotting
VlnPlot(brain_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
        pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot1 <- FeatureScatter(brain_data, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(brain_data, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + NoLegend()
plot1 + plot2

# Subsetting the data based on QC metrics
brain_subset <- subset(
  brain_data, 
  subset = nFeature_Spatial < 8000 & nFeature_Spatial > 1000 & 
    nCount_Spatial < 50000 & percent.mt < 30
)

print(paste("Filter out", ncol(brain_data) - ncol(brain_subset), 
            "samples because of the outlier QC metrics, with", ncol(brain_subset),
            "samples left."))

SpatialFeaturePlot(brain_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "bottom")  

brain_norm <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)
print(names(brain_norm)) dim(brain_norm@assays$SCT@counts)
dim(brain_norm@assays$SCT@scale.data) 
brain_obj <- RunPCA(brain_norm, assay = "SCT", verbose = FALSE)
# compute K nearest neighbors (KNN)
brain_obj <- FindNeighbors(brain_obj, reduction = "pca", dims = 1:30)
# Leiden algorithm for community detection
brain_obj <- FindClusters(brain_obj, verbose = FALSE)
# PCA result is the default UMAP input, use dimensions 1:30 as input features
brain_obj <- RunUMAP(brain_obj, reduction = "pca", dims = 1:30)

plot3 <- DimPlot(brain_obj, reduction = "umap", label = TRUE) + NoLegend()
plot4 <- SpatialDimPlot(brain_obj, label = TRUE, label.size = 3) + NoLegend()
plot3 + plot4

brain_obj@reductions
# identity class of each sample
table(brain_obj@active.ident)

# find all markers of cluster 10
cluster10_markers <- FindMarkers(brain_obj, ident.1 = 10, min.pct = 0.25)
head(cluster10_markers, n = 5)

SpatialFeaturePlot(object = brain_obj, 
                   features = rownames(cluster10_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
# this code chunk is not evaluated for now because of time constraints
brain_obj_markers <- FindAllMarkers(brain_obj, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
brain_obj_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top10 <- brain_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(brain_obj, features = top10$gene) + NoLegend()


brain_moransi <- FindSpatiallyVariableFeatures(
  brain_obj, assay = "SCT", 
  features = VariableFeatures(brain_obj)[1:10],
  selection.method = "moransi") 
moransi_output_df <- brain_moransi@assays$SCT@meta.features %>%
  na.exclude
head(moransi_output_df[order(moransi_output_df$MoransI_observed, decreasing = T), ])

top_features_moransi <- head(
  SpatiallyVariableFeatures(brain_moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(brain_moransi, 
                   features = top_features_moransi, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "Top 3 genes with the largest Moran's I",
    subtitle = "among 10 top variable genes for illustration purposes")
bottom_features_moransi <- tail(
  SpatiallyVariableFeatures(brain_moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(brain_moransi, 
                   features = bottom_features_moransi, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "Bottom 3 genes with the smallest Moran's I",
    subtitle = "among 10 top variable genes for illustration purposes")

brain_variogram <- FindSpatiallyVariableFeatures(
  brain_obj, assay = "SCT", 
  features = VariableFeatures(brain_obj)[1:10],
  selection.method = "markvariogram")  
variogram_output_df <- brain_variogram@assays$SCT@meta.features %>%
  na.exclude # there are NA rows b/c we only calculated the variogram for 10 genes
head(variogram_output_df[order(variogram_output_df$r.metric.5), ])

top_features_variogram <- head(
  SpatiallyVariableFeatures(brain_variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(brain_variogram, 
                   features = top_features_variogram, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "3 genes with the top spatially variable rank (by mark-variogram)",
    subtitle = "among 10 top variable genes for illustration purposes")
bottom_features_variogram <- tail(
  SpatiallyVariableFeatures(brain_variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(brain_variogram, 
                   features = bottom_features_variogram, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "3 genes with the bottom spatially variale rank (by mark-variogram)",
    subtitle = "among 10 top variable genes for illustration purposes")
# Load required libraries
library(SingleR)
library(celldex)  # Provides reference datasets for SingleR

# Load reference data for mouse (you can replace this with a specific reference dataset if needed)
ref_data <- celldex::MouseRNAseqData()

# Extract the counts and metadata from the Seurat object
seurat_counts <- brain_obj@assays$SCT@counts
seurat_metadata <- brain_obj@meta.data

# Run SingleR for cell type annotation
singleR_results <- SingleR(
  test = as.matrix(seurat_counts),  # Counts matrix from Seurat object
  ref = ref_data,                  # Reference dataset
  labels = ref_data$label.main     # Use main cell type labels
)

# Add SingleR predicted labels to Seurat metadata
brain_obj$SingleR_label <- singleR_results$labels

# Annotate the UMAP with SingleR labels
umap_annotated_plot <- DimPlot(
  brain_obj,
  reduction = "umap",
  group.by = "SingleR_label",  # Use SingleR labels for grouping
  label = TRUE,
  repel = TRUE
) + ggtitle("UMAP with SingleR Cell Type Annotation")

# Display the annotated UMAP
print(umap_annotated_plot)
