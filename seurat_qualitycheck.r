# Function to perform quality check on a Seurat object and generate quality control plots.
# This function takes a Seurat object as input, calculates the percentage of mitochondrial genes
# based on a user-defined pattern, and generates violin plots and feature scatter plots for quality assessment.
#
# Args:
#   seurat_object: A Seurat object.
#   mt_pattern:  Pattern to identify mitochondrial genes. Default is "^MT-" (for human).
#                You can change this to "^mt-" for mouse or other species-specific patterns.
#
# Returns:
#   A list containing three ggplot2 plots:
#     - Violin plot of nFeature_RNA, nCount_RNA, and percent.mt.
#     - Feature scatter plot of nCount_RNA vs percent.mt.
#     - Feature scatter plot of nCount_RNA vs nFeature_RNA.
#
# Example Usage for Human data:
# Assuming you have a Seurat object named 'seurat_obj_human'.
# quality_plots_human <- seurat_qualitycheck(seurat_obj_human)
# quality_plots_human[[1]] # To display the violin plot
# quality_plots_human[[2]] # To display the scatter plot of nCount_RNA vs percent.mt
# quality_plots_human[[3]] # To display the scatter plot of nCount_RNA vs nFeature_RNA
#
# Example Usage for Mouse data:
# Assuming you have a Seurat object named 'seurat_obj_mouse' and mitochondrial genes are prefixed with "mt-".
# quality_plots_mouse <- seurat_qualitycheck(seurat_obj_mouse, mt_pattern = "^mt-")
# quality_plots_mouse[[1]] # To display the violin plot for mouse data
# quality_plots_mouse[[2]] # To display the scatter plot of nCount_RNA vs percent.mt for mouse data
# quality_plots_mouse[[3]] # To display the scatter plot of nCount_RNA vs nFeature_RNA for mouse data


seurat_qualitycheck <- function(seurat_object, mt_pattern = "^MT-"){
    # Initialize an empty list to store the plots.
    plot_list <- list()

    # Set the default assay to "RNA".
    # This ensures that the operations are performed on the RNA assay of the Seurat object.
    DefaultAssay(seurat_object) <- "RNA"

    # Calculate the percentage of reads that map to mitochondrial genes.
    # PercentageFeatureSet() calculates the percentage of counts originating from a set of features.
    # Here, we use 'pattern = mt_pattern' to identify mitochondrial genes based on the provided pattern.
    # 'col.name = "percent.mt"' specifies the name of the metadata column to store these percentages.
    # By default, mt_pattern is set to "^MT-" (for human mitochondrial genes).
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = mt_pattern)

    # Generate a violin plot to visualize quality metrics.
    # VlnPlot() creates violin plots for specified features.
    # Here, we plot 'nFeature_RNA' (number of features), 'nCount_RNA' (number of counts), and 'percent.mt' (mitochondrial gene percentage).
    # 'ncol = 3' arranges the plots in 3 columns.
    VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) -> p1

    # Generate a feature scatter plot of nCount_RNA vs percent.mt.
    # FeatureScatter() creates scatter plots to visualize the relationship between two features.
    # This plot helps to see if there is a correlation between total counts and mitochondrial gene percentage.
    p2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")

    # Generate a feature scatter plot of nCount_RNA vs nFeature_RNA.
    # This plot helps to see the relationship between the number of counts and the number of detected genes per cell.
    # It's useful for identifying potential low-quality cells or doublets.
    p3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    # Store the generated plots in the list.
    plot_list <- list(p1, p2, p3)

    # Return the list of plots.
    return (plot_list)
}
