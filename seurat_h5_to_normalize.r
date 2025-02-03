# Function to create a Seurat object from 10X Genomics data.
# This function takes the directory containing the 10X Genomics data as input,
# reads the data, creates a Seurat object, calculates mitochondrial gene percentage,
# and returns the Seurat object.
#
# Args:
#   dir: Path to the directory containing the 10X Genomics output files (e.g., barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz or filtered_feature_bc_matrix.h5).
#
# Returns:
#   A Seurat object.
#
# Example Usage:
# Assuming your 10X Genomics data is in a directory named "raw_feature_bc_matrix".
# You can create a Seurat object like this:
#
# seurat_object <- createSeuratObjectFrom10X("raw_feature_bc_matrix")
#
# Now 'seurat_object' is a Seurat object ready for further analysis.

createSeuratObjectFrom10X <- function (dir)
{
    # Load the Seurat library.
    # Seurat is a package in R designed for single-cell RNA-seq data analysis.
    library(Seurat)

    # Read 10X Genomics formatted data from the specified directory.
    # Read10X_h5() function reads 10X Genomics Cell Ranger output in HDF5 format (typically 'filtered_feature_bc_matrix.h5').
    # It returns a sparse matrix.
    object <- Read10X_h5(dir)

    # Check if the object is a list and contains "Gene Expression" as a named element.
    # Recent versions of Cell Ranger output h5 files that can contain multiple datatypes (e.g., Gene Expression, Antibody Capture, CRISPR Guide Capture).
    # This condition checks if the loaded object is a list and if it has an element named "Gene Expression".
    if (!is.null(names(object)) & ("Gene Expression" %in% names(object))) {
        # If "Gene Expression" is present, extract that specific data.
        # This ensures that we are working with the Gene Expression data if multiple datatypes are present in the 10X data.
        object <- object[["Gene Expression"]]
    }

    # Create a Seurat object from the loaded data.
    # CreateSeuratObject() function takes a count matrix as input and creates a Seurat object.
    # 'project' argument is used to name the project, here we use the directory name for simplicity.
    object <- CreateSeuratObject(object, project = dir)

    # Calculate the percentage of reads that map to mitochondrial genes.
    # PercentageFeatureSet() calculates the percentage of counts originating from a set of features.
    # Here, we use 'pattern = "^MT-"' to identify mitochondrial genes (genes starting with "MT-").
    # 'col.name = "percent.mt"' specifies the name of the metadata column to store these percentages.
    object <- PercentageFeatureSet(object, pattern = "^MT-",
        col.name = "percent.mt")

    # Return the created Seurat object.
    # This Seurat object now contains gene expression data and metadata including the percentage of mitochondrial genes.
    return(object)
}
