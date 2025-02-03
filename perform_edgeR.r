#!/usr/bin/env Rscript

# Function to perform generalized edgeR analysis for differential gene expression.
# This function takes count data, group information, and a list of comparisons as input,
# and performs edgeR analysis to identify differentially expressed genes for each comparison.
#
# Usage:
# 1. Prepare your count data as a CSV file (e.g., counts.csv) where rows are genes and columns are samples.
#    The first column should be gene names and the first row should be sample names.
#    Example:
#       gene,sample1,sample2,sample3,sample4,sample5,sample6
#       gene1,100,120,50,60,200,220
#       gene2,50,60,20,30,100,110
#       ...
#
# 2. Prepare your group information as a CSV file (e.g., group.csv) with one column indicating the group for each sample.
#    The order of groups in this file must correspond to the column order in your count data.
#    Example:
#       control
#       control
#       treat1
#       treat1
#       treat2
#       treat2
#
# 3. Prepare a comparisons list as a CSV file (e.g., comparisons.csv) with two columns: 'target' and 'reference' groups for each comparison.
#    Each row defines a comparison you want to perform (target vs reference).
#    Example:
#       treat1,control
#       treat2,control
#       treat2,treat1
#
# 4. Run this script from the command line using Docker (or Rscript directly if you have dependencies installed):
#    docker run -v $(pwd):/data your-dockerhub-username/edger-generalized \
#        Rscript script.R \
#        -c /data/counts.csv \
#        -g /data/group.csv \
#        -v /data/comparisons.csv \
#        -o /data/output
#
#    Replace 'your-dockerhub-username/edger-generalized' with your Docker Hub image name.
#    Make sure 'counts.csv', 'group.csv', and 'comparisons.csv' are in the same directory where you run the command.
#    Results will be saved in the 'output' directory within the same directory.
#
# Arguments:
#   countdata: Count data frame (rows: gene names, columns: sample names).
#   group:     Factor vector of sample groups (corresponding to the column order of countdata).
#   comparisons: A list of comparisons. Each element is a vector of two group names to compare (e.g., list(c("groupA", "groupB"), c("groupA", "groupC"))).

library(edgeR)
library(optparse) # For command line argument parsing

# Define command line arguments
option_list <- list(
  make_option(c("-c", "--countdata"), type="character", default=NULL,
              help="Count data file (CSV format)", metavar="character"),
  make_option(c("-g", "--groupfile"), type="character", default=NULL,
              help="Group information file (CSV format, one column factor)", metavar="character"),
  make_option(c("-v", "--comparisonsfile"), type="character", default=NULL,
              help="Comparisons list file (CSV format, two columns: target, reference)", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default="output",
              help="Output directory", metavar="character")
)

opt_parser <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt_parser$countdata) || is.null(opt_parser$groupfile) || is.null(opt_parser$comparisonsfile)){
  print_help(opt_parser)
  stop("Count data, group file, and comparisons file must be provided.", call.=FALSE)
}

# Get filenames from arguments
countdata_file <- opt_parser$countdata
groupfile <- opt_parser$groupfile
comparisonsfile <- opt_parser$comparisonsfile
output_dir <- opt_parser$outputdir

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load data and configuration files
countdata <- read.csv(countdata_file, row.names=1)
group <- factor(read.csv(groupfile, header=FALSE)[,1])
comparisons_df <- read.csv(comparisonsfile, header=FALSE)
comparisons_list <- list()
for (i in 1:nrow(comparisons_df)) {
  comparisons_list[[i]] <- c(as.character(comparisons_df[i,1]), as.character(comparisons_df[i,2]))
}


perform_edgeR_generalized <- function(countdata, group, comparisons) {
  # 1. Create DGEList object
  y <- DGEList(counts=countdata, group=group)

  # 2. Filtering (remove low expression genes)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]

  # 3. Normalization (TMM normalization)
  y <- calcNormFactors(y)

  # 4. Dispersion estimation
  design <- model.matrix(~0+group) # Create model matrix for group comparison
  colnames(design) <- levels(group) # Set column names
  y <- estimateDisp(y, design)

  # 5. BCV plot (check dispersion estimation)
  plotBCV(y)

  # 6. MDS plot (check sample relationships)
  plotMDS(y, col=as.numeric(group)) # Color-code by condition
  legend("topright", legend=levels(group), col=1:nlevels(group), pch=16) # Add legend

  # List to store results
  qlf_results <- list()

  # 7. Differential gene expression analysis (GLM QLF test) - Run for each comparison
  for (comparison in comparisons) {
    target_group <- comparison[1]
    reference_group <- comparison[2]
    comparison_name <- paste0(target_group, "_vs_", reference_group)

    con_target_vs_ref <- makeContrasts(contrasts = paste0(target_group, " - ", reference_group), levels=design)
    fit_target_vs_ref <- glmQLFit(y, design)
    qlf_target_vs_ref <- glmQLFTest(fit_target_vs_ref, contrast=con_target_vs_ref)

    # 8. Summary of results and MD plot
    print(paste("Summary for comparison:", comparison_name)) # Display comparison name
    print(summary(decideTests(qlf_target_vs_ref)))
    plotMD(qlf_target_vs_ref, main = comparison_name) # Add comparison name to MD plot title

    # 9. Save results (RDS file)
    date_prefix <- format(Sys.Date(), "%y%m%d")
    filename <- paste0(output_dir, '/', date_prefix, '_qlf_', comparison_name, '.rds') # Include comparison name in filename
    saveRDS(qlf_target_vs_ref, filename)
    qlf_results[[comparison_name]] <- qlf_target_vs_ref # Store results in list (comparison name as key)
  }

  return(qlf_results) # Return list of results for all comparisons
}


# Run the function
results <- perform_edgeR_generalized(countdata, group, comparisons_list)

# Completion message
cat("edgeR analysis completed. Results saved in:", output_dir, "\n")
