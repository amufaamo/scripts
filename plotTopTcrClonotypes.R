# 必要なライブラリ
library(Seurat)
library(dplyr)
library(ggplot2)
library(forcats)

# ==============================================================================
# TCR上位クローンの頻度をプロットする関数
# ==============================================================================
#' @title 上位TCRクローンの可視化
#' @description Seuratオブジェクトから指定されたCDR3カラムの情報を集計し、
#'              出現頻度上位のクローンを棒グラフで表示します。
#'
#' @param seurat_object Seuratオブジェクト。
#' @param cdr3_column (character) プロットしたいメタデータ内のCDR3カラム名。
#'                    例: `"TRB_cdr3"` や `"TRA_cdr3"`
#' @param top_n (numeric) 上位何件のクローンを表示するか。デフォルトは15。
#' @param plot_by (character) "count"（細胞数）または "percentage"（割合）のどちらで
#'                プロットするかを指定。デフォルトは "count"。
#' @param plot_color (character) 棒グラフの色を指定。
#'
#' @return ggplotオブジェクト。
#'
plotTopTcrClonotypes <- function(seurat_object,
                                 cdr3_column = "TRB_cdr3",
                                 top_n = 15,
                                 plot_by = "count",
                                 plot_color = "#2a9d8f") {

  # --- 1. 入力チェック ---
  if (!inherits(seurat_object, "Seurat")) {
    stop("Error: 'seurat_object' must be a Seurat object.")
  }
  if (!cdr3_column %in% colnames(seurat_object@meta.data)) {
    stop("Error: '", cdr3_column, "' not found in the metadata of the Seurat object.")
  }
  if (!plot_by %in% c("count", "percentage")) {
    stop("Error: 'plot_by' must be either 'count' or 'percentage'.")
  }

  # --- 2. データの準備と集計 ---
  cdr3_data <- seurat_object[[cdr3_column, drop = TRUE]]
  cdr3_filtered <- na.omit(cdr3_data)

  if (length(cdr3_filtered) == 0) {
    warning("No valid CDR3 data found in the specified column after removing NAs.")
    return(NULL) # プロットせずにNULLを返す
  }
  
  cdr3_counts <- data.frame(cdr3_sequence = cdr3_filtered) %>%
    count(cdr3_sequence, sort = TRUE, name = "count")

  # --- 3. プロットするデータの整形 ---
  if (plot_by == "percentage") {
    plot_data <- cdr3_counts %>%
      mutate(value = (count / sum(count)) * 100) %>%
      slice_head(n = top_n)
    y_label <- "Percentage of Total TCR+ Cells (%)"
  } else {
    plot_data <- cdr3_counts %>%
      mutate(value = count) %>%
      slice_head(n = top_n)
    y_label <- "Cell Count"
  }
  
  # --- 4. ggplot2でプロット ---
  p <- ggplot(plot_data, aes(x = fct_reorder(cdr3_sequence, -value), y = value)) +
    geom_col(fill = plot_color, alpha = 0.9) +
    labs(
      title = paste("Top", top_n, "Clonotypes from", cdr3_column),
      x = "CDR3 Amino Acid Sequence",
      y = y_label
    ) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(p)
}
