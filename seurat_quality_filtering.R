# --- 必要なライブラリ ---
library(Seurat)
library(ggplot2)
library(patchwork)


# ==============================================================================
# 【完成版】Seuratオブジェクトをフィルタリングし、結果を可視化する関数
# - nFeature_RNA, nCount_RNA, percent.mt の3指標に対応！
# ==============================================================================
#' @title SeuratオブジェクトのQCフィルタリングと可視化
#' @description 3つのQC指標（nFeature_RNA, nCount_RNA, percent.mt）の
#'              閾値に基づいて細胞をフィルタリングし、処理前後のプロットを表示します。
#'
#' @param seurat_obj フィルタリング対象のSeuratオブジェクト。
#' @param min_features `nFeature_RNA`の下限値。
#' @param max_features `nFeature_RNA`の上限値。
#' @param min_count `nCount_RNA`の下限値。デフォルトは-Inf（下限なし）。
#' @param max_count `nCount_RNA`の上限値。デフォルトはInf（上限なし）。
#' @param max_percent_mt `percent.mt`の上限値。
#' @param show_plots `TRUE`の場合、フィルタリング前後のプロットを表示します。
#'
#' @return フィルタリング後のSeuratオブジェクト。
#'
seurat_quality_filtering <- function(seurat_obj, 
                               min_features = -Inf, 
                               max_features = Inf, 
                               min_count = -Inf, # ★nCount_RNAの下限を追加！
                               max_count = Inf,  # ★nCount_RNAの上限を追加！
                               max_percent_mt = 100,
                               show_plots = TRUE) {
  
  # フィルタリング前の細胞数を記録します
  n_cells_before <- ncol(seurat_obj)
  cat("--- フィルタリング前の細胞数:", n_cells_before, "---\n")
  
  # subset関数にnCount_RNAの条件も追加します！
  filtered_obj <- subset(seurat_obj, 
                         subset = nFeature_RNA > min_features & 
                                  nFeature_RNA < max_features & 
                                  nCount_RNA > min_count &      # ★条件を追加
                                  nCount_RNA < max_count &      # ★条件を追加
                                  percent.mt < max_percent_mt)
  
  # フィルタリング後の情報を表示します
  n_cells_after <- ncol(filtered_obj)
  cat("--- フィルタリング後の細胞数:", n_cells_after, "---\n")
  n_filtered_out <- n_cells_before - n_cells_after
  percent_filtered_out <- round((n_filtered_out / n_cells_before) * 100, 2)
  cat("--- 除外された細胞数:", n_filtered_out, paste0("(", percent_filtered_out, "%)"), "---\n")
  
  # プロット表示機能（ここは変更ありません）
  if (show_plots) {
    cat("--- フィルタリング前後のQCプロットを表示します ---\n\n")
    
    #
  
  return(filtered_obj)
}
