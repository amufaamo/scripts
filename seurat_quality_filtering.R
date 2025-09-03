# --- 必要なライブラリを再度確認！ ---
# patchworkはプロットを素敵に並べるのに使います
library(Seurat)
library(ggplot2)
library(patchwork)


# ==============================================================================
# 【改良版】Seuratオブジェクトをフィルタリングし、結果を可視化する関数
# ==============================================================================
#' @title SeuratオブジェクトのQCフィルタリングと可視化
#' @description nFeature_RNAとpercent.mtの閾値で細胞をフィルタリングし、
#'              処理前後のQCプロットを並べて表示します。
#'
#' @param seurat_obj フィルタリング対象のSeuratオブジェクト。
#' @param min_features `nFeature_RNA`の下限値。
#' @param max_features `nFeature_RNA`の上限値。
#' @param max_percent_mt `percent.mt`の上限値。
#' @param show_plots `TRUE`の場合、フィルタリング前後のプロットを表示します。
#'
#' @return フィルタリング後のSeuratオブジェクト。
#'
seurat_quality_filtering <- function(seurat_obj, 
                               min_features = 200, 
                               max_features = Inf, 
                               max_percent_mt = 5,
                               show_plots = TRUE) {
  
  # フィルタリング前の情報を記録します
  n_cells_before <- ncol(seurat_obj)
  cat("--- フィルタリング前の細胞数:", n_cells_before, "---\n")
  
  # subset関数でフィルタリングを実行します
  filtered_obj <- subset(seurat_obj, 
                         subset = nFeature_RNA > min_features & 
                                  nFeature_RNA < max_features & 
                                  percent.mt < max_percent_mt)
  
  # フィルタリング後の情報を表示します
  n_cells_after <- ncol(filtered_obj)
  cat("--- フィルタリング後の細胞数:", n_cells_after, "---\n")
  n_filtered_out <- n_cells_before - n_cells_after
  percent_filtered_out <- round((n_filtered_out / n_cells_before) * 100, 2)
  cat("--- 除外された細胞数:", n_filtered_out, paste0("(", percent_filtered_out, "%)"), "---\n")
  
  # show_plotsがTRUEの場合、プロットを表示します
  if (show_plots) {
    cat("--- フィルタリング前後のQCプロットを表示します ---\n\n")
    
    # 1. ヴァイオリンプロットの比較
    p_vln_before <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) + 
                    patchwork::plot_annotation(title = "Before Filtering")
    
    p_vln_after <- VlnPlot(filtered_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) + 
                   patchwork::plot_annotation(title = "After Filtering")
    
    # patchworkを使って上下に並べて表示！
    print(p_vln_before / p_vln_after)
    
    # 2. 散布図の比較
    # フィルタリング前の全データをグレーで、フィルタリング後のデータをカラーで重ねてプロットします
    # これで、どの点が除外されたか一目瞭然です！
    plot_data_before <- seurat_obj@meta.data
    plot_data_after <- filtered_obj@meta.data
    
    p_scatter1 <- ggplot(plot_data_before, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(alpha = 0.4, color = "grey80", size = 0.5) + # 全ての点を薄いグレーで描画
      geom_point(data = plot_data_after, alpha = 0.4, color = "#F8766D", size = 0.5) + # 残った点だけを色付きで上描き
      labs(title = "nFeature_RNA vs nCount_RNA") +
      theme_classic(base_size = 14)

    p_scatter2 <- ggplot(plot_data_before, aes(x = nCount_RNA, y = percent.mt)) +
      geom_point(alpha = 0.4, color = "grey80", size = 0.5) +
      geom_point(data = plot_data_after, alpha = 0.4, color = "#F8766D", size = 0.5) +
      labs(title = "percent.mt vs nCount_RNA") +
      theme_classic(base_size = 14)
      
    print(p_scatter1 + p_scatter2)
  }
  
  return(filtered_obj)
}
