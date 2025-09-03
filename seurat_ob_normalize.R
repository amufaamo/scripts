# --- 必要なライブラリの確認 ---
library(Seurat)
library(ggplot2)
library(patchwork)
library(glmGamPoi)


# ==============================================================================
# Seurat解析パイプライン（正規化〜クラスタリング）を実行する関数
# ==============================================================================
#' @title Seurat標準解析パイプラインの実行
#' @description フィルタリング済みのSeuratオブジェクトに対して、SCTransform、PCA、
#'              クラスタリング、UMAP可視化までの一連の解析を実行します。
#'
#' @param seurat_obj フィルタリング済みのSeuratオブジェクト。
#' @param dims_to_use (numeric vector) `FindNeighbors`と`RunUMAP`で使用するPCの次元数。
#'                    例: `1:20`。ElbowPlotを見て決めるのがおすすめです。
#' @param cluster_resolution (numeric) `FindClusters`で使用する解像度。
#'                           値が大きいほどクラスター数は多くなります。
#'
#' @return 解析が完了し、UMAPとクラスター情報が追加されたSeuratオブジェクト。
#'
seurat_ob_normalize <- function(seurat_obj, 
                              dims_to_use = 1:20, 
                              cluster_resolution = 0.5) {
  
  cat("--- Seurat解析パイプラインを開始します ---\n")
  
  # --- Step 1: SCTransformによる正規化 ---
  cat("Step 1: SCTransformによる正規化を実行中...\n")
  seurat_obj <- SCTransform(seurat_obj, 
                            method = "glmGamPoi", 
                            vars.to.regress = "percent.mt", 
                            verbose = FALSE)
  
  # --- Step 2: PCAと次元の決定 ---
  cat("Step 2: PCAを実行し、ElbowPlotを表示します...\n")
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  
  # ElbowPlotを表示して、ユーザーがdims_to_useの参考にできるようにします
  print(ElbowPlot(seurat_obj, ndims = 50))
  cat(">>> ElbowPlotを確認し、必要であればdims_to_useの引数を調整してください。\n")
  
  # --- Step 3: クラスタリング ---
  cat(paste0("Step 3: ", length(dims_to_use), "次元、解像度", cluster_resolution, "でクラスタリングを実行中...\n"))
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_to_use, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = cluster_resolution, verbose = FALSE)
  
  # --- Step 4: UMAPによる可視化 ---
  cat("Step 4: UMAPを実行し、結果をプロットします...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_to_use, verbose = FALSE)
  
  # 最終的なUMAPプロットを表示
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE))
  
  cat("--- 解析パイプラインが正常に完了しました！ ---\n")
  
  return(seurat_obj)
}
