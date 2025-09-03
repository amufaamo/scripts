# --- 最初に、必要なライブラリを読み込んでおいてくださいね！ ---
# まだインストールしていない場合は install.packages("Seurat") などでインストールしてください。
library(Seurat)
library(ggplot2)
library(patchwork)


# ==============================================================================
# 10X GenomicsのH5データからSeuratオブジェクトを作成しQCを行う関数
# ==============================================================================
#' @title Seuratオブジェクトの作成とQC
#' @description 10X GenomicsのH5ファイルを読み込み、Seuratオブジェクトを作成し、
#'              基本的なQCメトリクスを計算してプロットを表示します。
#'
#' @param h5_file_path (character) 入力する`filtered_feature_bc_matrix.h5`ファイルへのパスです。
#' @param project_name (character) プロジェクト名（Seuratオブジェクトに格納されます）。
#' @param save_plots (logical) `TRUE`の場合、QCプロットをPDFファイルとして保存します。デフォルトは`FALSE`です。
#' @param output_dir (character) プロットを保存する場合の出力先ディレクトリです。デフォルトは現在の作業ディレクトリです。
#'
#' @return QC情報（ミトコンドリア遺伝子割合など）が追加されたSeuratオブジェクトを返します。
#'
seurat_h5_qc <- function(h5_file_path, 
                              project_name, 
                              save_plots = FALSE, 
                              output_dir = ".") {
  
  # --- 1. 入力ファイルの存在をチェックします ---
  if (!file.exists(h5_file_path)) {
    stop("指定されたH5ファイルが見つかりません: ", h5_file_path)
  }
  cat("--- ファイルを読み込みます: ", h5_file_path, "---\n")

  # --- 2. データを読み込み、Seuratオブジェクトを作成します ---
  input_data <- Read10X_h5(filename = h5_file_path)
  
  # マルチモーダルデータの場合は遺伝子発現データを抽出
  if (inherits(input_data, "list") && "Gene Expression" %in% names(input_data)) {
    counts <- input_data[["Gene Expression"]]
  } else {
    counts <- input_data
  }
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = project_name)
  cat("--- Seuratオブジェクトを作成しました ---\n")
  
  # --- 3. ミトコンドリア遺伝子の割合を計算します ---
  # ヒトの場合は "^MT-"、マウスの場合は "^mt-" を使ってくださいね！
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  cat("--- QCメトリクスを計算しました ---\n")
  
  # --- 4. QCプロットを作成します ---
  cat("--- QCプロットを作成しています... ---\n")
  vln_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  
  scatter_plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter_plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  # --- 5. プロットを表示し、オプションで保存もします ---
  # RStudioのPlotsペインにプロットを表示！
  print(vln_plot)
  print(scatter_plot1 + scatter_plot2)
  
  # もし save_plots = TRUE なら、PDFにも保存します
  if (save_plots) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    pdf_path <- file.path(output_dir, paste0(project_name, "_QC_plots.pdf"))
    
    pdf(pdf_path, width = 12, height = 5)
    print(vln_plot)
    print(scatter_plot1 + scatter_plot2)
    dev.off() # PDFファイルを閉じます
    
    cat("--- プロットをPDFに保存しました:", pdf_path, "---\n")
  }
  
  # --- 6. 処理が完了したSeuratオブジェクトを返します ---
  cat("--- 処理が完了しました！ ---\n\n")
  return(seurat_obj)
}
