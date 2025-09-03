# 必要なライブラリ
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)
library(readr)

# ==============================================================================
# 【最終版・デバッグ情報強化】TCRペアを手動で作成し、情報を追加する関数
# ==============================================================================
add_seurat_oj_tcr <- function(seurat_object, tcr_csv_path) {

  # (内部のヘルパー関数定義は変更なし)
  # --- START: 内部ヘルパー関数の定義 ---
  tcr_csv_to_dataframe <- function(csv_path) {
    tcr_raw <- readr::read_csv(csv_path, show_col_types = FALSE)
    tcr_productive <- tcr_raw %>% mutate(productive = (tolower(productive) == "true")) %>% filter(productive == TRUE)
    tra_representative <- tcr_productive %>% filter(chain == "TRA") %>% group_by(barcode) %>% slice_max(order_by = umis, n = 1, with_ties = FALSE) %>% ungroup() %>% select(barcode, TRA_v_gene = v_gene, TRA_j_gene = j_gene, TRA_cdr3 = cdr3, TRA_cdr3_nt = cdr3_nt, TRA_umis = umis)
    trb_representative <- tcr_productive %>% filter(chain == "TRB") %>% group_by(barcode) %>% slice_max(order_by = umis, n = 1, with_ties = FALSE) %>% ungroup() %>% select(barcode, TRB_v_gene = v_gene, TRB_j_gene = j_gene, TRB_cdr3 = cdr3, TRB_cdr3_nt = cdr3_nt, TRB_umis = umis)
    tcr_paired_manual <- dplyr::inner_join(tra_representative, trb_representative, by = "barcode")
    return(tcr_paired_manual)
  }
  # --- END: 内部ヘルパー関数の定義 ---

  ### ここからがメインの実行ロジックです ###
  if (!inherits(seurat_object, "Seurat")) stop("Error: 'seurat_object' must be a Seurat object.")
  if (!file.exists(tcr_csv_path)) stop("Error: TCR CSV file not found at: ", tcr_csv_path)

  cat("--- Step 1: TCRデータを手動でペアリングします ---\n")
  tcr_df <- tcr_csv_to_dataframe(tcr_csv_path)
  cat(">>> 作成されたTCRペアの数:", nrow(tcr_df), "\n")
  if (nrow(tcr_df) == 0) { warning("Warning: 有効なTCRペアが見つかりませんでした。"); return(seurat_object) }

  cat("--- Step 2: SeuratとTCRの共通バーコードを検索します ---\n")
  seurat_barcodes <- str_replace(colnames(seurat_object), "-[0-9]+$", "")
  tcr_barcodes <- str_replace(tcr_df$barcode, "-[0-9]+$", "")
  common_barcodes <- intersect(seurat_barcodes, tcr_barcodes)
  cat(">>> 共通のバーコード数:", length(common_barcodes), "\n")
  if (length(common_barcodes) == 0) { warning("Warning: SeuratとTCRデータで共通のバーコードが見つかりませんでした。"); return(seurat_object) }

  cat("--- Step 3: メタデータにTCR情報を結合します ---\n")
  seurat_metadata <- seurat_object@meta.data %>%
    tibble::rownames_to_column(var = "seurat_barcode") %>%
    mutate(barcode_key = str_replace(seurat_barcode, "-[0-9]+$", ""))
  tcr_to_join <- tcr_df %>%
    mutate(barcode_key = str_replace(barcode, "-[0-9]+$", ""))
  combined_metadata <- dplyr::left_join(seurat_metadata, tcr_to_join, by = "barcode_key")
  rownames(combined_metadata) <- combined_metadata$seurat_barcode
  
  cat("--- Step 4: Seuratオブジェクトを更新し、結果を表示します ---\n")
  seurat_object@meta.data <- combined_metadata %>% select(-seurat_barcode, -barcode_key, -barcode)

  # 成功した行を抽出
  successful_rows <- seurat_object@meta.data %>% filter(!is.na(TRA_v_gene))
  n_matched <- nrow(successful_rows)
  percent_matched <- round((n_matched / ncol(seurat_object)) * 100, 2)
  
  cat(paste0("--- 完了！🎉 ", n_matched, "個の細胞 (", percent_matched, "%) にTCR情報が正常に追加されました ---\n\n"))
  cat(">>> 追加されたTCR情報のプレビュー（最初の数件）:\n")
  print(head(successful_rows))

  return(seurat_object)
}
