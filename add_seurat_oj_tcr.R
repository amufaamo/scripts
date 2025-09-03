# 必要なライブラリを最初に読み込んでください
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)
library(scRepertoire) # combineTCR を使うので必要です
library(readr)

# ==============================================================================
# 【自己完結版】SeuratオブジェクトにTCR情報を追加する最終関数
# ==============================================================================
#' @title SeuratオブジェクトへのTCR情報の統合（自己完結版）
#' @description Cell RangerのTCRアノテーションファイル（CSV）を読み込み、処理し、
#'              指定されたSeuratオブジェクトのメタデータに結合します。
#'              必要なヘルパー関数は全てこの関数内に含まれています。
#'
#' @param seurat_object TCR情報を追加したいSeuratオブジェクト。
#' @param tcr_csv_path (character) `filtered_contig_annotations.csv` ファイルへのパス。
#'
#' @return TCR情報がメタデータに追加された、更新済みのSeuratオブジェクト。
#'
add_seurat_oj_tcr <- function(seurat_object,
                               tcr_csv_path,
                               filter_multi_chain = TRUE) { # ★新しい引数を追加！

  # --- (内部ヘルパー関数の定義は長いので、ここでは省略しますが、前の回答と同じものが含まれます) ---
  # --- START: 内部ヘルパー関数の定義 ---
  PREFIX_TCR_PAIR  <- "TCR_pair_"; PREFIX_TCR_ALPHA <- "TCR_TRA_"; PREFIX_TCR_BETA  <- "TCR_TRB_"
  paste_existing_cols <- function(df, cols) {
    existing_cols <- intersect(cols, names(df)); if (length(existing_cols) == 0) { return(rep(NA_character_, nrow(df))) }
    apply(df[, existing_cols, drop = FALSE], 1, function(row_values) { paste0(ifelse(is.na(row_values), "", row_values), collapse = "") })
  }
  csv_to_tcr_pair_dataframe <- function(csv_path, sample_name = "temp_sample", id_name = "temp_id") {
    barcode_prefix <- paste0(sample_name, "_", id_name, "_")
    tcr_raw <- tryCatch({ readr::read_csv(csv_path, show_col_types = FALSE) }, error = function(e) { stop(e$message); NULL })
    if (is.null(tcr_raw)) return(NULL)
    # ★★★ ここで新しい引数を使います！ ★★★
    pair_list <- tryCatch({ scRepertoire::combineTCR(tcr_raw, samples = sample_name, ID = id_name, filterMulti = filter_multi_chain, removeNA = TRUE) }, error = function(e) { stop(e$message); NULL })
    if (is.null(pair_list) || length(pair_list) == 0 || is.null(pair_list[[1]])) { return(data.frame()) }
    pair <- pair_list[[1]]
    pair <- pair %>% dplyr::mutate(barcode = stringr::str_replace_all(barcode, pattern = barcode_prefix, replacement = "")) %>% dplyr::select(-dplyr::any_of(c("sample", "ID"))) %>% dplyr::rename_with(~ stringr::str_c(PREFIX_TCR_PAIR, .), dplyr::everything())
    return(pair)
  }
  csv_to_tcr_tra_dataframe <- function(csv_path){
    tcr_raw <- tryCatch({ readr::read_csv(csv_path, show_col_types = FALSE) }, error = function(e) { stop(e$message); NULL }); if (is.null(tcr_raw) || !"chain" %in% names(tcr_raw)) { return(data.frame()) }
    TRA <- tcr_raw %>% dplyr::filter(chain == 'TRA') %>% dplyr::distinct(barcode, .keep_all = TRUE) %>% dplyr::rename_with(~ stringr::str_c(PREFIX_TCR_ALPHA, .), dplyr::everything()); return(TRA)
  }
  csv_to_tcr_trb_dataframe <- function(csv_path){
    tcr_raw <- tryCatch({ readr::read_csv(csv_path, show_col_types = FALSE) }, error = function(e) { stop(e$message); NULL }); if (is.null(tcr_raw) || !"chain" %in% names(tcr_raw)) { return(data.frame()) }
    TRB <- tcr_raw %>% dplyr::filter(chain == 'TRB') %>% dplyr::distinct(barcode, .keep_all = TRUE) %>% dplyr::rename_with(~ stringr::str_c(PREFIX_TCR_BETA, .), dplyr::everything()); return(TRB)
  }
  merge_tcr <- function(pair, TRA, TRB){
    pair_barcode_col <- paste0(PREFIX_TCR_PAIR, "barcode"); tra_barcode_col  <- paste0(PREFIX_TCR_ALPHA, "barcode"); trb_barcode_col  <- paste0(PREFIX_TCR_BETA, "barcode"); if (nrow(pair) == 0) { return(data.frame()) }
    tcr_merged <- dplyr::full_join(pair, TRA, by = stats::setNames(tra_barcode_col, pair_barcode_col)); tcr_merged <- dplyr::full_join(tcr_merged, TRB, by = stats::setNames(trb_barcode_col, pair_barcode_col)); return(tcr_merged)
  }
  tcr_csv_to_dataframe <- function(csv_path){
    pair_data <- csv_to_tcr_pair_dataframe(csv_path); tra_data  <- csv_to_tcr_tra_dataframe(csv_path); trb_data  <- csv_to_tcr_trb_dataframe(csv_path); if (nrow(pair_data) == 0) { return(data.frame()) }
    merged_df <- merge_tcr(pair_data, tra_data, trb_data); if (nrow(merged_df) == 0) { return(merged_df) }
    trb_clonotype_col <- paste0(PREFIX_TCR_BETA, "raw_clonotype_id"); trb_subclonotype_col <- paste0(PREFIX_TCR_BETA, "exact_subclonotype_id"); cols_to_remove <- c(paste0(PREFIX_TCR_ALPHA, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")), paste0(PREFIX_TCR_BETA, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")), paste0(PREFIX_TCR_PAIR, c("raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id"))); cols_to_remove_existing <- intersect(cols_to_remove, names(merged_df))
    processed_df <- merged_df %>% dplyr::mutate(raw_clonotype_id = if (trb_clonotype_col %in% names(.)) .data[[trb_clonotype_col]] else NA_character_, exact_subclonotype_id = dplyr::case_when(trb_clonotype_col %in% names(.) & trb_subclonotype_col %in% names(.) & !is.na(.data[[trb_clonotype_col]]) & !is.na(.data[[trb_subclonotype_col]]) ~ stringr::str_c(.data[[trb_clonotype_col]], .data[[trb_subclonotype_col]], sep = "_"), TRUE ~ NA_character_)) %>% dplyr::select(-dplyr::any_of(cols_to_remove_existing), -dplyr::any_of("raw_consensus_id")) %>% dplyr::relocate(dplyr::any_of(c("raw_clonotype_id", "exact_subclonotype_id")), dplyr::starts_with(PREFIX_TCR_PAIR), dplyr::starts_with(PREFIX_TCR_ALPHA), dplyr::starts_with(PREFIX_TCR_BETA), .before = tidyselect::everything())
    regions_nt <- c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt"); regions_aa <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
    processed_df <- processed_df %>% dplyr::mutate(!!paste0(PREFIX_TCR_ALPHA, "full_length_nt") := paste_existing_cols(., paste0(PREFIX_TCR_ALPHA, regions_nt)), !!paste0(PREFIX_TCR_BETA, "full_length_nt")  := paste_existing_cols(., paste0(PREFIX_TCR_BETA, regions_nt)), !!paste0(PREFIX_TCR_ALPHA, "full_length_aa") := paste_existing_cols(., paste0(PREFIX_TCR_ALPHA, regions_aa)), !!paste0(PREFIX_TCR_BETA, "full_length_aa")  := paste_existing_cols(., paste0(PREFIX_TCR_BETA, regions_aa))) %>% dplyr::mutate(dplyr::across(dplyr::ends_with("full_length_nt") | dplyr::ends_with("full_length_aa"), ~ dplyr::na_if(., "")))
    pair_barcode_col_rename <- paste0(PREFIX_TCR_PAIR, "barcode"); if (!"barcode" %in% names(processed_df) && pair_barcode_col_rename %in% names(processed_df)) { processed_df <- processed_df %>% dplyr::rename(barcode = !!rlang::sym(pair_barcode_col_rename)) }; if ("barcode" %in% names(processed_df)) { processed_df <- processed_df %>% dplyr::relocate(barcode) }
    tra_key_col <- paste0(PREFIX_TCR_ALPHA, "chain"); trb_key_col <- paste0(PREFIX_TCR_BETA, "chain"); if (tra_key_col %in% names(processed_df) && trb_key_col %in% names(processed_df)) { final_df <- processed_df %>% dplyr::filter(!is.na(.data[[tra_key_col]]) & !is.na(.data[[trb_key_col]])) } else { final_df <- processed_df }; return(final_df)
  }
  # --- END: 内部ヘルパー関数の定義 ---

  ### ここからがメインの実行ロジックです ###
  if (!inherits(seurat_object, "Seurat")) stop("Error: 'seurat_object' must be a Seurat object.")
  if (!file.exists(tcr_csv_path)) stop("Error: TCR CSV file not found at: ", tcr_csv_path)
  cat("--- Step 1: TCRアノテーションファイルを読み込んで整形します ---\n")
  tcr_df <- tcr_csv_to_dataframe(tcr_csv_path)
  if (nrow(tcr_df) == 0) { warning("Warning: TCR data frame is empty. No TCR info will be added."); return(seurat_object) }
  cat("--- Step 2: バーコードをマッチングのために準備します ---\n")
  seurat_metadata <- seurat_object@meta.data %>% tibble::rownames_to_column(var = "seurat_barcode") %>% mutate(barcode_key = str_replace(seurat_barcode, "-[0-9]+$", ""))
  tcr_df_prepared <- tcr_df %>% mutate(barcode_key = str_replace(barcode, "-[0-9]+$", "")) %>% distinct(barcode_key, .keep_all = TRUE)
  cat(paste0(">>> Seuratオブジェクトの細胞数: ", nrow(seurat_metadata), "\n")); cat(paste0(">>> 読み込んだTCRエントリー数: ", nrow(tcr_df_prepared), "\n"))
  cat("--- Step 3: メタデータにTCR情報を左結合（left_join）します ---\n")
  combined_metadata <- dplyr::left_join(seurat_metadata, tcr_df_prepared, by = "barcode_key")
  if(nrow(combined_metadata) != nrow(seurat_metadata)) stop("Error: Merging metadata changed the number of cells. Aborting.")
  rownames(combined_metadata) <- combined_metadata$seurat_barcode; combined_metadata <- combined_metadata %>% select(-seurat_barcode, -barcode_key)
  cat("--- Step 4: Seuratオブジェクトのメタデータを更新します ---\n")
  seurat_object@meta.data <- combined_metadata
  match_col <- paste0(PREFIX_TCR_PAIR, "CTgene"); n_matched <- if (match_col %in% names(seurat_object@meta.data)) sum(!is.na(seurat_object@meta.data[[match_col]])) else 0
  percent_matched <- round((n_matched / ncol(seurat_object)) * 100, 2)
  cat(paste0("--- 完了！ ", n_matched, "個の細胞 (", percent_matched, "%) にTCR情報が追加されました ---\n"))
  return(seurat_object)
}
