# 必要ライブラリ (パッケージとしてビルドする場合の推奨記述)
# スクリプトとして直接実行する場合は、冒頭に library(dplyr), library(stringr), etc. の記述が必要です。
#' @importFrom dplyr %>% mutate case_when relocate starts_with coalesce across ends_with na_if select all_of any_of filter distinct full_join rename_with everything
#' @importFrom stringr str_c str_replace_all
#' @importFrom tidyselect everything
#' @importFrom readr read_csv
#' @importFrom scRepertoire combineBCR
#' @importFrom stats setNames

# --- 定数定義 ---
PREFIX_PAIR <- "BCR_pair_"
PREFIX_IGH <- "BCR_IGH_"
PREFIX_IGL <- "BCR_IGL_" # IGK/IGL共通のプレフィックスとして使用

# --- ヘルパー関数定義 ---

#' @title CSVからBCRペア鎖データを抽出・整形
#' @description 指定されたCSVファイルからscRepertoire::combineBCRを使用して
#'              ペア鎖データを抽出し、整形します。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @param sample_name `character(1)`. combineBCRで使用する一時的なサンプル名。
#'                    内部で除去されるため任意の値で良い。デフォルトは "temp_sample"。
#' @param id_name `character(1)`. combineBCRで使用する一時的なID名。
#'                 内部で除去されるため任意の値で良い。デフォルトは "temp_id"。
#' @return `data.frame`. 整形されたBCRペア鎖情報。
#'         列名は `PREFIX_PAIR` で指定されたプレフィックスが付与される。
csv_to_bcr_pair_dataframe <- function(csv_path,
                                      sample_name = "temp_sample",
                                      id_name = "temp_id") {
  
  # combineBCRはサンプル名とID名をbarcodeのプレフィックスに使う
  barcode_prefix <- paste0(sample_name, "_", id_name, "_")
  
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL) # エラー時はNULLを返す (通常はstopで停止)
  })
  
  if (is.null(bcr_raw)) return(NULL)
  
  # combineBCRはリストを返す仕様
  pair_list <- tryCatch({
    scRepertoire::combineBCR(bcr_raw, samples = sample_name, ID = id_name)
  }, error = function(e) {
    stop("Error in combineBCR: ", e$message)
    return(NULL)
  })
  
  # リストの最初の要素が存在するか確認
  if (is.null(pair_list) || length(pair_list) == 0 || is.null(pair_list[[1]])) {
    warning("combineBCR did not return valid pair data.")
    # 空のデータフレームを返すか、エラーにするか選択
    return(data.frame())
  }
  pair <- pair_list[[1]]
  
  pair <- pair %>%
    # combineBCRが付加したプレフィックスを除去
    dplyr::mutate(barcode = stringr::str_replace_all(barcode, pattern = barcode_prefix, replacement = "")) %>%
    # 不要な sample, ID 列を削除
    dplyr::select(-dplyr::any_of(c("sample", "ID"))) %>% # any_ofで列が存在しなくてもエラーにならない
    # 全ての列名にプレフィックスを追加
    dplyr::rename_with(~ stringr::str_c(PREFIX_PAIR, .), dplyr::everything())
  
  return(pair)
}

#' @title CSVからBCR重鎖(IGH)データを抽出・整形
#' @description 指定されたCSVファイルからIGH鎖のデータを抽出し、
#'              barcodeごとにユニークにして整形します。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 整形されたBCR IGH情報。
#'         列名は `PREFIX_IGH` で指定されたプレフィックスが付与される。
csv_to_bcr_igh_dataframe <- function(csv_path){
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  
  if (is.null(bcr_raw) || !"chain" %in% names(bcr_raw) || !"barcode" %in% names(bcr_raw)) {
    stop("CSV file must contain 'chain' and 'barcode' columns.")
    return(NULL)
  }
  
  IGH <- bcr_raw %>%
    dplyr::filter(chain == 'IGH') %>%
    # barcodeごとに最初の行のみを保持
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    # 全ての列名にプレフィックスを追加
    dplyr::rename_with(~ stringr::str_c(PREFIX_IGH, .), dplyr::everything())
  
  return(IGH)
}

#' @title CSVからBCR軽鎖(IGK/IGL)データを抽出・整形
#' @description 指定されたCSVファイルからIGKまたはIGL鎖のデータを抽出し、
#'              barcodeごとにユニークにして整形します。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 整形されたBCR IGK/IGL情報。
#'         列名は `PREFIX_IGL` で指定されたプレフィックスが付与される。
csv_to_bcr_igl_dataframe <- function(csv_path){
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  
  if (is.null(bcr_raw) || !"chain" %in% names(bcr_raw) || !"barcode" %in% names(bcr_raw)) {
    stop("CSV file must contain 'chain' and 'barcode' columns.")
    return(NULL)
  }
  
  # IGK または IGL 鎖を抽出
  IGL <- bcr_raw %>%
    dplyr::filter(chain %in% c("IGL", "IGK")) %>%
    # barcodeごとに最初の行のみを保持
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    # 全ての列名にプレフィックスを追加
    dplyr::rename_with(~ stringr::str_c(PREFIX_IGL, .), dplyr::everything())
  
  return(IGL)
}

#' @title BCR各鎖データを結合
#' @description ペア鎖、重鎖、軽鎖のデータフレームをbarcodeをキーにして結合します。
#'              `pair` データフレームのbarcode列を基準に結合します。
#' @param pair `data.frame`. `csv_to_bcr_pair_dataframe` の出力。
#' @param IGH `data.frame`. `csv_to_bcr_igh_dataframe` の出力。
#' @param IGL `data.frame`. `csv_to_bcr_igl_dataframe` の出力。
#' @return `data.frame`. 結合されたBCRデータフレーム。
merge_bcr <- function(pair, IGH, IGL){
  
  # barcode列名を特定 (各DFでプレフィックスが異なるため)
  pair_barcode_col <- paste0(PREFIX_PAIR, "barcode")
  igh_barcode_col <- paste0(PREFIX_IGH, "barcode")
  igl_barcode_col <- paste0(PREFIX_IGL, "barcode")
  
  # pairとIGHを結合
  # by句で結合キーの列名を明示的に指定 (左のDFのキー名 = 右のDFのキー名)
  # pairデータフレームに存在するbarcodeを基準とするため、left_joinが適切かもしれないが、
  # 元のコードに合わせてfull_joinを使用。
  # `setNames` は c("右のDFのキー名" = "左のDFのキー名") の形にする
  bcr_merged <- dplyr::full_join(pair, IGH, by = stats::setNames(igh_barcode_col, pair_barcode_col))
  
  # 上記結果とIGLを結合
  bcr_merged <- dplyr::full_join(bcr_merged, IGL, by = stats::setNames(igl_barcode_col, pair_barcode_col))
  
  # 注意: full_joinの結果、IGHやIGLのみに存在するbarcodeの行も含まれる可能性があります。
  # ペアとして存在するbarcodeのみに限定したい場合は、最初のpairデータフレームを基準に
  # left_joinを使うか、最後に filter(!is.na(.data[[pair_barcode_col]])) を行うことを検討します。
  # また、結合によって元のbarcode列（IGH_barcode, IGL_barcode）が残るため、
  # 必要に応じて select() で削除することも可能です。現状は維持されます。
  
  return(bcr_merged)
}


# --- メイン関数定義 ---

#' @title CSVからBCRデータを処理するメイン関数
#' @description 指定されたCSVファイルからBCRデータを読み込み、ペア鎖、重鎖、軽鎖に
#'              分割・整形した後、結合し、代表ID(raw_clonotype_idなど)や
#'              全長配列を生成し、列順序を整えます。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 全ての処理が完了したBCRデータフレーム。
bcr_csv_to_dataframe <- function(csv_path){
  
  # --- 入力チェック ---
  if (!file.exists(csv_path)) {
    stop("Input CSV file not found: ", csv_path)
    return(NULL)
  }
  
  # --- データ抽出・整形 ---
  pair <- csv_to_bcr_pair_dataframe(csv_path)
  IGH <- csv_to_bcr_igh_dataframe(csv_path)
  IGL <- csv_to_bcr_igl_dataframe(csv_path) # IGK/IGLを含む
  
  if (is.null(pair) || is.null(IGH) || is.null(IGL)) {
    stop("Failed to process one or more BCR components.")
    return(NULL)
  }
  
  # --- データ結合 ---
  df <- merge_bcr(pair, IGH, IGL)
  
  # --- 代表ID生成と不要列削除 ---
  # 列名を事前に定義
  igh_clonotype_col <- paste0(PREFIX_IGH, "raw_clonotype_id")
  igl_clonotype_col <- paste0(PREFIX_IGL, "raw_clonotype_id")
  igh_consensus_col <- paste0(PREFIX_IGH, "raw_consensus_id")
  igl_consensus_col <- paste0(PREFIX_IGL, "raw_consensus_id")
  igh_subclonotype_col <- paste0(PREFIX_IGH, "exact_subclonotype_id")
  igl_subclonotype_col <- paste0(PREFIX_IGL, "exact_subclonotype_id")
  
  # 削除対象の列名をベクトル化 (元のIGH/IGLのID列)
  cols_to_remove <- c(igh_clonotype_col, igl_clonotype_col,
                      igh_consensus_col, igl_consensus_col,
                      igh_subclonotype_col, igl_subclonotype_col)
  # データフレームに実際に存在する列のみを削除対象とする
  cols_to_remove_existing <- intersect(cols_to_remove, names(df))
  
  df <- df %>%
    dplyr::mutate(
      # 新しい代表ID列を対応するBCR_IGH_* 列から作成
      # 元の列が存在しない場合は NA を代入
      # .data[[]] を使用して列名を安全に参照
      raw_clonotype_id = if (igh_clonotype_col %in% names(.)) .data[[igh_clonotype_col]] else NA_character_,
      raw_consensus_id = if (igh_consensus_col %in% names(.)) .data[[igh_consensus_col]] else NA_character_,
      exact_subclonotype_id = if (igh_subclonotype_col %in% names(.)) .data[[igh_subclonotype_col]] else NA_character_
      # 注意: ここで以前あった BCR_pair_*_id の生成ロジックは削除されています（ユーザーの最初のコードブロックに基づき）
      # もしペアIDが必要な場合は、コメントアウトされた古いコードのロジックを追加してください。
    ) %>%
    # 不要なIGH/IGLの元ID列を削除
    dplyr::select(-dplyr::all_of(cols_to_remove_existing))
  
  # --- 列順序の整理 ---
  # 新しい代表ID列を先頭に配置し、次に各プレフィックスを持つ列を配置
  df <- df %>%
    dplyr::relocate(raw_clonotype_id, raw_consensus_id, exact_subclonotype_id, # 新しい代表ID列
                    dplyr::starts_with(PREFIX_PAIR),
                    dplyr::starts_with(PREFIX_IGH),
                    dplyr::starts_with(PREFIX_IGL),
                    .before = tidyselect::everything()) # それ以外の列の前に配置
  
  # --- 全長配列の生成 ---
  # 必要な列名を動的に生成
  regions_nt <- c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt")
  regions_aa <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
  
  igh_nt_cols <- paste0(PREFIX_IGH, regions_nt)
  igl_nt_cols <- paste0(PREFIX_IGL, regions_nt)
  igh_aa_cols <- paste0(PREFIX_IGH, regions_aa)
  igl_aa_cols <- paste0(PREFIX_IGL, regions_aa)
  
  # データフレームに存在する列のみを結合対象とする関数
  paste_existing_cols <- function(df, cols) {
    existing_cols <- intersect(cols, names(df))
    if (length(existing_cols) == 0) return(NA_character_)
    # 各行ごとに、存在する列の値を結合
    apply(df[, existing_cols, drop = FALSE], 1, function(row_values) {
      # NAを空文字に変換してから結合
      paste0(dplyr::coalesce(row_values, ""), collapse = "")
    })
  }
  
  df <- df %>%
    dplyr::mutate(
      # IGH ヌクレオチド全長配列
      !!paste0(PREFIX_IGH, "full_length_nt") := paste_existing_cols(., igh_nt_cols),
      # IGL/IGK ヌクレオチド全長配列
      !!paste0(PREFIX_IGL, "full_length_nt") := paste_existing_cols(., igl_nt_cols),
      # IGH アミノ酸全長配列
      !!paste0(PREFIX_IGH, "full_length_aa") := paste_existing_cols(., igh_aa_cols),
      # IGL/IGK アミノ酸全長配列
      !!paste0(PREFIX_IGL, "full_length_aa") := paste_existing_cols(., igl_aa_cols)
    ) %>%
    # もし生成された全長配列が空文字列になった場合（全ての構成要素がNAまたは存在しなかった場合）は NA に置換
    dplyr::mutate(
      dplyr::across(dplyr::ends_with("full_length_nt") | dplyr::ends_with("full_length_aa"), ~dplyr::na_if(., ""))
    )
  
  # --- 最後の列整理 (任意) ---
  # 例えば、結合によって生成されたIGHやIGLのbarcode列が不要であれば削除
  # df <- df %>% dplyr::select(-any_of(c(paste0(PREFIX_IGH, "barcode"), paste0(PREFIX_IGL, "barcode"))))
  # pairのbarcode列 (PREFIX_PAIR + "barcode") を "barcode" にリネームすることも可能
  # df <- df %>% dplyr::rename(barcode = !!paste0(PREFIX_PAIR, "barcode"))
  
  return(df)
}

# --- 関数の使用例 ---
# library(dplyr) # 必要に応じて実行前に読み込み
# library(stringr)
# library(readr)
# library(scRepertoire)
# library(tidyselect)

# csv_file <- "path/to/your/bcr_data.csv" # 実際のCSVファイルのパスに置き換えてください
# if (file.exists(csv_file)) {
#   processed_bcr_data <- bcr_csv_to_dataframe(csv_file)
#   print(head(processed_bcr_data))
#   # 必要に応じて結果を確認
#   # print(colnames(processed_bcr_data))
# } else {
#   warning("Example CSV file not found: ", csv_file)
# }
