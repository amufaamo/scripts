# ==============================================================================
# BCR データ処理スクリプト
# ==============================================================================
#
# 説明:
#   指定されたCSVファイルからBCR (B Cell Receptor) データを読み込み、
#   ペア鎖、重鎖(IGH)、軽鎖(IGK/IGL)に分割・整形した後、結合します。
#   さらに、代表IDの生成、全長配列の生成、列順序の整理を行います。
#
# 依存ライブラリ:
#   dplyr, stringr, tidyselect, readr, scRepertoire
#
# 注意:
#   スクリプトとして直接実行する場合は、事前に以下のコマンドなどで
#   必要なライブラリをインストールし、スクリプト冒頭のコメントアウトされた
#   library() 関数を実行してください。
#   例: install.packages(c("dplyr", "stringr", "tidyselect", "readr", "devtools"))
#       devtools::install_github("ncborcherding/scRepertoire")
#
# ==============================================================================

# --- ライブラリ読み込み ---
# (スクリプトとして実行する場合)
# library(dplyr)
# library(stringr)
# library(tidyselect)
# library(readr)
# library(scRepertoire)

# (パッケージとしてビルドする場合の推奨記述)
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
  
  # CSV読み込み
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL) # エラー時はNULLを返す
  })
  
  if (is.null(bcr_raw)) return(NULL)
  
  # combineBCRでペア鎖情報を抽出 (リストで返る)
  pair_list <- tryCatch({
    scRepertoire::combineBCR(bcr_raw, samples = sample_name, ID = id_name)
  }, error = function(e) {
    stop("Error in combineBCR: ", e$message)
    return(NULL)
  })
  
  # 抽出結果のチェック
  if (is.null(pair_list) || length(pair_list) == 0 || is.null(pair_list[[1]])) {
    warning("combineBCR did not return valid pair data.")
    return(data.frame()) # 空のデータフレームを返す
  }
  pair <- pair_list[[1]]
  
  # 整形処理
  pair <- pair %>%
    # combineBCRが付加したプレフィックスを除去
    dplyr::mutate(barcode = stringr::str_replace_all(barcode, pattern = barcode_prefix, replacement = "")) %>%
    # 不要な sample, ID 列を削除
    dplyr::select(-dplyr::any_of(c("sample", "ID"))) %>% # 列が存在しなくてもエラーにならない
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
  
  # CSV読み込み
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  
  # 必須列のチェック
  if (is.null(bcr_raw) || !"chain" %in% names(bcr_raw) || !"barcode" %in% names(bcr_raw)) {
    stop("CSV file must contain 'chain' and 'barcode' columns.")
    return(NULL)
  }
  
  # IGH鎖の抽出と整形
  IGH <- bcr_raw %>%
    dplyr::filter(chain == 'IGH') %>%
    # barcodeごとに最初の行のみを保持 (重複除去)
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
  
  # CSV読み込み
  bcr_raw <- tryCatch({
    readr::read_csv(csv_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading CSV file: ", csv_path, "\n", e$message)
    return(NULL)
  })
  
  # 必須列のチェック
  if (is.null(bcr_raw) || !"chain" %in% names(bcr_raw) || !"barcode" %in% names(bcr_raw)) {
    stop("CSV file must contain 'chain' and 'barcode' columns.")
    return(NULL)
  }
  
  # IGK または IGL 鎖の抽出と整形
  IGL <- bcr_raw %>%
    dplyr::filter(chain %in% c("IGL", "IGK")) %>%
    # barcodeごとに最初の行のみを保持 (重複除去)
    dplyr::distinct(barcode, .keep_all = TRUE) %>%
    # 全ての列名にプレフィックスを追加
    dplyr::rename_with(~ stringr::str_c(PREFIX_IGL, .), dplyr::everything())
  
  return(IGL)
}


#' @title BCR各鎖データを結合
#' @description ペア鎖、重鎖、軽鎖のデータフレームをbarcodeをキーにして結合します。
#'              `pair` データフレームのbarcode列 (`BCR_pair_barcode`) を基準に結合します。
#' @param pair `data.frame`. `csv_to_bcr_pair_dataframe` の出力。
#' @param IGH `data.frame`. `csv_to_bcr_igh_dataframe` の出力。
#' @param IGL `data.frame`. `csv_to_bcr_igl_dataframe` の出力。
#' @return `data.frame`. 結合されたBCRデータフレーム。
merge_bcr <- function(pair, IGH, IGL){
  
  # barcode列名を特定 (各DFでプレフィックスが異なるため)
  pair_barcode_col <- paste0(PREFIX_PAIR, "barcode")
  igh_barcode_col <- paste0(PREFIX_IGH, "barcode")
  igl_barcode_col <- paste0(PREFIX_IGL, "barcode")
  
  # `pair` データフレームを基準に `IGH` を結合
  # `by` の指定: c("左のDFのキー名" = "右のDFのキー名")
  # `full_join` は両方のDFに存在するキーと、片方にしか存在しないキーの行をすべて保持
  bcr_merged <- dplyr::full_join(pair, IGH, by = stats::setNames(igh_barcode_col, pair_barcode_col))
  
  # 上記結果に `IGL` を結合
  bcr_merged <- dplyr::full_join(bcr_merged, IGL, by = stats::setNames(igl_barcode_col, pair_barcode_col))
  
  # 注意点:
  # - `full_join` を使用しているため、IGH鎖のみ、IGL鎖のみ、またはペア情報のみを持つbarcodeの行も結果に含まれる可能性があります。
  #   もしscRepertoire::combineBCRが出力したペア情報を持つbarcodeのみに関心がある場合は、
  #   `left_join` を使用するか、最後に `dplyr::filter(!is.na(.data[[pair_barcode_col]]))` でフィルタリングすることを検討してください。
  # - 結合により、元の `BCR_IGH_barcode`, `BCR_IGL_barcode` 列もデータフレームに残ります。不要であれば `select()` で除外できます。
  
  return(bcr_merged)
}


#' @title 存在する列の値を連結する内部ヘルパー関数
#' @description データフレームと列名のベクトルを受け取り、データフレーム内に
#'              存在する列の値を行ごとに連結します。存在しない列は無視されます。
#'              全ての指定列が存在しないか、値が全てNAの場合は空文字列を返します。
#' @param df `data.frame`. 対象のデータフレーム。
#' @param cols `character`. 連結したい列名のベクトル。
#' @return `character`. 各行について連結された文字列ベクトル。
paste_existing_cols <- function(df, cols) {
  # データフレーム内に存在する列名のみを抽出
  existing_cols <- intersect(cols, names(df))
  
  # 連結対象の列が一つも存在しない場合
  if (length(existing_cols) == 0) {
    return(rep(NA_character_, nrow(df))) # または rep("", nrow(df))
  }
  
  # 各行ごとに処理
  apply(df[, existing_cols, drop = FALSE], 1, function(row_values) {
    # NAを空文字列に置換してから連結
    paste0(dplyr::coalesce(row_values, ""), collapse = "")
  })
}


# --- メイン関数定義 ---

#' @title CSVからBCRデータを処理するメイン関数
#' @description 指定されたCSVファイルからBCRデータを読み込み、ペア鎖、重鎖、軽鎖に
#'              分割・整形した後、結合し、代表ID(raw_clonotype_idなど)や
#'              全長配列を生成し、列順序を整えます。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 全ての処理が完了したBCRデータフレーム。
bcr_csv_to_dataframe <- function(csv_path){
  
  # --- 1. 入力チェック ---
  if (!file.exists(csv_path)) {
    stop("Input CSV file not found: ", csv_path)
    return(NULL)
  }
  
  # --- 2. データ抽出・整形 ---
  # 各コンポーネント (ペア、重鎖、軽鎖) のデータを抽出・整形
  pair_data <- csv_to_bcr_pair_dataframe(csv_path)
  igh_data <- csv_to_bcr_igh_dataframe(csv_path)
  igl_data <- csv_to_bcr_igl_dataframe(csv_path) # IGK/IGLを含む
  
  # いずれかの処理でエラーが発生した場合 (NULLが返された場合) は停止
  if (is.null(pair_data) || is.null(igh_data) || is.null(igl_data)) {
    stop("Failed to process one or more BCR components. Check previous errors or warnings.")
    return(NULL)
  }
  
  # --- 3. データ結合 ---
  # 抽出したデータをbarcodeをキーに結合
  merged_df <- merge_bcr(pair_data, igh_data, igl_data)
  
  # --- 4. 代表ID生成と不要列削除 ---
  # 代表IDとして使用する元の列名を定義
  igh_clonotype_col <- paste0(PREFIX_IGH, "raw_clonotype_id")
  igh_consensus_col <- paste0(PREFIX_IGH, "raw_consensus_id")
  igh_subclonotype_col <- paste0(PREFIX_IGH, "exact_subclonotype_id")
  
  # 代表ID生成後に削除する元のID列名を定義 (IGHとIGL/IGKの両方)
  igl_clonotype_col <- paste0(PREFIX_IGL, "raw_clonotype_id")
  igl_consensus_col <- paste0(PREFIX_IGL, "raw_consensus_id")
  igl_subclonotype_col <- paste0(PREFIX_IGL, "exact_subclonotype_id")
  
  cols_to_remove <- c(
    igh_clonotype_col, igl_clonotype_col,
    igh_consensus_col, igl_consensus_col,
    igh_subclonotype_col, igl_subclonotype_col
  )
  # データフレームに実際に存在する列のみを削除対象とする
  cols_to_remove_existing <- intersect(cols_to_remove, names(merged_df))
  
  # 代表ID列を生成し、元のID列を削除
  processed_df <- merged_df %>%
    dplyr::mutate(
      # 新しい代表ID列をBCR_IGH_* 列から作成 (存在しない場合はNA)
      # .data[[]] を使用して列名を安全に参照
      raw_clonotype_id = if (igh_clonotype_col %in% names(.)) .data[[igh_clonotype_col]] else NA_character_,
      raw_consensus_id = if (igh_consensus_col %in% names(.)) .data[[igh_consensus_col]] else NA_character_,
      exact_subclonotype_id = if (igh_subclonotype_col %in% names(.)) .data[[igh_subclonotype_col]] else NA_character_
      # 注意: 以前のコードにあった BCR_pair_*_id (IGHとIGLのIDを結合) の生成ロジックは、
      #       ユーザー提示の最初のコードに基づき、ここでは実装されていません。
      #       もしペアIDが必要な場合は、別途 mutate() で追加してください。
      #       例:
      #       pair_raw_clonotype_id = dplyr::case_when(
      #         !is.na(.data[[igh_clonotype_col]]) & !is.na(.data[[igl_clonotype_col]])
      #         ~ stringr::str_c(.data[[igh_clonotype_col]], .data[[igl_clonotype_col]], sep = "_"),
      #         TRUE ~ NA_character_
      #       ),
      #       ... (他のペアIDも同様) ...
    ) %>%
    # 不要になった元のIGH/IGLのID列を削除
    dplyr::select(-dplyr::all_of(cols_to_remove_existing))
  
  # --- 5. 列順序の整理 ---
  # 代表ID列、ペア鎖列、重鎖列、軽鎖列の順に並べ替え
  processed_df <- processed_df %>%
    dplyr::relocate(
      raw_clonotype_id, raw_consensus_id, exact_subclonotype_id, # 新しい代表ID列
      dplyr::starts_with(PREFIX_PAIR),                          # ペア鎖関連列
      dplyr::starts_with(PREFIX_IGH),                           # 重鎖関連列
      dplyr::starts_with(PREFIX_IGL),                           # 軽鎖関連列
      .before = tidyselect::everything()                        # その他の列の前に配置
    )
  
  # --- 6. 全長配列の生成 ---
  # ヌクレオチド配列とアミノ酸配列のリージョン名を定義
  regions_nt <- c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt")
  regions_aa <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
  
  # 各鎖・各配列タイプの列名を生成
  igh_nt_cols <- paste0(PREFIX_IGH, regions_nt)
  igl_nt_cols <- paste0(PREFIX_IGL, regions_nt)
  igh_aa_cols <- paste0(PREFIX_IGH, regions_aa)
  igl_aa_cols <- paste0(PREFIX_IGL, regions_aa)
  
  # 全長配列を生成 (内部関数 paste_existing_cols を使用)
  processed_df <- processed_df %>%
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
    # 生成された全長配列が空文字列の場合 (元の列が存在しない or 全てNA) は NA に置換
    dplyr::mutate(
      dplyr::across(
        dplyr::ends_with("full_length_nt") | dplyr::ends_with("full_length_aa"),
        ~ dplyr::na_if(., "")
      )
    )
  
  # --- 7. 最終的な列整理 (任意) ---
  # 例えば、結合によって残ったIGHやIGLのbarcode列が不要な場合:
  # processed_df <- processed_df %>%
  #   dplyr::select(-dplyr::any_of(c(paste0(PREFIX_IGH, "barcode"), paste0(PREFIX_IGL, "barcode"))))
  
  # ペアのbarcode列 (例: BCR_pair_barcode) を単に "barcode" にリネームする場合:
  # pair_barcode_col <- paste0(PREFIX_PAIR, "barcode")
  # if (pair_barcode_col %in% names(processed_df)) {
  #   processed_df <- processed_df %>%
  #     dplyr::rename(barcode = !!rlang::sym(pair_barcode_col))
  # }
  
  # --- 8. 結果を返す ---
  return(processed_df)
}

# --- 関数の使用例 ---
# (必要なライブラリを読み込んだ上で実行)

# # 解析対象のCSVファイルのパスを指定
# csv_file <- "path/to/your/bcr_data.csv" # <--- 実際のパスに置き換えてください

# # ファイルが存在するか確認してから実行
# if (file.exists(csv_file)) {
#     # メイン関数を実行してBCRデータを処理
#     processed_bcr_data <- bcr_csv_to_dataframe(csv_file)

#     # 結果の最初の数行を表示
#     print(head(processed_bcr_data))

#     # 結果の列名を確認 (任意)
#     # print(colnames(processed_bcr_data))

#     # 結果の次元数を確認 (任意)
#     # print(dim(processed_bcr_data))

# } else {
#     warning("Example CSV file not found: ", csv_file)
# }
