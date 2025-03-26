#' @title CSVからBCRデータを処理するメイン関数
#' @description 指定されたCSVファイルからBCRデータを読み込み、ペア鎖、重鎖、軽鎖に
#'              分割・整形した後、結合し、代表ID(raw_clonotype_idなど)や
#'              全長配列を生成し、列順序を整えます。
#' @param csv_path `character(1)`. 入力CSVファイルのパス。
#' @return `data.frame`. 全ての処理が完了したBCRデータフレーム。
#' @importFrom dplyr %>% mutate case_when relocate starts_with coalesce across ends_with na_if select all_of
#' @importFrom stringr str_c
#' @importFrom tidyselect everything
bcr_csv_to_dataframe <- function(csv_path){
  
  # --- 入力チェック ---
  if (!file.exists(csv_path)) {
    stop("Input CSV file not found: ", csv_path)
    return(NULL)
  }
  
  # --- データ抽出・整形 ---
  pair <- csv_to_bcr_pair_dataframe(csv_path)
  IGH <- csv_to_bcr_igh_dataframe(csv_path)
  IGL <- csv_to_bcr_igl_dataframe(csv_path)
  
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
      raw_clonotype_id = if (igh_clonotype_col %in% names(.)) .data[[igh_clonotype_col]] else NA_character_,
      raw_consensus_id = if (igh_consensus_col %in% names(.)) .data[[igh_consensus_col]] else NA_character_,
      exact_subclonotype_id = if (igh_subclonotype_col %in% names(.)) .data[[igh_subclonotype_col]] else NA_character_
      # 注意: ここで以前あった BCR_pair_*_id の生成ロジックは削除されています
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
  # 必要な列名を定義
  igh_fwr1_nt_col <- paste0(PREFIX_IGH, "fwr1_nt")
  igh_cdr1_nt_col <- paste0(PREFIX_IGH, "cdr1_nt")
  igh_fwr2_nt_col <- paste0(PREFIX_IGH, "fwr2_nt")
  igh_cdr2_nt_col <- paste0(PREFIX_IGH, "cdr2_nt")
  igh_fwr3_nt_col <- paste0(PREFIX_IGH, "fwr3_nt")
  igh_cdr3_nt_col <- paste0(PREFIX_IGH, "cdr3_nt")
  igh_fwr4_nt_col <- paste0(PREFIX_IGH, "fwr4_nt")
  # ... (IGLのnt, IGHのaa, IGLのaaも同様に定義) ...
  igl_fwr1_nt_col <- paste0(PREFIX_IGL, "fwr1_nt")
  igl_cdr1_nt_col <- paste0(PREFIX_IGL, "cdr1_nt")
  # ... (以下同様)
  
  igh_fwr1_col <- paste0(PREFIX_IGH, "fwr1")
  igh_cdr1_col <- paste0(PREFIX_IGH, "cdr1")
  # ... (以下同様)
  
  igl_fwr1_col <- paste0(PREFIX_IGL, "fwr1")
  igl_cdr1_col <- paste0(PREFIX_IGL, "cdr1")
  # ... (以下同様)
  
  df <- df %>%
    dplyr::mutate(
      # IGH ヌクレオチド全長配列
      !!paste0(PREFIX_IGH, "full_length_nt") := paste0(
        dplyr::coalesce(.data[[igh_fwr1_nt_col]], ""), dplyr::coalesce(.data[[igh_cdr1_nt_col]], ""),
        dplyr::coalesce(.data[[igh_fwr2_nt_col]], ""), dplyr::coalesce(.data[[igh_cdr2_nt_col]], ""),
        dplyr::coalesce(.data[[igh_fwr3_nt_col]], ""), dplyr::coalesce(.data[[igh_cdr3_nt_col]], ""),
        dplyr::coalesce(.data[[igh_fwr4_nt_col]], "")
      ),
      # IGL/IGK ヌクレオチド全長配列 (列名の変数を適宜修正)
      !!paste0(PREFIX_IGL, "full_length_nt") := paste0(
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr1_nt")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGL, "cdr1_nt")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr2_nt")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGL, "cdr2_nt")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr3_nt")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGL, "cdr3_nt")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr4_nt")]], "")
      ),
      # IGH アミノ酸全長配列 (列名の変数を適宜修正)
      !!paste0(PREFIX_IGH, "full_length_aa") := paste0(
        dplyr::coalesce(.data[[paste0(PREFIX_IGH, "fwr1")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGH, "cdr1")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGH, "fwr2")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGH, "cdr2")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGH, "fwr3")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGH, "cdr3")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGH, "fwr4")]], "")
      ),
      # IGL/IGK アミノ酸全長配列 (列名の変数を適宜修正)
      !!paste0(PREFIX_IGL, "full_length_aa") := paste0(
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr1")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGL, "cdr1")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr2")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGL, "cdr2")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr3")]], ""), dplyr::coalesce(.data[[paste0(PREFIX_IGL, "cdr3")]], ""),
        dplyr::coalesce(.data[[paste0(PREFIX_IGL, "fwr4")]], "")
      )
    ) %>%
    # もし生成された全長配列が空文字列になった場合（全ての構成要素がNAだった場合）は NA に置換
    dplyr::mutate(
      across(ends_with("full_length_nt") | ends_with("full_length_aa"), ~na_if(., ""))
    )
  
  return(df)
}
