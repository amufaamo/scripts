merge_featurecount_data <- function(file_fullpath_list) {
  # 結果を格納するデータフレームを初期化
  merged_df <- NULL
  
  # ファイルリストをループ処理
  for (file_path in file_fullpath_list) {
    # ファイルを読み込み (header, comment.char, sep は featurecount の出力形式に合わせて調整)
    df <- read.table(file_path, header = TRUE, comment.char = "#", sep = "\t")
    
    # サンプル名をファイル名から取得
    sample_name <- gsub(".txt$", "", basename(file_path))  # 拡張子を削除
    sample_name <- gsub("^241225_", "", sample_name)  # 接頭辞を削除

    # "Geneid" カラムの存在確認
    if (!"Geneid" %in% colnames(df)) {
      stop(paste("Error: Geneid列が見つかりません:", file_path))
    }
    
    # 7列目のカラム名を取得（異なる場合に備えて自動判定）
    count_col <- colnames(df)[7]

    # 必要な列（Geneidとカウントデータ）を抽出
    count_data <- df[, c("Geneid", count_col)]
    colnames(count_data) <- c("Geneid", sample_name)  # カラム名を変更

    # 最初のファイルの場合
    if (is.null(merged_df)) {
      merged_df <- count_data
    } else {
      # Geneid をキーとして結合
      merged_df <- merge(merged_df, count_data, by = "Geneid", all = TRUE)
    }
  }

  return(merged_df)
  
  # # 結果をファイルに出力
  # write.table(merged_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # cat("結合されたデータフレームを", output_file, "に出力しました。\n")
}
