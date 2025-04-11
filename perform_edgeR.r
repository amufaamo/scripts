#!/usr/bin/env Rscript
# ==============================================================================
# R Script: process_featurecounts.R
#
# Description:
#   featureCounts の出力ファイルをマージし、フィルタリング後、
#   CPM および スケール化された logCPM (Z-score) を計算するスクリプト。
#
# Usage:
#   source("process_featurecounts.R")
#   results <- process_featurecounts(
#                 input_dir = "path/to/featurecounts/output",
#                 file_pattern = "\\.txt$", # featureCountsファイルのパターン
#                 sample_metadata_df = your_metadata_dataframe
#              )
#   cpm_data <- results$cpm
#   scaled_logcpm_data <- results$scaled_logcpm
#
# Requirements: edgeR, dplyr, matrixStats
# Last Updated: 2025-04-09
# ==============================================================================

# --- 0. 必要なライブラリの読み込み ---
# stopifnot はエラー時に停止させる堅牢な方法
stopifnot(requireNamespace("edgeR", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE)) # データ操作に使用
stopifnot(requireNamespace("matrixStats", quietly = TRUE)) # rowVarsに使用

library(edgeR)
library(dplyr)
library(matrixStats)

cat("ライブラリ 'edgeR', 'dplyr', 'matrixStats' を読み込みました。\n")

# --- 1. featureCounts ファイルをマージする関数 ---
# 提供された merge_featurecount_data 関数を基にする
merge_featurecounts <- function(input_dir, pattern) {
  # ディレクトリが存在するか確認
  if (!dir.exists(input_dir)) {
    stop("エラー: 指定された入力ディレクトリが存在しません: ", input_dir)
  }

  # ディレクトリ内のファイルリストを取得
  all_files <- list.files(input_dir, pattern = pattern, full.names = TRUE, recursive = FALSE) # recursive = FALSE を明示

  if (length(all_files) == 0) {
     warning(paste("指定されたパターン '", pattern, "' に一致するファイルが",
                   input_dir, "に見つかりませんでした。"))
     return(NULL)
  }

  # ファイルリストから "ensembl" を含むものを除外
  files_to_process <- all_files[!grepl("ensembl", basename(all_files), perl = TRUE)]

  if (length(files_to_process) == 0) {
    warning(paste("パターン '", pattern, "' に一致し、かつファイル名に 'ensembl' を含まないファイルが",
                  input_dir, "に見つかりませんでした。"))
    return(NULL)
  }

  cat("処理対象のファイル (", length(files_to_process), "件):\n", sep="")
  print(basename(files_to_process))

  # 結果を格納するデータフレームを初期化
  merged_df <- NULL

  # フィルタリングされたファイルリストをループ処理
  for (file_path in files_to_process) {
    cat("  読み込み中:", basename(file_path), "...\n")
    # エラーハンドリングを追加
    df <- tryCatch({
        read.table(file_path, header = TRUE, comment.char = "#", sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    }, error = function(e) {
        warning(paste("ファイルの読み込みに失敗しました:", file_path, "-", e$message))
        return(NULL) # エラー時はNULLを返す
    })

    if (is.null(df)) next # 読み込み失敗時はスキップ

    # "Geneid" カラムの存在確認
    if (!"Geneid" %in% colnames(df)) {
      warning(paste("警告: Geneid列が見つかりません。スキップします:", basename(file_path)))
      next
    }

    # カウント列を特定 (featureCounts の標準出力は通常7列目だが、より堅牢にする)
    # 最後の列をカウントデータと仮定する (より安全な方法は列名で指定すること)
    # ここでは提供されたスクリプトに合わせて7列目を仮定するが、注意が必要
    count_col_index <- 7
    if (ncol(df) < count_col_index) {
       warning(paste("警告: 列数が", count_col_index, "未満です。スキップします:", basename(file_path)))
       next
    }
    count_col_name <- colnames(df)[count_col_index]
    cat("    カウント列として使用:", count_col_name, "\n")

    # サンプル名をファイル名から取得 (より堅牢な方法を推奨)
    # ここでは拡張子を除去するだけ
    sample_name <- gsub(pattern, "", basename(file_path))
    # 必要に応じて接頭辞なども除去: sample_name <- gsub("^prefix_", "", sample_name)
    cat("    サンプル名として使用:", sample_name, "\n")


    # 必要な列（Geneidとカウントデータ）を抽出
    # check.names=FALSE を read.table で指定したので、特殊文字が含まれていても大丈夫なはず
    count_data <- df[, c("Geneid", count_col_name), drop = FALSE]
    # カウント列が数値でない場合は警告し、数値に変換を試みる
    if (!is.numeric(count_data[[2]])) {
        warning(paste("警告:", basename(file_path), "のカウント列", count_col_name, "が数値ではありません。数値への変換を試みます。"))
        count_data[[2]] <- suppressWarnings(as.numeric(count_data[[2]]))
    }
    colnames(count_data) <- c("Geneid", sample_name) # カラム名を変更

    # 結合
    if (is.null(merged_df)) {
      merged_df <- count_data
    } else {
      merged_df <- merge(merged_df, count_data, by = "Geneid", all = TRUE) # all=TRUE で outer join
    }
  } # ループ終了

  # マージ後に発生したNAを0で置換（Outer join で片方にしか存在しない遺伝子のため）
  if (!is.null(merged_df)) {
      # Geneid列以外を対象にNAを0に置換
      numeric_cols <- setdiff(colnames(merged_df), "Geneid")
      for(col in numeric_cols) {
          merged_df[[col]][is.na(merged_df[[col]])] <- 0
      }
      cat("データのマージが完了しました。最終的な次元:", paste(dim(merged_df), collapse=" x "), "\n")
  } else {
      cat("有効なデータがマージされませんでした。\n")
  }

  return(merged_df)
}


# --- 2. CPM と Scaled logCPM を計算するメイン関数 ---
process_featurecounts <- function(
    input_dir,                    # featureCounts出力ファイルがあるディレクトリ
    file_pattern,                 # featureCounts出力ファイルのパターン (例: "\\.txt$")
    sample_metadata_df,           # サンプルメタデータ (必須: サンプル名列, group列)
    sample_name_col = "SampleName", # メタデータ内のサンプル名列名
    group_col = "group",            # メタデータ内のグループ列名
    prior_count_logcpm = 2        # logCPM計算時のprior count
    ) {

  cat("\n--- データ処理プロセス開始 ---\n")

  # --- 2a. データのマージ ---
  cat("Step 1: featureCounts ファイルのマージを開始します...\n")
  merged_counts_df <- merge_featurecounts(input_dir = input_dir, pattern = file_pattern)

  if (is.null(merged_counts_df) || nrow(merged_counts_df) == 0 || ncol(merged_counts_df) <= 1) {
    stop("エラー: featureCounts ファイルのマージに失敗したか、有効なデータが含まれていません。")
  }

  # --- 2b. メタデータとカウントデータの整合性チェックと準備 ---
  cat("Step 2: メタデータとカウントデータの整合性を確認し、準備します...\n")

  # メタデータに必要な列が存在するか確認
  required_meta_cols <- c(sample_name_col, group_col)
  if (!all(required_meta_cols %in% colnames(sample_metadata_df))) {
    stop("エラー: メタデータに必須の列 (", paste(required_meta_cols, collapse=", "), ") が存在しません。")
  }

  # マージされたカウントデータのサンプル名を取得 (Geneid列を除く)
  count_sample_names <- colnames(merged_counts_df)[-1]

  # メタデータからカウントデータに存在するサンプルを抽出
  metadata_subset <- sample_metadata_df[sample_metadata_df[[sample_name_col]] %in% count_sample_names, ]

  if (nrow(metadata_subset) == 0) {
      stop("エラー: メタデータに、マージされたカウントデータのサンプル名が見つかりません。")
  }
  if (nrow(metadata_subset) != length(count_sample_names)) {
      warning("警告: メタデータ内のサンプル数とカウントデータのサンプル数が一致しません。カウントデータに存在するサンプルのみを使用します。")
      # カウントデータ側もメタデータに存在するサンプルに絞る
      samples_to_keep <- metadata_subset[[sample_name_col]]
      merged_counts_df <- merged_counts_df[, c("Geneid", samples_to_keep), drop = FALSE]
      count_sample_names <- samples_to_keep # 更新
  }

  # カウントデータの列順序に合わせてメタデータを並び替え
  metadata_ordered <- metadata_subset[match(count_sample_names, metadata_subset[[sample_name_col]]), ]

  # グループ情報をファクターとして取得
  analysis_group <- factor(metadata_ordered[[group_col]])
  cat("使用するグループ情報:\n")
  print(table(analysis_group))

  # カウントデータを数値行列に変換 (Geneidをrownamesに設定)
  count_matrix <- as.matrix(merged_counts_df[, -1, drop = FALSE])
  rownames(count_matrix) <- merged_counts_df$Geneid
  storage.mode(count_matrix) <- "numeric" # Ensure numeric storage

  # --- 2c. edgeR オブジェクト作成とフィルタリング ---
  cat("Step 3: DGEListオブジェクトを作成し、低発現遺伝子をフィルタリングします...\n")
  y <- DGEList(counts = count_matrix, group = analysis_group)

  # フィルタリングのためのデザイン行列を作成 (グループ情報のみを使用)
  design <- model.matrix(~ analysis_group)
  cat("フィルタリングに使用するデザイン行列 (最初の数行):\n")
  print(head(design))

  # filterByExpr を使用してフィルタリング対象を決定
  keep <- filterByExpr(y, design = design)
  y_filtered <- y[keep, , keep.lib.sizes = FALSE]

  n_genes_before <- nrow(y)
  n_genes_after <- nrow(y_filtered)
  cat("フィルタリング完了。\n")
  cat("  フィルタリング前の遺伝子数:", n_genes_before, "\n")
  cat("  フィルタリング後の遺伝子数:", n_genes_after, "\n")

  if (n_genes_after == 0) {
      stop("エラー: フィルタリング後に遺伝子が残りませんでした。フィルタリング条件を確認してください。")
  }

  # --- 2d. 正規化係数の計算 ---
  cat("Step 4: TMM正規化係数を計算します...\n")
  y_filtered <- calcNormFactors(y_filtered)
  cat("正規化係数の計算完了。\n")

  # --- 2e. CPMの計算 ---
  cat("Step 5: CPM (Counts Per Million) を計算します...\n")
  cpm_matrix <- edgeR::cpm(y_filtered, log = FALSE)
  # データフレームに変換し、Geneid列を追加
  cpm_df <- as.data.frame(round(cpm_matrix, 3)) %>%
              tibble::rownames_to_column("Geneid")
  cat("CPM計算完了。\n")

  # --- 2f. Scaled logCPM (Z-score) の計算 ---
  cat("Step 6: スケール化された logCPM (Z-score) を計算します...\n")
  logcpm_matrix <- edgeR::cpm(y_filtered, log = TRUE, prior.count = prior_count_logcpm)

  # 遺伝子ごとの分散を計算 (ゼロ分散のチェック)
  gene_vars <- matrixStats::rowVars(logcpm_matrix, na.rm = TRUE)
  zero_var_genes <- sum(gene_vars < 1e-8, na.rm = TRUE) # 非常に小さい分散もゼロとみなす
  if (zero_var_genes > 0) {
    showNotification(paste("スケール化警告: 分散がほぼゼロの遺伝子が", zero_var_genes,
                           "個あります。これらの遺伝子のZ-scoreは NaN になります。"),
                     type="warning", duration=10)
  }

  # scale関数で行ごとにZ-score化 (行列を転置して適用し、再度転置して元に戻す)
  # na.rm=TRUE で欠損値があっても計算できるようにする (通常logCPMではないはずだが念のため)
  scaled_logcpm_matrix <- t(scale(t(logcpm_matrix), center = TRUE, scale = TRUE))
  scaled_logcpm_matrix[is.nan(scaled_logcpm_matrix)] <- 0 # NaNが発生した場合0に置換（ゼロ分散遺伝子など）

  # データフレームに変換し、Geneid列を追加
  scaled_logcpm_df <- as.data.frame(round(scaled_logcpm_matrix, 3)) %>%
                        tibble::rownames_to_column("Geneid")
  cat("Scaled logCPM計算完了。\n")

  # --- 2g. 結果の返却 ---
  cat("--- データ処理プロセス終了 ---\n")
  results <- list(
    cpm = cpm_df,
    scaled_logcpm = scaled_logcpm_df,
    n_genes_filtered = n_genes_after,
    n_genes_initial = n_genes_before
  )
  return(results)
}

# --- スクリプト読み込み完了メッセージ ---
cat("関数 'merge_featurecounts' と 'process_featurecounts' が定義されました。\n")
cat("使用法: source(\"process_featurecounts.R\") してから関数を呼び出してください。\n")
cat("例:\n")
cat("  # sample_meta <- data.frame(SampleName=paste0('Sample', 1:4), group=rep(c('A','B'), each=2))\n")
cat("  # results <- process_featurecounts(input_dir='path/to/counts', file_pattern='\\.txt$', sample_metadata_df=sample_meta)\n")
cat("  # cpm_table <- results$cpm\n")
cat("  # scaled_logcpm_table <- results$scaled_logcpm\n")


# --- (オプション) スクリプトとして直接実行された場合のサンプル実行 ---
# この部分は source() で読み込んだ場合は実行されません
if (!interactive() && sys.nframe() == 0) {
  cat("\n--- スクリプト直接実行時のサンプル解析開始 ---\n")

  # サンプルデータとディレクトリの準備
  set.seed(123)
  temp_dir <- tempdir()
  sample_output_dir <- file.path(temp_dir, "sample_featurecounts_output")
  if (!dir.exists(sample_output_dir)) dir.create(sample_output_dir)
  cat("サンプルデータ出力先:", sample_output_dir, "\n")

  # サンプルメタデータ作成
  sample_names <- paste0("Sample", 1:6)
  sample_metadata <- data.frame(
    SampleName = sample_names,
    group = rep(c("Control", "Treatment"), each = 3),
    batch = rep(c("B1", "B2", "B3"), times = 2),
    stringsAsFactors = FALSE
  )
  cat("サンプルメタデータ:\n")
  print(sample_metadata)

  # ダミーのfeatureCountsファイルを作成
  gene_ids <- paste0("Gene", 1:500)
  for (s_name in sample_names) {
    # ランダムなカウントデータ生成 (一部ゼロや低カウントを含むように)
    counts <- rnbinom(500, mu = sample(c(0, 5, 50, 150), 500, replace = TRUE, prob=c(0.1, 0.3, 0.4, 0.2)), size = 10)
    # featureCountsのような形式のデータフレームを作成
    dummy_fc_df <- data.frame(
      Geneid = gene_ids,
      Chr = "chr1", Start = 1, End = 1000, Strand = "+", Length = 1000,
      CountCol = counts, # 7列目にカウントデータ
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(dummy_fc_df)[7] <- s_name # 7列目の名前をサンプル名に (注意: 実際のfeatureCounts出力に合わせる必要あり)
    # ファイルに書き出し (タブ区切り、ヘッダー付き、コメント行付き)
    file_path <- file.path(sample_output_dir, paste0(s_name, ".featureCounts.txt"))
    cat("# Program: featureCounts vX.Y.Z\n", file = file_path)
    cat("# Command: featureCounts ...\n", file = file_path, append = TRUE)
    write.table(dummy_fc_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
  }
  # "ensembl"を含むダミーファイルも作成（除外されることを確認するため）
   dummy_fc_df_ensembl <- data.frame( Geneid = gene_ids[1:10], CountCol = 1:10 ); colnames(dummy_fc_df_ensembl)[2] <- "Sample_ensembl"
   write.table(dummy_fc_df_ensembl, file = file.path(sample_output_dir, "Sample_with_ensembl_name.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


  cat("\n--- サンプルデータで process_featurecounts を実行 ---\n")
  sample_results <- NULL
  tryCatch({
    sample_results <- process_featurecounts(
      input_dir = sample_output_dir,
      file_pattern = "\\.featureCounts\\.txt$", # 作成したファイルのパターン
      sample_metadata_df = sample_metadata,
      sample_name_col = "SampleName", # メタデータの列名
      group_col = "group"             # メタデータの列名
    )
  }, error = function(e) {
    cat("\nサンプル解析中にエラーが発生しました:\n", e$message, "\n")
  })

  # 結果の表示
  if (!is.null(sample_results)) {
    cat("\n--- サンプル結果の概要 ---\n")
    cat("フィルタリング前後の遺伝子数:", sample_results$n_genes_initial, "->", sample_results$n_genes_filtered, "\n")

    cat("\nCPM データフレーム (最初の数行):\n")
    print(head(sample_results$cpm))

    cat("\nScaled logCPM データフレーム (最初の数行):\n")
    print(head(sample_results$scaled_logcpm))
  } else {
    cat("\nサンプル解析で結果が得られませんでした。\n")
  }

  # クリーンアップ（不要ならコメントアウト）
  # unlink(sample_output_dir, recursive = TRUE)
  # cat("\nサンプルディレクトリを削除しました:", sample_output_dir, "\n")

  cat("\n--- サンプル解析終了 ---\n")
}
