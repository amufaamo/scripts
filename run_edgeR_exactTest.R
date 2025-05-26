run_edgeR_exactTest <- function(counts_df, group_vector, comparison_pair) {

  # 引数のチェック (簡易版)
  if (!is.data.frame(counts_df) && !is.matrix(counts_df)) { #行列も許容するように修正
    stop("counts_df はデータフレームまたは行列である必要があります。")
  }
  if (ncol(counts_df) != length(group_vector)) {
    stop("counts_df の列数と group_vector の長さが一致しません。")
  }
  if (length(comparison_pair) != 2) {
    stop("comparison_pair は2つのグループ名を含むベクトルである必要があります。")
  }
  if (!all(comparison_pair %in% unique(group_vector))) {
    stop("comparison_pair で指定されたグループ名が group_vector に存在しません。")
  }
  if (comparison_pair[1] == comparison_pair[2]) {
    stop("comparison_pair で指定された2つのグループ名は異なるものである必要があります。")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("パッケージ 'edgeR' がインストールされていません。インストールしてください。")
  }

  cat("比較対象:", comparison_pair[1], "vs", comparison_pair[2], "\n")

  # --- ここから修正 ---
  # 指定された2つのグループに属するサンプルのみを抽出
  group1_name <- comparison_pair[1]
  group2_name <- comparison_pair[2]

  samples_to_keep_logical <- group_vector == group1_name | group_vector == group2_name
  
  if (sum(samples_to_keep_logical) < 2) {
      stop(paste("指定された比較ペア (", group1_name, ", ", group2_name, ") に該当するサンプルが2つ未満です。解析を続行できません。", sep=""))
  }

  counts_df_subset <- counts_df[, samples_to_keep_logical, drop = FALSE]
  group_vector_subset <- group_vector[samples_to_keep_logical]
  
  cat("サブセット化後:", ncol(counts_df_subset), "サンプルが選択されました。\n")
  cat("選択されたグループ:", paste(unique(group_vector_subset), collapse=", "), "\n")

  # サブセット化後、各グループに最低1サンプルが存在することを確認
  table_subset_groups <- table(group_vector_subset)
  if (!all(comparison_pair %in% names(table_subset_groups)) || any(table_subset_groups[comparison_pair] == 0)) {
      stop(paste("サブセット化の結果、比較ペア (", group1_name, ", ", group2_name, ") の一方または両方のグループのサンプル数が0になりました。", sep=""))
  }
  # --- ここまで修正 ---

  # 1. カウントデータの行列への変換
  # 関数内では、入力データフレームが行名に遺伝子IDを持つ数値データであると仮定
  count_data_matrix <- as.matrix(counts_df_subset) # 修正: サブセット化されたデータを使用

  # 2. グループ情報の定義
  # group_vector を factor に変換し、比較の順番を comparison_pair に合わせる
  group <- factor(group_vector_subset, levels = comparison_pair) # 修正: サブセット化されたグループ情報を使用

  # 3. DGEListオブジェクトの作成
  d <- edgeR::DGEList(counts = count_data_matrix, group = group)
  cat("DGEListオブジェクト作成完了。\n")
  cat("サンプル数:", ncol(d), "初期遺伝子数:", nrow(d), "\n")

  # 4. 低発現遺伝子のフィルタリング
  # filterByExpr は group 情報を使ってフィルタリング閾値を自動で決定
  keep <- edgeR::filterByExpr(d, group = group) # group引数はDGEListオブジェクトから取得されるので、ここではdだけで良い場合もある
  d <- d[keep, , keep.lib.sizes = FALSE]
  cat("フィルタリング後の遺伝子数:", nrow(d), "\n")

  if (nrow(d) == 0) {
    stop("フィルタリングの結果、解析対象の遺伝子が0になりました。フィルタリング基準を確認してください。")
  }

  # 5. 正規化 (TMM)
  d <- edgeR::calcNormFactors(d)
  cat("TMM正規化完了。\n")

  # 6. 分散の推定
  # コモン分散の推定
  d <- edgeR::estimateCommonDisp(d)
  cat("コモン分散推定完了。\n")
  # トレンド分散の推定 (オプションだが一般的に推奨)
  d <- edgeR::estimateTrendedDisp(d)
  cat("トレンド分散推定完了。\n")
  # タグワイズ分散の推定
  d <- edgeR::estimateTagwiseDisp(d)
  cat("タグワイズ分散推定完了。\n")

  # 7. 差次発現検定 (exactTest)
  # comparison_pair の順番で比較 (例: comparison_pair = c("control", "treatment") なら treatment vs control)
  # DGEList 作成時に group の factor の levels を comparison_pair で指定しておけば、
  # 自動的に意図した比較になる。
  result_test <- edgeR::exactTest(d) # pair は DGEList の group レベルに基づき自動設定
  cat("exactTest実行完了。\n")

  # 8. 結果の取得とデータフレームへの変換
  # フィルタリング後の全遺伝子の結果を取得
  result_table <- edgeR::topTags(result_test, n = nrow(d$counts))
  result_df <- as.data.frame(result_table)
  cat("結果テーブル作成完了。\n")

  return(result_df)
}
