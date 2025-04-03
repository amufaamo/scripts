# --- 3. edgeR解析を実行する関数 (1比較ペア版、QC/Norm/Dispのサンプル選択オプション付き) ---

perform_edgeR_one_comparison <- function(
    countdata,
    group,
    comparison_pair,
    output_dir = NULL,
    use_all_samples_for_qc_norm_disp = TRUE # ★新しい引数 (デフォルトはTRUE)
    ) {

  # --- 入力データの検証 ---
  # ... (前回の検証コードはそのまま) ...
   if (!is.data.frame(countdata) && !is.matrix(countdata)) { stop("エラー: countdata は data.frame または matrix である必要があります。") }
   if (!is.vector(group) && !is.factor(group)) { stop("エラー: group は vector または factor である必要があります。") }
   if (!is.character(comparison_pair) || length(comparison_pair) != 2) { stop("エラー: comparison_pair は c('target', 'reference') 形式の長さ2の文字ベクトルである必要があります。") }
   if (ncol(countdata) != length(group)) { stop("エラー: カウントデータの列数とグループ情報のサンプル数が一致しません。") }

  original_group_factor <- factor(group) # 元のグループ情報 (因子)
  target_group <- comparison_pair[1]
  reference_group <- comparison_pair[2]

  if (!target_group %in% levels(original_group_factor) || !reference_group %in% levels(original_group_factor)) { stop(paste("エラー: comparison_pair", paste(comparison_pair, collapse=" vs "), "に含まれるグループが group 情報に存在しません。")) }
  if (target_group == reference_group) { stop("エラー: 比較ペアの target と reference が同じグループです。") }
  comparison_name <- paste0(target_group, "_vs_", reference_group)
  cat("入力データの検証OK\n")
  cat("比較対象:", comparison_name, "\n")


  # --- ★使用するサンプルを決定 ---
  if (use_all_samples_for_qc_norm_disp) {
      cat("初期ステップ (QC, Norm, Disp) で「全サンプル」を使用します。\n")
      analysis_countdata <- countdata
      analysis_group <- original_group_factor
  } else {
      cat("初期ステップ (QC, Norm, Disp) で「比較ペアのサンプルのみ」を使用します。\n")
      samples_to_keep_logical <- original_group_factor %in% comparison_pair
      if (sum(samples_to_keep_logical) < 2) { # 最低2サンプル必要
          stop("エラー: 比較ペアに属するサンプル数が2未満のため、ペアのみでの解析は実行できません。")
      }
      analysis_countdata <- countdata[, samples_to_keep_logical, drop = FALSE]
      # 使用しないレベルを削除
      analysis_group <- droplevels(original_group_factor[samples_to_keep_logical])

      # サブセット後に両方のグループが存在するか確認
      if (nlevels(analysis_group) < 2) {
           stop("エラー: サブセット後に比較に必要な2つのグループが存在しません。各グループにサンプルがあるか確認してください。")
      }
      cat("サブセット後のサンプル数:", ncol(analysis_countdata), "\n")
      print("サブセット後のグループ情報:")
      print(table(analysis_group))
  }
  cat("解析に使用するグループ情報:\n"); print(table(analysis_group))


  # --- 出力ディレクトリとファイル保存設定 ---
  save_files <- !is.null(output_dir)
  # ... (保存設定のコードはそのまま) ...
  if (save_files) {
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE); cat("出力ディレクトリを作成しました:", output_dir, "\n") }
     cat("結果ファイルは次のディレクトリに保存されます:", output_dir, "\n")
  } else {
    cat("ファイル保存はスキップされます。\n")
  }


  # --- edgeR 解析ステップ (analysis_countdata と analysis_group を使用) ---

  # 1. DGEListオブジェクトの作成 (★analysis_***を使用)
  y <- DGEList(counts = analysis_countdata, group = analysis_group)
  cat("DGEListオブジェクトを作成しました。\n")

  # 2. フィルタリング (★analysis_groupを使用)
  keep <- filterByExpr(y, group = analysis_group) # デザイン行列を使う場合は design= を指定
   if(sum(keep) == 0) { stop("エラー: フィルタリング後に遺伝子が残りませんでした。filterByExprの基準を確認してください。") }
  y <- y[keep, , keep.lib.sizes = FALSE]
  cat("低発現遺伝子をフィルタリングしました。残りの遺伝子数:", nrow(y), "\n")

  # 3. 正規化 (TMM法)
  y <- calcNormFactors(y)
  cat("TMM正規化を実行しました。\n")

  # 4. 分散推定 (★analysis_groupを使用)
  design <- model.matrix(~ 0 + analysis_group)
  colnames(design) <- levels(analysis_group) # designの列名は analysis_group のレベルに依存
  print("デザイン行列:")
  print(design)
   y <- tryCatch({ estimateDisp(y, design) }, error = function(e) { stop(paste("分散推定中にエラー:", e$message)) })
  cat("分散推定を実行しました。\n")

  # 5. BCVプロット
  print("BCVプロットを表示します...")
  plotBCV(y)
  title("Biological Coefficient of Variation (BCV) Plot")

  # 6. MDSプロット
  print("MDSプロットを表示します...")
  # ... (MDSプロットのコードはそのまま or ggplot版) ...
   mds_data <- plotMDS(y, plot=FALSE)
   plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=16, main="MDS Plot", xlab="Dim 1", ylab="Dim 2", las=1)
   text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7)
   legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8)


  # 7. 差次発現解析 (GLM QL F-test) - design と比較ペアを使用
  cat("\n--- 差次発現解析を開始:", comparison_name, "---\n")
  fit_glm <- glmQLFit(y, design)

  # コントラスト作成 (target/reference は元の comparison_pair の名前を使用)
  # makeContrasts は design 行列の列名 (levels=design) を参照する
  contrast_formula <- paste0(target_group, " - ", reference_group)
  con_target_vs_ref <- makeContrasts(contrasts = contrast_formula, levels = design)
  print("コントラスト行列:")
  print(con_target_vs_ref)

  # GLM QL F-test 実行
  qlf_result_single <- glmQLFTest(fit_glm, contrast = con_target_vs_ref)

  # 8. 結果要約とMDプロット
  cat("結果の要約 (FDR < 0.05):\n")
  print(summary(decideTests(qlf_result_single)))
  print(paste("MDプロットを表示します (", comparison_name, ")..."))
  plotMD(qlf_result_single, main = comparison_name)
  abline(h = c(-1, 1), col = "blue", lty = 2)

  # 9. 結果をRDSファイルとして保存 (Optional)
  # ... (保存コードはそのまま) ...
   if (save_files) {
     date_prefix <- format(Sys.Date(), "%y%m%d")
     filename <- file.path(output_dir, paste0(date_prefix, '_qlf_', comparison_name, '.rds'))
     tryCatch({ saveRDS(qlf_result_single, file = filename); cat("結果をRDSとして保存しました:", filename, "\n") }, error = function(e) { warning("RDSファイルの保存中にエラー:", e$message) })
   }


  cat("\n--- 解析完了:", comparison_name, "---\n")

  # 単一の比較結果 (qlfオブジェクト) を返す
return(topTags(qlf_result_single, n = Inf)$table)
#  return(qlf_result_single$table)
}

cat("単一比較用解析関数 'perform_edgeR_one_comparison' を定義しました (サンプル選択オプション付き)。\n")
