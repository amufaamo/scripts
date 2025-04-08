#!/usr/bin/env Rscript

# ==============================================================================
# R Script: perform_edgeR.R
#
# Description:
#   edgeRパッケージを用いて、RNA-Seqカウントデータから差次発現遺伝子解析を
#   行うための関数 `perform_edgeR` を定義します。
#   この関数は、単一の比較ペア、または複数の比較ペアのリストを入力として受け付け、
#   結果として `topTags` の結果テーブル (data.frame) またはそのリストを返します。
#   バッチ効果の補正や、QC/正規化/分散推定ステップで使用するサンプルの選択が可能です。
#   オプションで各種プロットやバッチ効果補正後のPCAプロットを生成・保存します。
#   有意性の判定に使用するFDRとlogFCの閾値を引数で指定できます。
#
# Usage:
#   source("perform_edgeR.R")
#   # 単一比較の場合: c("Case", "Control") を渡す -> 結果テーブル(df)が返る
#   result_table_single <- perform_edgeR(..., comparisons = c("Case", "Control"), ...)
#   # 複数比較の場合: list(c("A", "B"), c("A", "C")) を渡す -> 結果テーブルのリストが返る
#   results_table_list <- perform_edgeR(..., comparisons = list(A_vs_B=c("A","B"), A_vs_C=c("A","C")), ...)
#
# Last Updated: 2025-04-08 (最終確認版)
# Requirements: edgeR, limma, matrixStats
# ==============================================================================

# --- 0. 必要なライブラリの読み込み ---
stopifnot(requireNamespace("edgeR", quietly = TRUE))
stopifnot(requireNamespace("limma", quietly = TRUE))
stopifnot(requireNamespace("matrixStats", quietly = TRUE))
library(edgeR)
library(limma)
library(matrixStats)

cat("ライブラリ 'edgeR', 'limma', 'matrixStats' を読み込みました。\n")


# --- 3. edgeR解析を実行する関数 (最終版) ---

perform_edgeR <- function(
    countdata,
    group,
    comparisons, # 単一ペアのベクトル または ペアのリスト
    batch_info = NULL,
    output_dir = NULL,
    use_all_samples_for_qc_norm_disp = TRUE,
    fdr_threshold = 0.05, # FDR閾値
    logfc_threshold = 0    # logFC閾値 (絶対値)
    ) {

  # --- 入力データの検証 ---
  cat("--- 入力データの検証を開始 ---\n")
  if (!is.data.frame(countdata) && !is.matrix(countdata)) { stop("エラー: countdata は data.frame または matrix である必要があります。") }
  if (!is.vector(group) && !is.factor(group)) { stop("エラー: group は vector または factor である必要があります。") }
  if (ncol(countdata) != length(group)) { stop("エラー: カウントデータの列数とグループ情報のサンプル数が一致しません。") }
  original_group_factor <- factor(group)
  if (!is.numeric(fdr_threshold) || fdr_threshold < 0 || fdr_threshold > 1) { stop("エラー: fdr_threshold は 0 から 1 の数値である必要があります。") }
  if (!is.numeric(logfc_threshold) || logfc_threshold < 0) { stop("エラー: logfc_threshold は 0 以上の数値である必要があります。") }
  cat(paste0("使用するFDR閾値: ", fdr_threshold, "\n"))
  cat(paste0("使用するlogFC閾値 (絶対値): ", logfc_threshold, "\n"))

  # comparisons 引数のチェックと処理
  is_list_input <- FALSE; internal_comparison_list <- list()
  if (is.list(comparisons) && !is.data.frame(comparisons)) { if(length(comparisons) == 0) stop("エラー: comparisons リストが空です。"); is_list_input <- TRUE; internal_comparison_list <- comparisons; cat("複数の比較ペアがリストとして指定されました。\n") }
  else if (is.character(comparisons) && length(comparisons) == 2) { is_list_input <- FALSE; internal_comparison_list <- list(comparisons); cat("単一の比較ペアがベクトルとして指定されました。\n") }
  else { stop("エラー: 'comparisons' 引数は、c('target', 'reference')形式のベクトル、またはそのリストである必要があります。") }

  valid_pairs_exist <- FALSE; processed_list <- list(); list_names <- names(internal_comparison_list); valid_names <- c()
  for(i in 1:length(internal_comparison_list)) {
      pair <- internal_comparison_list[[i]]
      # リスト要素に名前があればそれを使う、なければペアから生成
      pair_name_raw <- if(!is.null(list_names) && nzchar(list_names[i])) list_names[i] else paste(pair[1], pair[2], sep="_vs_")

      if(is.character(pair) && length(pair) == 2 && pair[1] != pair[2] && all(pair %in% levels(original_group_factor))) {
          valid_pairs_exist <- TRUE
          processed_list[[length(processed_list) + 1]] <- pair # 有効なペアのみ追加
          valid_names <- c(valid_names, pair_name_raw) # 対応する名前を保持
          if(is_list_input) cat("  - OK:", paste(pair, collapse=" vs "), "( Name:", pair_name_raw, ")\n")
      } else {
          warning("リスト内の無効なペア:", paste(pair, collapse=" / "), " - スキップされます。グループ名がgroup引数に存在するか、ペアが同一でないか確認してください。")
      }
  }
   internal_comparison_list <- processed_list # 有効なペアのみのリストで上書き
   names(internal_comparison_list) <- valid_names # リストに名前を再設定

   if (length(internal_comparison_list) == 0) {
       stop("エラー: 有効な比較ペアが指定されませんでした。")
   }

  # バッチ情報の検証
  analysis_batch <- NULL; original_batch_factor <- NULL
  if (!is.null(batch_info)) {
      if (length(batch_info) != ncol(countdata)) { stop("エラー: バッチ情報の長さがカウントデータの列数と一致しません。") }
      original_batch_factor <- factor(batch_info)
      cat("バッチ情報が提供されました。\n")
  } else {
      cat("バッチ情報は提供されていません。\n")
  }

  # グループとバッチのクロス集計を表示
  if (!is.null(original_batch_factor)) {
      cat("\n--- グループとバッチのクロス集計 ---\n")
      if(length(original_group_factor) == length(original_batch_factor)) {
           cross_tab <- table(Group = original_group_factor, Batch = original_batch_factor)
           print(cross_tab)
           if(any(cross_tab == 0)){
                warning("警告: グループとバッチの組み合わせにサンプル数が0のセルがあります。交絡の可能性があります。")
           }
      } else {
          warning("警告: クロス集計スキップ。グループとバッチの長さが一致しません。")
      }
      cat("---------------------------------\n\n")
  }

  # --- 使用するサンプルを決定 ---
  if (use_all_samples_for_qc_norm_disp) {
      cat("初期ステップで「全サンプル」を使用します。\n")
      analysis_countdata <- countdata
      analysis_group <- original_group_factor
      if (!is.null(original_batch_factor)) { analysis_batch <- original_batch_factor }
  } else {
      cat("初期ステップで「比較ペアのサンプルのみ」を使用します (注意: 複数比較リストの場合、最初の有効なペアに基づきます)。\n")
      # 最初の有効なペアを取得
      current_pair_for_subset <- internal_comparison_list[[1]]
      samples_to_keep_logical <- original_group_factor %in% current_pair_for_subset
      if (sum(samples_to_keep_logical) < 2) { stop("エラー: 比較ペアに属するサンプル数が2未満のため、ペアのみでの解析は実行できません。") }
      analysis_countdata <- countdata[, samples_to_keep_logical, drop = FALSE]
      analysis_group <- droplevels(original_group_factor[samples_to_keep_logical])
      if (nlevels(analysis_group) < 2) { stop("エラー: サブセット後に比較に必要な2つのグループが存在しません。") }
      if (!is.null(original_batch_factor)) {
          analysis_batch <- droplevels(original_batch_factor[samples_to_keep_logical])
          if(nlevels(analysis_batch) <= 1) { analysis_batch <- NULL } # 1レベル以下ならモデルに含めない
      }
      cat("サブセット後のサンプル数:", ncol(analysis_countdata), "\n")
      cat("サブセット後のグループ情報:\n"); print(table(analysis_group))
      if (!is.null(analysis_batch)) { cat("サブセット後のバッチ情報:\n"); print(table(analysis_batch)) } else { cat("サブセット後のバッチ情報: なし\n")}
  }
  cat("解析に使用するグループ情報:\n"); print(table(analysis_group))
  if (!is.null(analysis_batch)) { cat("解析に使用するバッチ情報:\n"); print(table(analysis_batch)) }

  # --- 出力ディレクトリとファイル保存設定 ---
  save_files <- !is.null(output_dir)
  if (save_files) {
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE); cat("出力ディレクトリを作成しました:", output_dir, "\n") }
     cat("結果ファイルは次のディレクトリに保存されます:", output_dir, "\n")
  } else {
    cat("ファイル保存はスキップされます。\n")
  }


  # --- edgeR 解析ステップ (共通部分) ---
  cat("\n--- edgeR 解析ステップ開始 (共通部分) ---\n")
  # 1. DGEList
  y <- DGEList(counts = analysis_countdata, group = analysis_group)
  if (!is.null(analysis_batch)) { y$samples$batch <- analysis_batch }
  cat("1. DGEListオブジェクト作成完了。\n")

  # 2. Design Matrix
  group_term_name <- "analysis_group"; batch_term_name <- "analysis_batch"
  if (!is.null(analysis_batch)) {
       design_formula <- as.formula(paste("~ 0 +", group_term_name, "+", batch_term_name))
       # model.matrixに渡すデータフレームを明示的に作成
       model_data <- data.frame(analysis_group=analysis_group)
       model_data[[batch_term_name]] <- analysis_batch # 列名を変数で指定
       design <- model.matrix(design_formula, data=model_data)
       cat("2. デザイン行列を 'group + batch' で作成。\n")
  } else {
       design_formula <- as.formula(paste("~ 0 +", group_term_name))
       model_data <- data.frame(analysis_group=analysis_group)
       design <- model.matrix(design_formula, data=model_data)
       # インターセプトなし、グループのみの場合、列名は自動でレベル名になることが多いが、念のため
       colnames(design) <- levels(analysis_group)
       cat("2. デザイン行列を 'group' のみで作成。\n")
  }
  # Rで有効な列名に変換
  colnames(design) <- make.names(colnames(design))
  print("デザイン行列:")
  print(design)

  # Rank check
  design_rank <- qr(design)$rank; num_columns <- ncol(design)
  if (design_rank < num_columns) { warning(paste("警告: デザイン行列がフルランクではありません (ランク =", design_rank, ", 列数 =", num_columns, "). 交絡の可能性があります。")) }
  else { cat("   デザイン行列はフルランクです。\n") }

  # 3. Filtering
  cat("3. フィルタリング中...\n")
  keep <- filterByExpr(y, design = design)
  if(sum(keep) == 0) { stop("エラー: フィルタリング後に遺伝子が残りませんでした。") }
  y <- y[keep, , keep.lib.sizes = FALSE]
  cat("   完了。残りの遺伝子数:", nrow(y), "\n")

  # 4. Normalization
  cat("4. TMM正規化中...\n"); y <- calcNormFactors(y); cat("   完了。\n")

  # 5. Dispersion Estimation
  cat("5. 分散推定中...\n"); y <- tryCatch({ estimateDisp(y, design) }, error = function(e) { stop(paste("分散推定中にエラー:", e$message)) }); cat("   完了。\n")

  # 6. BCV Plot
  cat("6. BCVプロット生成中...\n")
  bcv_filename_base <- "Overall_BCV_Plot"
  tryCatch({
      plotBCV(y); title("BCV Plot")
      if (save_files) {
          pdf_filename_bcv <- file.path(output_dir, paste0(bcv_filename_base, ".pdf"))
          pdf(pdf_filename_bcv); plotBCV(y); title("BCV Plot"); dev.off()
          cat("   BCVプロットを保存しました:", pdf_filename_bcv, "\n")
      }
  }, error=function(e){ warning("BCVプロットの生成/保存中にエラー:", e$message)})

  # 7. MDS Plot
  cat("7. MDSプロット生成中...\n")
  mds_filename_base <- "Overall_MDS_Plot"
  tryCatch({
      mds_data <- plotMDS(y, plot=FALSE)
      # Plot colored by group
      plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=16, main="MDS Plot (Color=Group)", xlab="Dim 1", ylab="Dim 2", las=1)
      text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7)
      legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group")
      if (save_files) {
          pdf_filename_mds_group <- file.path(output_dir, paste0(mds_filename_base, "_GroupColor.pdf"))
          pdf(pdf_filename_mds_group); plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=16, main="MDS Plot (Color=Group)", xlab="Dim 1", ylab="Dim 2", las=1); text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group"); dev.off()
          cat("   MDSプロット(Group Color)を保存しました:", pdf_filename_mds_group, "\n")
      }
       # Plot colored/shaped by batch (if available)
       if (!is.null(y$samples$batch)) {
          batch_levels <- levels(y$samples$batch)
          batch_pch <- (15:(15 + length(batch_levels) - 1)) %% 26 # Use different point shapes
          plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=batch_pch[as.numeric(y$samples$batch)], main="MDS Plot (Color=Group, Shape=Batch)", xlab="Dim 1", ylab="Dim 2", las=1)
          text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7)
          legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group")
          legend("bottomright", legend=batch_levels, pch=batch_pch, border=NA, cex=0.8, title="Batch")
          if (save_files) {
              pdf_filename_mds_batch <- file.path(output_dir, paste0(mds_filename_base, "_BatchShape.pdf"))
              pdf(pdf_filename_mds_batch); plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=batch_pch[as.numeric(y$samples$batch)], main="MDS Plot (Color=Group, Shape=Batch)", xlab="Dim 1", ylab="Dim 2", las=1); text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group"); legend("bottomright", legend=batch_levels, pch=batch_pch, border=NA, cex=0.8, title="Batch"); dev.off()
              cat("   MDSプロット(Batch Shape)を保存しました:", pdf_filename_mds_batch, "\n")
          }
       }
  }, error=function(e){ warning("MDSプロットの生成/保存中にエラー:", e$message)})

  # 8. GLM Fit
  cat("8. GLMモデルフィッティング中...\n"); fit_glm <- glmQLFit(y, design); cat("   GLMフィット完了。\n")

  # --- 比較ごとの処理 (ループ) ---
  results_output_tables <- list()
  cat("\n--- 各比較の差次発現解析を開始 ---\n")

  comp_names_in_list <- names(internal_comparison_list)
  if(is.null(comp_names_in_list)) comp_names_in_list <- 1:length(internal_comparison_list) # Use index if no names

  for (i in 1:length(internal_comparison_list)) {
      current_pair <- internal_comparison_list[[i]]
      list_element_name <- comp_names_in_list[i]
      # Use provided list name if valid, otherwise generate one
      if(is.na(list_element_name) || is.null(list_element_name) || !nzchar(list_element_name)){
           comparison_name_key <- paste0(current_pair[1], "_vs_", current_pair[2])
      } else {
           comparison_name_key <- list_element_name
      }
      target_group <- current_pair[1]; reference_group <- current_pair[2]
      safe_comparison_name <- gsub("[^a-zA-Z0-9_.-]", "_", comparison_name_key) # Sanitize name for filenames
      cat("\n解析中:", comparison_name_key, "(", i, "/", length(internal_comparison_list), ")\n")

      # 8a. コントラスト作成
      target_design_colname <- make.names(paste0(group_term_name, target_group))
      reference_design_colname <- make.names(paste0(group_term_name, reference_group))
      if (!target_design_colname %in% colnames(design) || !reference_design_colname %in% colnames(design)) {
          warning(paste0("スキップ: コントラスト作成失敗 (", comparison_name_key, ")\n   デザイン行列に `", target_design_colname, "` または `", reference_design_colname, "` が見つかりません。\n   列名: ", paste(colnames(design), collapse = ", ")))
          results_output_tables[[comparison_name_key]] <- NULL # Store NULL for this failed comparison
          next # Skip to the next comparison
      }
      contrast_formula <- paste0(target_design_colname, " - ", reference_design_colname)
      con_target_vs_ref <- makeContrasts(contrasts = contrast_formula, levels = design)
      cat("   コントラスト:\n"); print(con_target_vs_ref)

      # 8b. GLM QL F-test 実行
      qlf_result_single <- NULL # Initialize
      tryCatch({
          qlf_result_single <- glmQLFTest(fit_glm, contrast = con_target_vs_ref)
          cat("   GLM QLF Test 完了。\n")
      }, error = function(e){
          warning(paste0("GLM QLF Test中にエラー (", comparison_name_key, "): ", e$message))
          results_output_tables[[comparison_name_key]] <- NULL # Store NULL if test fails
          # Use return() within tryCatch's error handler is tricky for loop control,
          # instead we'll check qlf_result_single is NULL later before topTags
      })
      # QLF Testでエラーが発生したら、この比較の残りをスキップ
      if(is.null(qlf_result_single)) next

      # 9. 結果要約とMDプロット
      cat("   結果要約とMDプロット生成中...\n")
      dt_summary <- decideTests(qlf_result_single, p.value = fdr_threshold, lfc = logfc_threshold, adjust.method = "BH")
      cat("     結果要約 (FDR < ", fdr_threshold, ", |logFC| > ", logfc_threshold, "):\n", sep="")
      print(summary(dt_summary))
      cat("     MDプロット生成中 (有意な遺伝子を強調表示)...\n")
      tryCatch({
          plotMD(qlf_result_single, status = dt_summary, main = comparison_name_key, hl.col=c("blue","red"))
          if (logfc_threshold > 0) { abline(h = c(-logfc_threshold, logfc_threshold), col = "grey", lty = 2) }
          if (save_files) {
              pdf_filename_md <- file.path(output_dir, paste0(safe_comparison_name, "_MD_Plot_Highlighted.pdf"))
              pdf(pdf_filename_md); plotMD(qlf_result_single, status = dt_summary, main = comparison_name_key, hl.col=c("red","blue")); if (logfc_threshold > 0) { abline(h = c(-logfc_threshold, logfc_threshold), col = "grey", lty = 2) }; dev.off()
              cat("     MDプロット (強調表示版) を保存しました:", pdf_filename_md, "\n")
          }
      }, error=function(e){ warning(paste0("MDプロット(", comparison_name_key, ")の生成/保存中にエラー:", e$message))})

      # 11b. 結果を TopTags テーブルとして取得・格納
      cat("   TopTagsテーブル生成中 (全遺伝子)...\n")
      current_result_table <- NULL
      tryCatch({
          tt <- topTags(qlf_result_single, n = Inf, sort.by = "PValue", adjust.method = "BH")
          current_result_table <- tt$table
          results_output_tables[[comparison_name_key]] <- current_result_table
          cat("     - TopTagsテーブル生成完了:", comparison_name_key, "\n")
      }, error = function(e){
          warning(paste0("TopTagsテーブルの生成中にエラー (", comparison_name_key, "): ", e$message))
          results_output_tables[[comparison_name_key]] <- NULL # Store NULL if topTags fails
      })

      # 11c. 結果をRDSファイルとして保存 (qlfオブジェクト) (Optional)
      if (save_files) {
        date_prefix <- format(Sys.Date(), "%y%m%d")
        rds_filename <- file.path(output_dir, paste0(date_prefix, '_qlf_', safe_comparison_name, '.rds'))
         tryCatch({ saveRDS(qlf_result_single, file = rds_filename); cat("   QLFオブジェクトをRDSとして保存しました:", rds_filename, "\n") }, error = function(e) { warning(paste0("RDSファイル(", comparison_name_key, ")の保存中にエラー:", e$message)) })
      }

  } # ループ終了
  cat("\n--- 全ての指定された比較の解析が完了 ---\n")


  # 10. バッチ補正後のPCA (共通部分の後)
  pca_results <- NULL
  if (!is.null(y$samples$batch)) {
      cat("10. バッチ補正後のPCAプロット生成中...\n")
      pca_filename_base <- "Overall_PCA_Plot_BatchCorrected"
      tryCatch({
          logcpm <- cpm(y, log=TRUE, prior.count=3)
          cat("   limma::removeBatchEffect を実行中...\n")
          # removeBatchEffect に渡す design 行列から、バッチ項に関連する列を除外する必要があるか確認
          # -> design 引数はバッチ効果推定に使うため、フルデザインで良い。
          #    covariates 引数で調整したい他の因子（例: group）を指定することも可能。
          #    ここではシンプルに batch と design を渡す。
          corrected_logcpm <- removeBatchEffect(logcpm, batch=y$samples$batch, design=design)
          cat("   PCAを計算中 (上位遺伝子選択)...\n")
          rv <- matrixStats::rowVars(corrected_logcpm)
          # Ensure rv doesn't contain NA/NaN which crashes order()
          valid_rv <- !is.na(rv) & is.finite(rv)
          rv_filtered <- rv[valid_rv]
          if(length(rv_filtered) < 2) {
              warning("PCA用の有効な分散を持つ遺伝子が2未満のため、PCAをスキップします。")
          } else {
               select <- order(rv_filtered, decreasing=TRUE)[seq_len(min(500, length(rv_filtered)))]
               # select はフィルタリング後のインデックスなので、元のデータにマップし直す
               original_indices <- which(valid_rv)[select]
               pca_results <- prcomp(t(corrected_logcpm[original_indices,]), scale. = TRUE)
               pca_data <- pca_results$x; percent_var <- round(100 * pca_results$sdev^2 / sum(pca_results$sdev^2), 1)
               batch_levels <- levels(y$samples$batch); batch_pch <- (15:(15 + length(batch_levels) - 1)) %% 26
               plot(pca_data[,1], pca_data[,2], col = as.numeric(y$samples$group), pch = batch_pch[as.numeric(y$samples$batch)], main = "PCA Plot (Batch Corrected logCPM, Top Var Genes)", xlab = paste0("PC1: ", percent_var[1], "% variance"), ylab = paste0("PC2: ", percent_var[2], "% variance"), las = 1)
               legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), title="Group", cex=0.8, border=NA)
               legend("bottomright", legend=batch_levels, pch=batch_pch, title="Batch", cex=0.8, border=NA)
               if (save_files) {
                   pdf_filename_pca <- file.path(output_dir, paste0(pca_filename_base, ".pdf"))
                   pdf(pdf_filename_pca); plot(pca_data[,1], pca_data[,2], col = as.numeric(y$samples$group), pch = batch_pch[as.numeric(y$samples$batch)], main = "PCA Plot (Batch Corrected logCPM, Top Var Genes)", xlab = paste0("PC1: ", percent_var[1], "% variance"), ylab = paste0("PC2: ", percent_var[2], "% variance"), las = 1); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), title="Group", cex=0.8, border=NA); legend("bottomright", legend=batch_levels, pch=batch_pch, title="Batch", cex=0.8, border=NA); dev.off()
                   cat("   バッチ補正PCAプロットを保存しました:", pdf_filename_pca, "\n")
               }
          }
      }, error=function(e){ warning(paste0("バッチ補正PCAの生成/保存中にエラー:", e$message))})
  } else {
      cat("10. バッチ情報がないため、バッチ補正PCAプロットはスキップされました。\n")
  }


  # --- 返り値の決定 ---
  cat("\n--- 解析プロセス終了 ---\n")
  # Remove NULL elements from results list (comparisons that failed)
  final_results_output_tables <- Filter(Negate(is.null), results_output_tables)

  if (is_list_input) {
      cat("入力がリストだったため、TopTagsテーブルのリストを返します (失敗した比較は除外)。\n")
      if(length(final_results_output_tables) == 0) warning("全ての比較で結果テーブルの生成に失敗しました。")
      return(final_results_output_tables)
  } else {
      if (length(final_results_output_tables) == 1) {
          cat("入力が単一ペアだったため、TopTagsテーブルを直接返します。\n")
          return(final_results_output_tables[[1]])
      } else {
          # This case should ideally not be reached if input validation & processing is correct
          warning("単一ペア入力処理中に予期せぬエラーまたは失敗が発生しました。リストまたはNULLを返します。")
          return(final_results_output_tables) # Return the list (might be empty)
      }
  }
}

# --- スクリプト読み込み完了メッセージ ---
cat("関数 'perform_edgeR' が定義されました。\n")
cat("使用法: source(\"perform_edgeR.R\") してから関数を呼び出してください。\n")

# --- (オプション) スクリプトとして直接実行された場合のサンプル実行 ---
if (sys.nframe() == 0 && !interactive()) {
  cat("\n--- スクリプト直接実行時のサンプル解析開始 ---\n")
  # --- サンプルデータ作成 ---
  set.seed(123)
  sample_counts <- matrix(rnbinom(2000, mu=50, size=5), ncol=8)
  rownames(sample_counts) <- paste0("Gene", 1:250)
  colnames(sample_counts) <- paste0("Sample", 1:8)
  sample_counts_df <- as.data.frame(sample_counts)
  sample_groups <- rep(c("GroupA", "GroupB"), each = 4)
  sample_batch <- factor(rep(c("BatchX", "BatchY"), times = 4))
  sample_comp_single <- c("GroupB", "GroupA")
  # Use named list for clarity in output
  sample_comp_list <- list(B_vs_A=c("GroupB", "GroupA"), A_vs_B=c("GroupA", "GroupB"))
  cat("サンプルデータ:\n")
  print(table(Group=sample_groups, Batch=sample_batch))

  # --- サンプル実行 (単一ペア) ---
  cat("\n--- サンプル実行 (単一ペア) ---\n")
  sample_result_single_table <- NULL
  tryCatch({
       sample_result_single_table <- perform_edgeR(
            countdata = sample_counts_df,
            group = sample_groups,
            comparisons = sample_comp_single,
            batch_info = sample_batch,
            fdr_threshold=0.05, # Use default FDR
            logfc_threshold=1    # Use logFC > 1 for summary
        )
   }, error = function(e){
       cat("\nサンプル解析(単一)中にエラー:\n", e$message, "\n")
   })
   cat("\n--- サンプル結果 (単一ペア, Top 3 テーブル) ---\n")
   if (!is.null(sample_result_single_table) && is.data.frame(sample_result_single_table)) {
       print(head(sample_result_single_table, 3))
   } else {
       cat("結果なし、またはデータフレームではありません。\n")
   }

  # --- サンプル実行 (複数ペアリスト) ---
  cat("\n--- サンプル実行 (複数ペアリスト) ---\n")
  sample_result_table_list <- NULL
  tryCatch({
       sample_result_table_list <- perform_edgeR(
           countdata = sample_counts_df,
           group = sample_groups,
           comparisons = sample_comp_list,
           batch_info = sample_batch,
           fdr_threshold=0.05,
           logfc_threshold = 0 # Default logFC for summary
        )
   }, error = function(e){
       cat("\nサンプル解析(複数)中にエラー:\n", e$message, "\n")
   })
   cat("\n--- サンプル結果 (複数ペアリスト) ---\n")
   if (!is.null(sample_result_table_list) && is.list(sample_result_table_list) && length(sample_result_table_list)>0) {
       cat("結果リストの名前:\n")
       print(names(sample_result_table_list))
       cat("\n最初の比較 (", names(sample_result_table_list)[1], ") の Top 3 テーブル:\n")
       if(!is.null(sample_result_table_list[[1]]) && is.data.frame(sample_result_table_list[[1]])) {
           print(head(sample_result_table_list[[1]], 3))
       } else {
            cat("最初の比較の結果テーブルがありません。\n")
       }
   } else {
       cat("結果なし、またはリストではありません。\n")
   }

  cat("\n--- サンプル解析終了 ---\n")
}
