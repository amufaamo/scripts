#!/usr/bin/env Rscript

# ==============================================================================
# R Script: perform_edgeR.R
# (Description, Usage, Requirements - 変更なし)
# Last Updated: 2025-04-08 (コントラスト作成ロジック修正)
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


# --- 3. edgeR解析を実行する関数 ---
perform_edgeR <- function(
    countdata,
    group,
    comparisons,
    batch_info = NULL,
    output_dir = NULL,
    use_all_samples_for_qc_norm_disp = TRUE,
    fdr_threshold = 0.05,
    logfc_threshold = 0
    ) {

  # --- 入力データの検証 ---
  cat("--- 入力データの検証を開始 ---\n")
  # ... (変更なし) ...
  if (!is.data.frame(countdata) && !is.matrix(countdata)) { stop("エラー: countdata は data.frame または matrix である必要があります。") }
  if (!is.vector(group) && !is.factor(group)) { stop("エラー: group は vector または factor である必要があります。") }
  if (ncol(countdata) != length(group)) { stop("エラー: カウントデータの列数とグループ情報のサンプル数が一致しません。") }
  original_group_factor <- factor(group)
  if (!is.numeric(fdr_threshold) || fdr_threshold < 0 || fdr_threshold > 1) { stop("エラー: fdr_threshold は 0 から 1 の数値である必要があります。") }
  if (!is.numeric(logfc_threshold) || logfc_threshold < 0) { stop("エラー: logfc_threshold は 0 以上の数値である必要があります。") }
  cat(paste0("使用するFDR閾値: ", fdr_threshold, "\n"))
  cat(paste0("使用するlogFC閾値 (絶対値): ", logfc_threshold, "\n"))
  is_list_input <- FALSE; internal_comparison_list <- list()
  if (is.list(comparisons) && !is.data.frame(comparisons)) { if(length(comparisons) == 0) stop("エラー: comparisons リストが空です。"); is_list_input <- TRUE; internal_comparison_list <- comparisons; cat("複数の比較ペアがリストとして指定されました。\n") } else if (is.character(comparisons) && length(comparisons) == 2) { is_list_input <- FALSE; internal_comparison_list <- list(comparisons); cat("単一の比較ペアがベクトルとして指定されました。\n") } else { stop("エラー: 'comparisons' 引数は、c('target', 'reference')形式のベクトル、またはそのリストである必要があります。") }
  valid_pairs_exist <- FALSE; processed_list <- list(); list_names <- names(internal_comparison_list); valid_names <- c()
  for(i in 1:length(internal_comparison_list)) { pair <- internal_comparison_list[[i]]; pair_name_raw <- if(!is.null(list_names) && nzchar(list_names[i])) list_names[i] else paste(pair, collapse="_vs_"); if(is.character(pair) && length(pair) == 2 && pair[1] != pair[2] && all(pair %in% levels(original_group_factor))) { valid_pairs_exist <- TRUE; processed_list[[length(processed_list) + 1]] <- pair; valid_names <- c(valid_names, pair_name_raw); if(is_list_input) cat("  - OK:", paste(pair, collapse=" vs "), "( Name:", pair_name_raw, ")\n") } else { warning("リスト内の無効なペア:", paste(pair, collapse=" / "), " - スキップされます。") } }; internal_comparison_list <- processed_list; names(internal_comparison_list) <- valid_names; if (length(internal_comparison_list) == 0) { stop("エラー: 有効な比較ペアが指定されませんでした。") }
  analysis_batch <- NULL; original_batch_factor <- NULL; if (!is.null(batch_info)) { if (length(batch_info) != ncol(countdata)) { stop("エラー: バッチ情報の長さがカウントデータの列数と一致しません。") }; original_batch_factor <- factor(batch_info); cat("バッチ情報が提供されました。\n") } else { cat("バッチ情報は提供されていません。\n") }
  if (!is.null(original_batch_factor)) { cat("\n--- グループとバッチのクロス集計 ---\n"); if(length(original_group_factor) == length(original_batch_factor)) { cross_tab <- table(Group = original_group_factor, Batch = original_batch_factor); print(cross_tab); if(any(cross_tab == 0)){ warning("警告: グループとバッチの組み合わせにサンプル数が0のセルがあります。交絡の可能性があります。") } } else { warning("警告: クロス集計スキップ。グループとバッチの長さが一致しません。") }; cat("---------------------------------\n\n") }

  # --- 使用するサンプルを決定 ---
  # ... (変更なし) ...
  if (use_all_samples_for_qc_norm_disp) { cat("初期ステップで「全サンプル」を使用します。\n"); analysis_countdata <- countdata; analysis_group <- original_group_factor; if (!is.null(original_batch_factor)) { analysis_batch <- original_batch_factor } } else { cat("初期ステップで「比較ペアのサンプルのみ」を使用します (注意: 複数比較リストの場合、最初の有効なペアに基づきます)。\n"); current_pair_for_subset <- internal_comparison_list[[1]]; samples_to_keep_logical <- original_group_factor %in% current_pair_for_subset; if (sum(samples_to_keep_logical) < 2) { stop("エラー: 比較ペアに属するサンプル数が2未満のため、ペアのみでの解析は実行できません。") }; analysis_countdata <- countdata[, samples_to_keep_logical, drop = FALSE]; analysis_group <- droplevels(original_group_factor[samples_to_keep_logical]); if (nlevels(analysis_group) < 2) { stop("エラー: サブセット後に比較に必要な2つのグループが存在しません。") }; if (!is.null(original_batch_factor)) { analysis_batch <- droplevels(original_batch_factor[samples_to_keep_logical]); if(nlevels(analysis_batch) <= 1) { analysis_batch <- NULL } }; cat("サブセット後のサンプル数:", ncol(analysis_countdata), "\n"); cat("サブセット後のグループ情報:\n"); print(table(analysis_group)); if (!is.null(analysis_batch)) { cat("サブセット後のバッチ情報:\n"); print(table(analysis_batch)) } else { cat("サブセット後のバッチ情報: なし\n")} }
  cat("解析に使用するグループ情報:\n"); print(table(analysis_group)); if (!is.null(analysis_batch)) { cat("解析に使用するバッチ情報:\n"); print(table(analysis_batch)) }
  save_files <- !is.null(output_dir); if (save_files) { if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE); cat("出力ディレクトリを作成しました:", output_dir, "\n") }; cat("結果ファイルは次のディレクトリに保存されます:", output_dir, "\n") } else { cat("ファイル保存はスキップされます。\n") }

  # --- edgeR 解析ステップ (共通部分) ---
  # ... (Step 1 DGEList から Step 8 GLM Fit までは変更なし) ...
  cat("\n--- edgeR 解析ステップ開始 (共通部分) ---\n"); y <- DGEList(counts = analysis_countdata, group = analysis_group); if (!is.null(analysis_batch)) { y$samples$batch <- analysis_batch }; cat("1. DGEList作成完了。\n"); group_term_name <- "analysis_group"; batch_term_name <- "analysis_batch"; if (!is.null(analysis_batch)) { design_formula <- as.formula(paste("~ 0 +", group_term_name, "+", batch_term_name)); model_data <- data.frame(analysis_group=analysis_group); model_data[[batch_term_name]] <- analysis_batch; design <- model.matrix(design_formula, data=model_data); cat("2. デザイン行列を 'group + batch' で作成。\n") } else { design_formula <- as.formula(paste("~ 0 +", group_term_name)); model_data <- data.frame(analysis_group=analysis_group); design <- model.matrix(design_formula, data=model_data); colnames(design) <- levels(analysis_group); cat("2. デザイン行列を 'group' のみで作成。\n") }; colnames(design) <- make.names(colnames(design)); print("デザイン行列:"); print(design); design_rank <- qr(design)$rank; num_columns <- ncol(design); if (design_rank < num_columns) { warning(paste("警告: デザイン行列がフルランクではありません (ランク =", design_rank, ", 列数 =", num_columns, ").")) } else { cat("   デザイン行列はフルランクです。\n") }; cat("3. フィルタリング中...\n"); keep <- filterByExpr(y, design = design); if(sum(keep) == 0) { stop("エラー: フィルタリング後に遺伝子が残りませんでした。") }; y <- y[keep, , keep.lib.sizes = FALSE]; cat("   完了。残りの遺伝子数:", nrow(y), "\n"); cat("4. TMM正規化中...\n"); y <- calcNormFactors(y); cat("   完了。\n"); cat("5. 分散推定中...\n"); y <- tryCatch({ estimateDisp(y, design) }, error = function(e) { stop(paste("分散推定中にエラー:", e$message)) }); cat("   完了。\n"); cat("6. BCVプロット生成中...\n"); bcv_filename_base <- "Overall_BCV_Plot"; tryCatch({ plotBCV(y); title("BCV Plot"); if (save_files) { pdf_filename_bcv <- file.path(output_dir, paste0(bcv_filename_base, ".pdf")); pdf(pdf_filename_bcv); plotBCV(y); title("BCV Plot"); dev.off(); cat("   BCVプロットを保存しました:", pdf_filename_bcv, "\n") } }, error=function(e){ warning("BCVプロットの生成/保存中にエラー:", e$message)}); cat("7. MDSプロット生成中...\n"); mds_filename_base <- "Overall_MDS_Plot"; tryCatch({ mds_data <- plotMDS(y, plot=FALSE); plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=16, main="MDS Plot (Color=Group)", xlab="Dim 1", ylab="Dim 2", las=1); text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group"); if (save_files) { pdf_filename_mds_group <- file.path(output_dir, paste0(mds_filename_base, "_GroupColor.pdf")); pdf(pdf_filename_mds_group); plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=16, main="MDS Plot (Color=Group)", xlab="Dim 1", ylab="Dim 2", las=1); text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group"); dev.off(); cat("   MDSプロット(Group Color)を保存しました:", pdf_filename_mds_group, "\n") }; if (!is.null(y$samples$batch)) { batch_levels <- levels(y$samples$batch); batch_pch <- (15:(15 + length(batch_levels) - 1)) %% 26; plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=batch_pch[as.numeric(y$samples$batch)], main="MDS Plot (Color=Group, Shape=Batch)", xlab="Dim 1", ylab="Dim 2", las=1); text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group"); legend("bottomright", legend=batch_levels, pch=batch_pch, border=NA, cex=0.8, title="Batch"); if (save_files) { pdf_filename_mds_batch <- file.path(output_dir, paste0(mds_filename_base, "_BatchShape.pdf")); pdf(pdf_filename_mds_batch); plot(mds_data$x, mds_data$y, col=as.numeric(y$samples$group), pch=batch_pch[as.numeric(y$samples$batch)], main="MDS Plot (Color=Group, Shape=Batch)", xlab="Dim 1", ylab="Dim 2", las=1); text(mds_data$x, mds_data$y, labels=colnames(y), pos=3, cex=0.7); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), border=NA, cex=0.8, title="Group"); legend("bottomright", legend=batch_levels, pch=batch_pch, border=NA, cex=0.8, title="Batch"); dev.off(); cat("   MDSプロット(Batch Shape)を保存しました:", pdf_filename_mds_batch, "\n") } } }, error=function(e){ warning("MDSプロットの生成/保存中にエラー:", e$message)}); cat("8. GLMモデルフィッティング中...\n"); fit_glm <- glmQLFit(y, design); cat("   GLMフィット完了。\n")

  # --- 比較ごとの処理 (ループ) ---
  results_output_tables <- list()
  cat("\n--- 各比較の差次発現解析を開始 ---\n")
  comp_names_in_list <- names(internal_comparison_list); if(is.null(comp_names_in_list)) comp_names_in_list <- 1:length(internal_comparison_list)
  for (i in 1:length(internal_comparison_list)) {
      current_pair <- internal_comparison_list[[i]]; list_element_name <- comp_names_in_list[i]; if(is.na(list_element_name) || is.null(list_element_name) || !nzchar(list_element_name)){ comparison_name_key <- paste0(current_pair[1], "_vs_", current_pair[2]) } else { comparison_name_key <- list_element_name }; target_group <- current_pair[1]; reference_group <- current_pair[2]; safe_comparison_name <- gsub("[^a-zA-Z0-9_.-]", "_", comparison_name_key); cat("\n解析中:", comparison_name_key, "(", i, "/", length(internal_comparison_list), ")\n")

      # 8a. コントラスト作成 (★修正箇所)
      if (!is.null(analysis_batch)) { # バッチ情報がモデルに含まれているかチェック
          # バッチあり -> "analysis_group" + "レベル名" の形式
          target_design_colname <- make.names(paste0(group_term_name, target_group))
          reference_design_colname <- make.names(paste0(group_term_name, reference_group))
      } else {
          # バッチなし -> "レベル名" そのもの
          target_design_colname <- make.names(target_group)
          reference_design_colname <- make.names(reference_group)
      }

      # 作成した列名が design 行列に存在するか確認
      if (!target_design_colname %in% colnames(design) || !reference_design_colname %in% colnames(design)) {
          warning(paste0("スキップ: コントラスト作成失敗 (", comparison_name_key, ")\n   デザイン行列に `",
                         target_design_colname, "` または `", reference_design_colname,
                         "` が見つかりません。\n   実際の列名: ",
                         paste(colnames(design), collapse = ", ")))
          results_output_tables[[comparison_name_key]] <- NULL # 結果リストにはNULLを入れる
          next # 次の比較へ
      }

      # 正しい列名を使ってコントラスト文字列を作成
      contrast_formula <- paste0(target_design_colname, " - ", reference_design_colname)
      con_target_vs_ref <- makeContrasts(contrasts = contrast_formula, levels = design)
      cat("   コントラスト:\n"); print(con_target_vs_ref)

      # 8b. GLM QL F-test 実行
      qlf_result_single <- NULL; tryCatch({ qlf_result_single <- glmQLFTest(fit_glm, contrast = con_target_vs_ref); cat("   GLM QLF Test 完了。\n") }, error = function(e){ warning(paste0("GLM QLF Test中にエラー (", comparison_name_key, "): ", e$message)); results_output_tables[[comparison_name_key]] <- NULL })
      if(is.null(qlf_result_single)) next

      # 9. 結果要約とMDプロット
      # ... (変更なし) ...
      cat("   結果要約とMDプロット生成中...\n"); dt_summary <- decideTests(qlf_result_single, p.value = fdr_threshold, lfc = logfc_threshold, adjust.method = "BH"); cat("     結果要約 (FDR < ", fdr_threshold, ", |logFC| > ", logfc_threshold, "):\n", sep=""); print(summary(dt_summary)); cat("     MDプロット生成中 (有意な遺伝子を強調表示)...\n"); tryCatch({ plotMD(qlf_result_single, status = dt_summary, main = comparison_name_key, hl.col=c("blue","red")); if (logfc_threshold > 0) { abline(h = c(-logfc_threshold, logfc_threshold), col = "grey", lty = 2) }; if (save_files) { pdf_filename_md <- file.path(output_dir, paste0(safe_comparison_name, "_MD_Plot_Highlighted.pdf")); pdf(pdf_filename_md); plotMD(qlf_result_single, status = dt_summary, main = comparison_name_key, hl.col=c("blue","red")); if (logfc_threshold > 0) { abline(h = c(-logfc_threshold, logfc_threshold), col = "grey", lty = 2) }; dev.off(); cat("     MDプロット (強調表示版) を保存しました:", pdf_filename_md, "\n") } }, error=function(e){ warning(paste0("MDプロット(", comparison_name_key, ")の生成/保存中にエラー:", e$message))})

      # 11b. 結果を TopTags テーブルとして取得・格納
      # ... (変更なし) ...
      cat("   TopTagsテーブル生成中 (全遺伝子)...\n"); current_result_table <- NULL; tryCatch({ tt <- topTags(qlf_result_single, n = Inf, sort.by = "PValue", adjust.method = "BH"); current_result_table <- tt$table; results_output_tables[[comparison_name_key]] <- current_result_table; cat("     - TopTagsテーブル生成完了:", comparison_name_key, "\n") }, error = function(e){ warning(paste0("TopTagsテーブルの生成中にエラー (", comparison_name_key, "): ", e$message)); results_output_tables[[comparison_name_key]] <- NULL })

      # 11c. 結果をRDSファイルとして保存 (qlfオブジェクト) (Optional)
      # ... (変更なし) ...
       if (save_files) { date_prefix <- format(Sys.Date(), "%y%m%d"); rds_filename <- file.path(output_dir, paste0(date_prefix, '_qlf_', safe_comparison_name, '.rds')); tryCatch({ saveRDS(qlf_result_single, file = rds_filename); cat("   QLFオブジェクトをRDSとして保存しました:", rds_filename, "\n") }, error = function(e) { warning(paste0("RDSファイル(", comparison_name_key, ")の保存中にエラー:", e$message)) }) }

  } # ループ終了
  cat("\n--- 全ての指定された比較の解析が完了 ---\n")

  # 10. バッチ補正後のPCA (共通部分の後)
  # ... (変更なし) ...
  pca_results <- NULL; if (!is.null(y$samples$batch)) { cat("10. バッチ補正後のPCAプロット生成中...\n"); pca_filename_base <- "Overall_PCA_Plot_BatchCorrected"; tryCatch({ logcpm <- cpm(y, log=TRUE, prior.count=3); cat("   limma::removeBatchEffect を実行中...\n"); corrected_logcpm <- removeBatchEffect(logcpm, batch=y$samples$batch, design=design); cat("   PCAを計算中 (上位遺伝子選択)...\n"); rv <- matrixStats::rowVars(corrected_logcpm); valid_rv <- !is.na(rv) & is.finite(rv); rv_filtered <- rv[valid_rv]; if(length(rv_filtered) < 2) { warning("PCA用の有効な分散を持つ遺伝子が2未満のため、PCAをスキップします。") } else { select <- order(rv_filtered, decreasing=TRUE)[seq_len(min(500, length(rv_filtered)))]; original_indices <- which(valid_rv)[select]; pca_results <- prcomp(t(corrected_logcpm[original_indices,]), scale. = TRUE); pca_data <- pca_results$x; percent_var <- round(100 * pca_results$sdev^2 / sum(pca_results$sdev^2), 1); batch_levels <- levels(y$samples$batch); batch_pch <- (15:(15 + length(batch_levels) - 1)) %% 26; plot(pca_data[,1], pca_data[,2], col = as.numeric(y$samples$group), pch = batch_pch[as.numeric(y$samples$batch)], main = "PCA Plot (Batch Corrected logCPM, Top Var Genes)", xlab = paste0("PC1: ", percent_var[1], "% variance"), ylab = paste0("PC2: ", percent_var[2], "% variance"), las = 1); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), title="Group", cex=0.8, border=NA); legend("bottomright", legend=batch_levels, pch=batch_pch, title="Batch", cex=0.8, border=NA); if (save_files) { pdf_filename_pca <- file.path(output_dir, paste0(pca_filename_base, ".pdf")); pdf(pdf_filename_pca); plot(pca_data[,1], pca_data[,2], col = as.numeric(y$samples$group), pch = batch_pch[as.numeric(y$samples$batch)], main = "PCA Plot (Batch Corrected logCPM, Top Var Genes)", xlab = paste0("PC1: ", percent_var[1], "% variance"), ylab = paste0("PC2: ", percent_var[2], "% variance"), las = 1); legend("topright", legend=levels(y$samples$group), fill=1:nlevels(y$samples$group), title="Group", cex=0.8, border=NA); legend("bottomright", legend=batch_levels, pch=batch_pch, title="Batch", cex=0.8, border=NA); dev.off(); cat("   バッチ補正PCAプロットを保存しました:", pdf_filename_pca, "\n") } } }, error=function(e){ warning(paste0("バッチ補正PCAの生成/保存中にエラー:", e$message))}) } else { cat("10. バッチ情報がないため、バッチ補正PCAプロットはスキップされました。\n") }

  # --- 返り値の決定 ---
  # ... (変更なし) ...
  cat("\n--- 解析プロセス終了 ---\n"); final_results_output_tables <- Filter(Negate(is.null), results_output_tables); if (is_list_input) { cat("入力がリストだったため、TopTagsテーブルのリストを返します (失敗した比較は除外)。\n"); if(length(final_results_output_tables) == 0) warning("全ての比較で結果テーブルの生成に失敗しました。"); return(final_results_output_tables) } else { if (length(final_results_output_tables) == 1) { cat("入力が単一ペアだったため、TopTagsテーブルを直接返します。\n"); return(final_results_output_tables[[1]]) } else { warning("単一ペア入力処理中に予期せぬエラーまたは失敗が発生しました。リストまたはNULLを返します。"); return(final_results_output_tables) } }
}

# --- スクリプト読み込み完了メッセージ ---
cat("関数 'perform_edgeR' が定義されました。\n")
cat("使用法: source(\"perform_edgeR.R\") してから関数を呼び出してください。\n")

# --- (オプション) スクリプトとして直接実行された場合のサンプル実行 ---
# ... (変更なし) ...
if (sys.nframe() == 0 && !interactive()) { cat("\n--- スクリプト直接実行時のサンプル解析開始 ---\n"); set.seed(123); sample_counts <- matrix(rnbinom(2000, mu=50, size=5), ncol=8); rownames(sample_counts) <- paste0("Gene", 1:250); colnames(sample_counts) <- paste0("Sample", 1:8); sample_counts_df <- as.data.frame(sample_counts); sample_groups <- rep(c("GroupA", "GroupB"), each = 4); sample_batch <- factor(rep(c("BatchX", "BatchY"), times = 4)); sample_comp_single <- c("GroupB", "GroupA"); sample_comp_list <- list(Comp1=c("GroupB", "GroupA"), Comp2=c("GroupA", "GroupB")); cat("サンプルデータ:\n"); print(table(Group=sample_groups, Batch=sample_batch)); cat("\n--- サンプル実行 (単一ペア) ---\n"); sample_result_single_table <- NULL; tryCatch({ sample_result_single_table <- perform_edgeR( countdata = sample_counts_df, group = sample_groups, comparisons = sample_comp_single, batch_info = sample_batch, fdr_threshold=0.1, logfc_threshold=1 ) }, error = function(e){ cat("\nサンプル解析(単一)中にエラー:\n", e$message, "\n") }); cat("\n--- サンプル結果 (単一ペア, Top 3 テーブル) ---\n"); if (!is.null(sample_result_single_table) && is.data.frame(sample_result_single_table)) { print(head(sample_result_single_table, 3)) } else { cat("結果なし\n") }; cat("\n--- サンプル実行 (複数ペアリスト) ---\n"); sample_result_table_list <- NULL; tryCatch({ sample_result_table_list <- perform_edgeR( countdata = sample_counts_df, group = sample_groups, comparisons = sample_comp_list, batch_info = sample_batch ) }, error = function(e){ cat("\nサンプル解析(複数)中にエラー:\n", e$message, "\n") }); cat("\n--- サンプル結果 (複数ペアリスト) ---\n"); if (!is.null(sample_result_table_list) && is.list(sample_result_table_list) && length(sample_result_table_list)>0) { cat("結果リストの名前:\n"); print(names(sample_result_table_list)); cat("\n最初の比較 (", names(sample_result_table_list)[1], ") の Top 3 テーブル:\n"); if(!is.null(sample_result_table_list[[1]]) && is.data.frame(sample_result_table_list[[1]])) { print(head(sample_result_table_list[[1]], 3)) } else { cat("最初の比較の結果テーブルがありません。\n") } } else { cat("結果なし\n") }; cat("\n--- サンプル解析終了 ---\n") }
