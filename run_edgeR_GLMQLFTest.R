# 必要なライブラリのロード（関数利用前にロードされていることを前提とします）
# library(edgeR)
# library(dplyr) # サンプルデータ作成の際に select を使う場合など

#' edgeR を用いた差次発現遺伝子解析 (GLM QLF Test)
#'
#' @param counts_df データフレーム。行が遺伝子、列がサンプルに対応する生のカウントデータ。
#'                  行名に遺伝子IDが設定されていることを想定。
#' @param group_vector 文字列ベクトル。各サンプルが属するグループを示す。
#'                     要素数は counts_df の列数と一致する必要がある。
#' @param comparison_pair 文字列ベクトル。比較したい2つのグループ名を指定 (例: c("control", "treatment"))。
#'                        最初の要素が対照群、2番目の要素が処理群となる (treatment vs control の logFC を得る場合)。
#'
#' @return topTags で得られる、差次発現遺伝子の情報を含むデータフレーム (F統計量などを含む)。
#'
#' @examples
#' \dontrun{
#' # ダミーデータの作成
#' set.seed(123)
#' counts_example <- data.frame(
#'   Gene1 = rnbinom(6, mu=100, size=10),
#'   Gene2 = rnbinom(6, mu=50, size=10),
#'   Gene3 = rnbinom(6, mu=200, size=10),
#'   Gene4 = rnbinom(6, mu=10, size=10),
#'   Gene5 = rnbinom(6, mu=5, size=10)
#' )
#' counts_example <- t(counts_example) # 行を遺伝子、列をサンプルにする
#' colnames(counts_example) <- paste0("Sample", 1:6)
#' rownames(counts_example) <- paste0("Gene", 1:5)
#'
#' group_info <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
#' comparison <- c("Control", "Treatment") # Treatment vs Control
#'
#' # 関数の実行
#' results_glm <- run_edgeR_GLMQLFTest(counts_df = as.data.frame(counts_example),
#'                                     group_vector = group_info,
#'                                     comparison_pair = comparison)
#' print(head(results_glm))
#' # 結果には 'F' (F統計量) カラムが含まれます。
#' }
run_edgeR_GLMQLFTest <- function(counts_df, group_vector, comparison_pair) {

  # --- 引数のチェック ---
  if (!is.data.frame(counts_df)) {
    stop("counts_df はデータフレームである必要があります。")
  }
  if (ncol(counts_df) != length(group_vector)) {
    stop("counts_df の列数と group_vector の長さが一致しません。")
  }
  if (length(comparison_pair) != 2) {
    stop("comparison_pair は2つのグループ名を含むベクトルである必要があります。")
  }
  # comparison_pairの要素がgroup_vectorに存在するか確認
  if (!all(comparison_pair %in% unique(as.character(group_vector)))) {
    stop("comparison_pair で指定されたグループ名が group_vector に存在しません。")
  }
  # group_vectorに2種類以上のグループが含まれているか確認 (GLMで比較するため)
  if (length(unique(as.character(group_vector))) < 2) {
    stop("group_vector には少なくとも2種類のグループが含まれている必要があります。")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("パッケージ 'edgeR' がインストールされていません。インストールしてください。")
  }

  # --- 1. カウントデータの行列への変換 ---
  count_data_matrix <- as.matrix(counts_df)

  # --- 2. グループ情報の定義 ---
  # group_vector を factor に変換。levelsは比較のベースラインを意識して設定。
  # comparison_pair[1] がベースライン (対照群) となるようにする。
  group <- factor(group_vector, levels = comparison_pair)
  
  # グループ情報が2レベルになっているか最終確認 (comparison_pairに含まれないグループがフィルタリングで消える可能性など)
  # DGEList作成後に再度確認する方がよいが、ここでも簡易チェック
  if (length(levels(droplevels(group))) < 2) {
       stop(paste("group_vectorから comparison_pair に基づいて有効なグループを作成した結果、比較可能な2つのグループになりませんでした。現在の有効なグループ:", paste(levels(droplevels(group)), collapse=", ")))
  }


  # --- 3. DGEListオブジェクトの作成 ---
  d <- edgeR::DGEList(counts = count_data_matrix, group = group)
  cat("DGEListオブジェクト作成完了。\n")
  cat("サンプル数:", ncol(d), "初期遺伝子数:", nrow(d), "\n")
  cat("グループ情報 (DGEListから):", levels(d$samples$group), "\n")


  # --- 4. 低発現遺伝子のフィルタリング ---
  # デザイン行列をフィルタリングの前に定義（filterByExprがdesignを要求する場合があるため）
  # ただし、このシンプルな2群比較では group のみで十分なことが多い
  design <- model.matrix(~group, data = d$samples)
  cat("デザイン行列作成 (フィルタリング前):\n")
  print(design)

  keep <- edgeR::filterByExpr(d, design = design) # group引数の代わりにdesign引数も使える
  d <- d[keep, , keep.lib.sizes = FALSE]
  cat("フィルタリング後の遺伝子数:", nrow(d), "\n")

  if (nrow(d) == 0) {
    stop("フィルタリングの結果、解析対象の遺伝子が0になりました。フィルタリング基準を確認してください。")
  }
  # フィルタリング後に再度groupのレベルを確認 (サンプルが全て落ちたグループがないか)
  d$samples$group <- droplevels(d$samples$group)
  if (length(levels(d$samples$group)) < 2) {
    stop(paste("フィルタリングの結果、比較可能な2つのグループが残存しませんでした。残存グループ:", paste(levels(d$samples$group), collapse=", ")))
  }


  # --- 5. 正規化 (TMM) ---
  d <- edgeR::calcNormFactors(d)
  cat("TMM正規化完了。\n")

  # --- 6. 分散の推定 (GLM用) ---
  # デザイン行列を再作成 (フィルタリング後のd$samplesを使うため)
  design <- model.matrix(~group, data = d$samples)
  cat("デザイン行列作成 (分散推定用):\n")
  print(design)
  if (ncol(design) < 2) {
      stop("デザイン行列の列が2未満です。通常、2群比較ではインターセプトとグループ効果の列が必要です。グループ設定を確認してください。")
  }


  d <- edgeR::estimateDisp(d, design, robust = TRUE) # robust=TRUE は外れ値に頑健
  cat("分散推定 (estimateDisp) 完了。\n")
  # (オプション) BCVプロットで分散推定を確認
  # plotBCV(d)

  # --- 7. モデルフィッティング (GLM QLF Test) ---
  fit <- edgeR::glmQLFit(d, design, robust = TRUE)
  cat("モデルフィッティング (glmQLFit) 完了。\n")

  # --- 8. 差次発現検定 (GLM QLF Test) ---
  # comparison_pair[2] (処理群) の効果を検定する。
  # design行列の列名を確認し、適切な係数を指定する。
  # model.matrix(~group) で group の levels が c("A", "B") の場合、
  # design の列は "(Intercept)" と "groupB" になる。
  # "groupB" が A に対する B の効果を示す。
  # coef引数には、この効果を示す列のインデックスまたは名前を指定。
  # comparison_pair[1]がベースラインなので、comparison_pair[2]に対応する係数を検定する。
  # デザイン行列の2列目がこの比較に該当する(最初の列はインターセプト)。
  coef_to_test <- ncol(design) # 通常はデザイン行列の最後の列
  
  # より頑健な方法として、列名を直接指定する
  # factor 'group' の2番目のレベル名 (comparison_pair[2]) を使って列名を生成
  target_colname <- paste0("group", comparison_pair[2]) # levels(group)[2] を使うとより動的
  if (!target_colname %in% colnames(design)) {
      cat("警告: 想定される係数名 '", target_colname, "' がデザイン行列の列名にありません。\n")
      cat("利用可能な列名:", paste(colnames(design), collapse=", "), "\n")
      cat("最後の係数 (列インデックス:", coef_to_test, ") を使って検定を試みます。\n")
      if (coef_to_test < 2) {
          stop("検定対象の係数を特定できませんでした。デザイン行列とcomparison_pairを確認してください。")
      }
  } else {
      coef_to_test <- target_colname
      cat("検定対象の係数として '", coef_to_test, "' を使用します。\n")
  }


  qlf_test <- edgeR::glmQLFTest(fit, coef = coef_to_test)
  cat("差次発現検定 (glmQLFTest) 実行完了。\n")

  # --- 9. 結果の取得とデータフレームへの変換 ---
  # フィルタリング後の全遺伝子の結果を取得
  result_table <- edgeR::topTags(qlf_test, n = nrow(d$counts), sort.by = "PValue")
  result_df <- as.data.frame(result_table)
  cat("結果テーブル作成完了。\n")

  # 結果テーブルにF統計量が含まれていることを確認
  if (!"F" %in% colnames(result_df)) {
    cat("警告: 結果テーブルにF統計量カラムが含まれていません。\n")
  }

  return(result_df)
}
