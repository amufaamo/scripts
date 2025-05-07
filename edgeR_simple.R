# 必要なライブラリのロード（関数内で明示的にロードするか、事前にロードされていることを前提とする）
# この例では、関数利用前にライブラリがロードされていることを前提とします。
# library(dplyr) # select関数を使う場合は必要だが、関数内では直接使わない形にする
# library(edgeR)

#' edgeR を用いた差次発現遺伝子解析 (exactTest)
#'
#' @param counts_df データフレーム。行が遺伝子、列がサンプルに対応する生のカウントデータ。
#'                  行名に遺伝子IDが設定されていることを想定。
#' @param group_vector 文字列ベクトル。各サンプルが属するグループを示す。
#'                     要素数は counts_df の列数と一致する必要がある。
#' @param comparison_pair 文字列ベクトル。比較したい2つのグループ名を指定 (例: c("control", "treatment"))。
#'                        最初の要素が対照群、2番目の要素が処理群となる (treatment vs control の logFC を得る場合)。
#'
#' @return topTags で得られる、差次発現遺伝子の情報を含むデータフレーム。
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
#' comparison <- c("Control", "Treatment")
#'
#' # 関数の実行
#' results <- run_edgeR_exactTest(counts_df = as.data.frame(counts_example), # データフレームとして渡す
#'                               group_vector = group_info,
#'                               comparison_pair = comparison)
#' print(head(results))
#' }
run_edgeR_exactTest <- function(counts_df, group_vector, comparison_pair) {

  # 引数のチェック (簡易版)
  if (!is.data.frame(counts_df)) {
    stop("counts_df はデータフレームである必要があります。")
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
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("パッケージ 'edgeR' がインストールされていません。インストールしてください。")
  }


  # 1. カウントデータの行列への変換
  # 関数内では、入力データフレームが行名に遺伝子IDを持つ数値データであると仮定
  count_data_matrix <- as.matrix(counts_df)

  # 2. グループ情報の定義
  # group_vector を factor に変換し、比較の順番を comparison_pair に合わせる
  group <- factor(group_vector, levels = comparison_pair)

  # 3. DGEListオブジェクトの作成
  d <- edgeR::DGEList(counts = count_data_matrix, group = group)
  cat("DGEListオブジェクト作成完了。\n")
  cat("サンプル数:", ncol(d), "初期遺伝子数:", nrow(d), "\n")

  # 4. 低発現遺伝子のフィルタリング
  # filterByExpr は group 情報を使ってフィルタリング閾値を自動で決定
  keep <- edgeR::filterByExpr(d, group = group)
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
  # トレンド分散の推定
  d <- edgeR::estimateTrendedDisp(d)
  cat("トレンド分散推定完了。\n")
  # タグワイズ分散の推定
  d <- edgeR::estimateTagwiseDisp(d)
  cat("タグワイズ分散推定完了。\n")

  # 7. 差次発現検定 (exactTest)
  # comparison_pair の順番で比較 (例: comparison_pair = c("control", "treatment") なら treatment vs control)
  # exactTest の pair 引数は levels(d$samples$group)[1:2] をデフォルトで使うので、
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

# --- 関数の使い方 (コメントアウト解除して実行可能) ---
# # ダミーデータの準備
# set.seed(123) # 再現性のためのシード設定
#
# # counts_df: 行が遺伝子、列がサンプル
# counts_example_df <- data.frame(
#   Sample1_Ctrl = rnbinom(100, mu=50, size=10),
#   Sample2_Ctrl = rnbinom(100, mu=60, size=10),
#   Sample3_Ctrl = rnbinom(100, mu=55, size=10),
#   Sample4_Treat = rnbinom(100, mu=150, size=10),
#   Sample5_Treat = rnbinom(100, mu=160, size=10),
#   Sample6_Treat = rnbinom(100, mu=155, size=10)
# )
# rownames(counts_example_df) <- paste0("Gene", 1:100)
#
# # group_vector: 各サンプルが属するグループ
# group_info_vector <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
#
# # comparison_pair: 比較するグループのペア (Treatment vs Control を見たい場合)
# comparison_groups <- c("Control", "Treatment") # 最初の要素が対照
#
# # 関数の実行
# if (requireNamespace("edgeR", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
#   cat("必要なパッケージがロードされています。\n")
#   tryCatch({
#     deg_results <- run_edgeR_exactTest(
#       counts_df = counts_example_df,
#       group_vector = group_info_vector,
#       comparison_pair = comparison_groups
#     )
#     # 結果の表示
#     print(head(deg_results))
#
#     # logFC が Treatment / Control を示しているか確認
#     # Treatment 群の平均カウントが高い遺伝子で logFC が正になるはず
#     # 例えば、Gene1 の平均: Control群 vs Treatment群
#     cat("\nGene1 の Control 群平均:", mean(as.numeric(counts_example_df["Gene1", 1:3])), "\n")
#     cat("Gene1 の Treatment 群平均:", mean(as.numeric(counts_example_df["Gene1", 4:6])), "\n")
#
#   }, error = function(e) {
#     cat("エラーが発生しました:", e$message, "\n")
#   })
# } else {
#   cat("edgeR または dplyr パッケージがインストールされていません。サンプルコードを実行できません。\n")
# }
