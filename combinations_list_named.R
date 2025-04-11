# 関数名を変更し、引数名を list から input_list に変更 (予約語 list との衝突回避)
combinations_list_named <- function(input_list) {
  # ユニークな要素を取得
  unique_groups <- unique(input_list)

  # ユニークな要素が2未満の場合は組み合わせを作れないので空リストを返す
  if (length(unique_groups) < 2) {
    warning("ユニークな要素が2未満のため、組み合わせを作成できません。")
    return(list())
  }

  # unique_groups から 2つの要素を選ぶ組み合わせを行列として生成
  combinations_matrix <- combn(unique_groups, 2)

  # 行列の各列をリストの要素に変換
  combination_pairs_list <- lapply(1:ncol(combinations_matrix), function(i) combinations_matrix[, i])

  # --- ここからが名前を付ける処理 ---
  # 1. 各リスト要素の名前を生成 (例: "要素1_vs_要素2")
  #    applyを使って行列の各列の要素を "_vs_" で連結
  element_names <- apply(combinations_matrix, 2, paste, collapse = "_vs_") # 区切り文字はお好みで変更可能

  # 2. 生成した名前をリストの各要素に設定
  names(combination_pairs_list) <- element_names
  # --- 名前付け処理ここまで ---

  # 名前付きのリストを返す
  return(combination_pairs_list)
}
