edgeR_modern <- function(df, group_list, reference_group = NULL){
  matrix <- as.matrix(df)
  group <- factor(group_list)
  # reference_groupが指定されていれば、それをリファレンスレベルに設定
  if (!is.null(reference_group)) {
    if (reference_group %in% levels(group)) {
      group <- relevel(group, ref = reference_group)
      cat("Reference group set to:", reference_group, "\n")
    } else {
      warning(paste("Specified reference_group '", reference_group, "' not found in group_list. Using default factor levels.", sep=""))
    }
  }
  
  cat("Factor levels for 'group' (first level is reference):\n")
  print(levels(group)) # レベルの順序を確認
  # print(group) # 実際のグループ割り当てとレベルを確認 (必要であればコメント解除)
  
  y <- DGEList(counts=matrix, group=group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  design <- model.matrix(~group)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  return(topTags(qlf, n = Inf))
}
