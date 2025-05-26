edgeR_modern <- function(df, group_list, target_group, reference_group) {
  cat("Starting edgeR analysis for comparison:", target_group, "vs", reference_group, "\n")
  
  matrix_data <- as.matrix(df)
  group <- factor(group_list)
  
  # DGEListオブジェクトの作成
  y <- DGEList(counts = matrix_data, group = group)
  
  # 低発現遺伝子のフィルタリング [cite: 142]
  # filterByExprは、実験デザインに基づいて賢明なフィルタリングを行う [cite: 150]
  # フィルタリングは、DE解析に関与する因子に基づいて行うべき [cite: 156]
  keep <- filterByExpr(y, group = group) 
  y <- y[keep, , keep.lib.sizes = FALSE]
  cat(sum(keep), "genes kept after filtering.\n")
  
  # ライブラリサイズの正規化 (TMM) [cite: 175]
  y <- normLibSizes(y)
  cat("Library sizes normalized using TMM.\n")
  
  # デザイン行列の作成 (インターセプトなし) [cite: 487]
  # 各グループに対する係数を作成する
  design <- model.matrix(~0 + group, data = y$samples)
  colnames(design) <- levels(group) # makeContrastsが正しく機能するように列名をレベル名に設定 [cite: 487]
  cat("Design matrix created with columns:\n")
  print(colnames(design))
  
  # コントラストの作成 [cite: 490, 491]
  # 例: "Trt1-Ctrl"
  contrast_str <- paste(target_group, "-", reference_group)
  cat("Creating contrast for the comparison:", contrast_str, "\n")
  
  my_contrast <- makeContrasts(contrasts = contrast_str, levels = design)
  cat("Contrast matrix:\n")
  print(my_contrast)
  
  # 分散の推定とGLMのフィッティング
  # glmQLFitはロバストな推定を推奨 [cite: 812]
  # estimateDispの実行はオプションであり、glmQLFit内で分散が推定される [cite: 331, 804]
  # y <- estimateDisp(y, design, robust = TRUE) # 必要に応じて実行
  # plotBCV(y) # BCVプロットの表示
  
  fit <- glmQLFit(y, design, robust = TRUE) 
  cat("QL dispersions estimated and GLM fitted.\n")
  # plotQLDisp(fit) # QL分散プロットの表示
  
  # 指定されたコントラストに対する準尤度F検定の実行 [cite: 489]
  qlf <- glmQLFTest(fit, contrast = my_contrast)
  cat("QL F-test performed for the contrast.\n")
  
  # 結果の抽出 (上位の変動遺伝子)
  cat("Extracting top tags...\n")
  results <- topTags(qlf, n = Inf)
  
  cat("edgeR analysis complete for:", target_group, "vs", reference_group, "\n")
  return(results)
}
