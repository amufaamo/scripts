heatmap_rnaseq <- function(count, gene_list){
    # 1. DGEListオブジェクトの作成
    y <- DGEList(counts=count, group=names(count))
    # 2. フィルタリング (発現量の低い遺伝子を除去)
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE] # フィルタリングされた遺伝子でDGEListオブジェクトを更新
    # 3. 正規化 (TMM正規化)
    y <- calcNormFactors(y)
    logcpm <- cpm(y, log = TRUE)
    logcpm <- logcpm[gene_list,]
    pheatmap::pheatmap(as.matrix(logcpm), scale = 'row') # pheatmapで図を作成
}
