volcano_plot_from_edgeR_table <- function(df){

    # FDRを-log10変換した列を追加
    df$negLog10FDR <- -log10(df$FDR)

    # ボルケーノプロットの作成 (PValueを使用する場合)
    volcano_plot <- ggplot(data = df, aes(x = logFC, y = negLog10FDR)) +
        geom_point(alpha = 0.4) + # 点の透明度を設定
        theme_classic() + # シンプルなテーマ
        # 有意差の閾値ラインを追加（例：PValue < 0.05）
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        # logFCの閾値ラインを追加（例：|logFC| > 1）
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")
    return(volcano_plot)
}
