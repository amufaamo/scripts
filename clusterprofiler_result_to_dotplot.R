# clusterprofiler result to dotplot

ggplot(as.data.frame(data), aes(x = RichFactor, y = reorder(Description, RichFactor))) +
    geom_point(aes(size = Count, color = qvalue)) +  # ドットの描画
    scale_color_gradient(low = "blue", high = "red") + # Qvalueに応じた色 (低い方が赤)
    scale_size_continuous(range = c(6, 8)) +  # サイズの範囲を3から10に設定
    labs(
        title = "Pathways enriched by genes belonging to the AD disease module", #必要に応じてタイトルを変更してください。
        x = "Rich Factor",
        y = NULL,  # y軸ラベルを削除
        size = "Gene Number",  # 凡例のタイトル
        color = "Qvalue"  # 凡例のタイトル
    ) 
