# 必要なパッケージをロード
library(ggplot2)
library(gridExtra)
library(stringr) # 文字列操作用パッケージ
read_mbias_data <- function(file_path) {
  text <- readLines(file_path)  
  contexts <- c("CpG", "CHG", "CHH")  
  read_groups <- c("R1", "R2")  
  
  data_list <- list()
  for (context in contexts) {
    for (read_group in read_groups) {
      start_line <- grep(paste0(context, " context \\(", read_group, "\\)"), text)
      if (length(start_line) > 0) {
        end_line <- grep("context \\(", text[(start_line+1):length(text)])  
        if (length(end_line) > 0) {
          data_text <- text[(start_line+2):(start_line + end_line[1] - 1)]  
        } else {
          data_text <- text[(start_line+2):length(text)]  
        }
        
        data <- read.table(textConnection(paste(data_text, collapse = "\n")), 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        
        # カラム名を修正
        colnames(data) <- c("position", "count_methylated", "count_unmethylated", "percentage_methylation", "coverage")
        
        data_list[[paste0(context, "_", read_group)]] <- data
      }
    }
  }
  return(data_list)
}
plot_mbias_r1 <- function(mbias_data_list, y_max_percent = 100, coverage_scale = 1000000) {
  ggplot() +
    geom_line(data = mbias_data_list$CpG_R1, aes(x = position, y = percentage_methylation, color = "CpG % Methylation")) +
    geom_line(data = mbias_data_list$CpG_R1, aes(x = position, y = coverage/coverage_scale, color = "CpG Coverage (scaled)")) +
    geom_line(data = mbias_data_list$CHG_R1, aes(x = position, y = percentage_methylation, color = "CHG % Methylation")) +
    geom_line(data = mbias_data_list$CHG_R1, aes(x = position, y = coverage/coverage_scale, color = "CHG Coverage (scaled)")) +
    geom_line(data = mbias_data_list$CHH_R1, aes(x = position, y = percentage_methylation, color = "CHH % Methylation")) +
    
    scale_y_continuous(
      name = "% Methylation",
      limits = c(0, y_max_percent),
      sec.axis = sec_axis(~.*coverage_scale, name="Coverage")
    ) +
    scale_x_continuous(name = "Position (bp)") +
    scale_color_manual(values = c(
      "CpG % Methylation" = "blue",
      "CHG % Methylation" = "forestgreen",
      "CHH % Methylation" = "purple",
      "CpG Coverage (scaled)" = "skyblue",
      "CHG Coverage (scaled)" = "lightgreen"
    )) +
    labs(title = "M-bias Plot (R1)", color = "Context") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

# R2用のグラフを作成する関数
plot_mbias_r2 <- function(mbias_data_list, y_max_percent = 100, coverage_scale = 1000000) {
  ggplot() +
    geom_line(data = mbias_data_list$CpG_R2, aes(x = position, y = percentage_methylation, color = "CpG % Methylation")) +
    geom_line(data = mbias_data_list$CpG_R2, aes(x = position, y = coverage/coverage_scale, color = "CpG Coverage (scaled)")) +
    geom_line(data = mbias_data_list$CHG_R2, aes(x = position, y = percentage_methylation, color = "CHG % Methylation")) +
    geom_line(data = mbias_data_list$CHG_R2, aes(x = position, y = coverage/coverage_scale, color = "CHG Coverage (scaled)")) +
    geom_line(data = mbias_data_list$CHH_R2, aes(x = position, y = percentage_methylation, color = "CHH % Methylation")) +
    
    scale_y_continuous(
      name = "% Methylation",
      limits = c(0, y_max_percent),
      sec.axis = sec_axis(~.*coverage_scale, name="Coverage")
    ) +
    scale_x_continuous(name = "Position (bp)") +
    scale_color_manual(values = c(
      "CpG % Methylation" = "blue",
      "CHG % Methylation" = "forestgreen",
      "CHH % Methylation" = "purple",
      "CpG Coverage (scaled)" = "skyblue",
      "CHG Coverage (scaled)" = "lightgreen"
    )) +
    labs(title = "M-bias Plot (R2)", color = "Context") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

# R1とR2を上下に配置する関数
plot_mbias_separate <- function(mbias_data_list, y_max_percent = 100, coverage_scale = 1000000) {
  p1 <- plot_mbias_r1(mbias_data_list, y_max_percent, coverage_scale)
  p2 <- plot_mbias_r2(mbias_data_list, y_max_percent, coverage_scale)
  
  # 上下に配置
  grid.arrange(p1, p2, nrow = 2)
}

# 実行関数
m_bias_plot <- function(file_path, y_max_percent = 100, coverage_scale = 1000000){
  data_list <- read_mbias_data(file_path)
  plot_mbias_separate(data_list, y_max_percent, coverage_scale)
}
