# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# 数据预处理 --------------------------------------------------------------
# 读取原始数据
raw_data <- read.delim("combined_metrics.txt", check.names = FALSE)

# 数据清洗转换
processed_data <- raw_data %>%
  mutate(
    across(
      c(Accuracy, Precision, Recall, F1_Score),
      ~ as.numeric(sub("%", "", .)) / 100  # 保持原始精度转换百分比
    )
  ) %>%
  pivot_longer(
    cols = -Model,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Accuracy", "Precision", "Recall", "F1_Score", "AUC")
    )
  )

# 可视化配置 --------------------------------------------------------------
# 颜色梯度配置
color_gradient <- scale_color_gradient(
  low = "#2166AC",
  high = "#B2182B",
  labels = percent_format(accuracy = 0.01),
  guide = guide_colorbar(
    order = 2,
    barwidth = unit(0.5, "cm"),
    barheight = unit(4, "cm"))
)

# 尺寸标度配置
size_scale <- scale_size_continuous(
  range = c(10, 22),
  breaks = seq(0.7, 1, 0.05),
  labels = percent_format(accuracy = 0.01),
  guide = guide_legend(
    order = 1,
    override.aes = list(color = "#B2182B"))
)

# 构建图形对象 ------------------------------------------------------------
performance_plot <- ggplot(
  processed_data,
  aes(x = Model, y = Metric)
) +
  # 几何对象层
  geom_point(
    aes(size = Value, color = Value, fill = Value),  # 添加 fill 映射
    alpha = 0.8,
    shape = 21,
    stroke = 0.5
  ) +
  geom_text(
    aes(label = ifelse(
      Metric == "AUC",
      sprintf("%.2f", Value),
      sprintf("%.2f%%", Value * 100)
    )),
    color = "black",
    size = 3.2,
    fontface = "bold"
  ) +
  # 标度配置
  scale_fill_gradient(  # 添加 fill 的颜色梯度
    low = "#2166AC",
    high = "#B2182B",
    labels = percent_format(accuracy = 0.01),
    guide = guide_colorbar(
      order = 2,
      barwidth = unit(0.5, "cm"),
      barheight = unit(4, "cm")
    )
  ) +
  scale_color_gradient(  # 确保 color 也有正确的梯度
    low = "#2166AC",
    high = "#B2182B",
    labels = percent_format(accuracy = 0.01),
    guide = guide_colorbar(
      order = 2,
      barwidth = unit(0.5, "cm"),
      barheight = unit(4, "cm")
    )
  ) +
  size_scale +
  
  # 坐标系统
  coord_fixed(ratio = 1.2) +  # 关键参数调整
  
  # 主题配置
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold",
      margin = margin(b = 10)
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 11,
      margin = margin(r = 5)
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.8
    ),
    legend.spacing = unit(0.5, "cm"),
    legend.box.margin = margin(l = 10),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    aspect.ratio = 1  # 强制1:1宽高比
  ) +
  
  # 标签配置
  labs(
    title = "Model Performance Metrics Comparison",
    size = "Metric Value",
    color = "Metric Value",
    fill = "Metric Value"  # 添加 fill 的标签
  )

# 图形输出 ----------------------------------------------------------------
ggsave(
  filename = "model_performance_square.pdf",
  plot = performance_plot,
  width = 8.27,  # A4纸宽度（英寸）
  height = 8.27, # 保持正方形比例
  units = "in",
  dpi = 600
)