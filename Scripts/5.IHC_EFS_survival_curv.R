setwd('path/to/data')
#
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)

EFS<-read.csv('KRAS病人EFS_test_50.csv')

proteins<-read.csv('PDE4D_CD8.csv')

proteins <- proteins %>%
  # 用 separate 分开 Percentage 和 Modifier
  separate(PDE4D, into = c("Percentage", "Modifier"), sep = "%", fill = "right") %>%
  # 将 Percentage 转换为数值
  mutate(
    Percentage = as.numeric(Percentage),
    # 计算 Modifier 的个数，默认空为 0
    Modifier_count = if_else(is.na(Modifier), 0, str_length(Modifier)),
    # 计算新列 Result
    PDE4D = Percentage * Modifier_count
  )
proteins <- proteins %>%
  mutate(CD8 = as.numeric(str_remove(CD8, "%")))

cor_test <- cor.test(proteins$PDE4D, proteins$CD8, method = "pearson")

library(ggpubr)
library(ggrepel)
########Fig6b#######
ggscatter(
  proteins,
  x = "PDE4D",
  y = "CD8",
  add = "reg.line",  # 添加回归线
  conf.int = T,   # 显示置信区间
  cor.coef = TRUE,   # 显示相关系数
  cor.method = "pearson",  # 相关性方法
  xlab = "PDE4D protein H-score",
  ylab = "CD8 T cells percentage (%)",
  cor.coef.coord = c(150, 31) # 指定R和P值的坐标
)#+geom_text_repel(aes(label = ID), size = 3, data = proteins,max.overlaps = 30) 


######
EFS$ID <- gsub("cg$", "", EFS$ID)
data<-merge(EFS,proteins,by='ID')

cutoff<-50
group <- ifelse(data$PDE4D> cutoff, "high", "low")  # 分组
data$group<-as.factor(group)
data<-data[order(data$PDE4D,decreasing = T),]
write.csv(data,file='integrated_data_49samples.csv',row.names = F)

data_summary <- data %>%
  group_by(group) %>%
  summarise(
    mean_CD8 = mean(CD8),
    se_CD8 = sd(CD8) / sqrt(n()),  # 标准误差
    .groups = 'drop'
  )


###画图
# Fig6c#####
ggplot(data, aes(x = group, y = CD8, fill = group)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 1.5, shape = 16,
              alpha = 0.8, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.size = 30, label.x = 1.5, label.y = 32) +
  scale_x_discrete(labels = c("high" = "PDE4D-high", "low" = "PDE4D-low")) +
  labs(x = "", y = "CD8+ T cell percentage") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Wilcoxon rank-sum test
wilcox_CD8 <- wilcox.test(CD8 ~ group, data = data, exact = FALSE)
print(wilcox_CD8)

data$PDE4D_value<-data$PDE4D
data$PDE4D<-data$group
surv<-Surv(time = data$EFS,event = data$Event)
fit<-survfit(surv~PDE4D,data=data)
#########Fig6d####
ggsurvplot(fit,data=data,pval = T,pval.method = F,risk.table = TRUE,risk.table.col = "strata",risk.table.height = 0.25,xlab="Months",break.time.by = 5,ylab='Event-free survival')
