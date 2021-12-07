# 1. Read data
library(dplyr)
merged_gc <- read.table("data/exp.ACRG/ACRG-GC-merged-meta.txt",header = T)
subset.m <- select(merged_gc, one_of(c("Mol.subtype","NLRP3","AIM2","NLRC4","NLRP1","NLRP6","CASP1","GSDMD",
                                       "ASC","IL1B","IL17A","CD36","CXCL20","MPO","LY6G5C",
                                       "LY6G5B","PAD4","PSGL1")))
# Warning message:
# Unknown columns: `ASC`, `CXCL20`, `LY6G5B`, `PAD4`, `PSGL1`

# 2. Boxplot
library(ggplot2)
library(glue) # 根据放入的变量生成图表标题
library(ggpubr) # 将图片保存为一个pdf文件
library(tidyverse) # map2函数用于循环
library(hrbrthemes) # theme
library(Cairo) # output
library(viridis) # color
library(showtext) 
library(sysfonts)

font_add("Arial Narrow", "../../fonts/Arial Narrow.TTF")
showtext_auto()

# "NLRP3","AIM2","NLRC4","NLRP1","NLRP6","CASP1","GSDMD",
# "ASC","IL1B","IL17A","CD36","CXCL20","MPO","LY6G5C",
# "LY6G5B","PAD4","PSGL1"

x <- names(subset.m)[1]
y <- names(subset.m)[-1]

plot_list <- map2(x, y, ~ subset.m %>% 
                    ggplot(aes_string(x = .x, y = .y, fill=.x)) +
                    geom_boxplot() +
                    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
                    theme_ipsum() +
                    theme(
                      legend.position="none",
                      plot.title = element_text(size=11)
                    ) +
                    labs(title = glue('{.x} ~ {.y}')))
ggexport(plotlist = plot_list, filename = "results/ACRG-boxplot-res.pdf")                   


# 3. Stats
library(ggstatsplot)

pdf(file = "results/ACRG-NLRP3-stats.pdf")
ggbetweenstats(
  data = subset.m,
  x = Mol.subtype,
  y = NLRP3,
  pairwise.display = "all",
  p.adjust.method = "fdr"
)
dev.off()
