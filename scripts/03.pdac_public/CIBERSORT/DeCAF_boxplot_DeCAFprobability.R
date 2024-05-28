# load libraries
library(gtable)
library(gridExtra)
library(grid)
library(gplots)
library(openxlsx)
library(tidyverse)
library(reshape2)
library(ggtext)

gc()

data <- read_excel("/work/users/c/h/changfei/CAPER_Paper/00_DeCAF_paper_summary/Fig3/Fig3_I/data/Supplementary Table S5 - clinical and molecular data - public.xlsx")

data <- as.data.frame(data)
data$Dataset <- factor(data$Dataset,levels = c("TCGA_PAAD",
                                               "CPTAC",
                                               "Moffitt_GEO_array",
                                               "PACA_AU_seq",
                                               "PACA_AU_array",
                                               "Linehan",
                                               "Puleo_array",
                                               "Olive",
                                               "Dijk",
                                               "Grunwald",
                                               "Hayashi") )


# Plot
point_size <- 0.5
p1 <- ggplot(data = data, aes(x = PurIST, y = DeCAF_prob)) +
  geom_boxplot(aes(color = PurIST), fill = "transparent", lwd = 0.8) +
  geom_point(aes(color = DeCAF), position = position_jitter(width = 0.2, height = 0), size = point_size, alpha = 0.5) +
  scale_color_manual(values = c("Basal-like" = "orange", "Classical" = "blue","permCAF" = "#FE46A5", "restCAF" = "#008080")) +
  # scale_color_manual(values = c("permCAF" = "#FE46A5", "restCAF" = "#008080")) +
    geom_signif(comparisons = list(c("Classical", "Basal-like")), size = 0.2, textsize = 8, margin_top = 0.075, 
              test = wilcox.test, color = "black", map_signif_level = FALSE, step_increase = 0.1) +
  
  theme(
    legend.text = element_text(size = 15),
    # axis.text = element_blank(),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    # strip.text = element_markdown(angle = 50, size = 25, hjust = 0.8, vjust = 0.8),
    # strip.background = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_line(color = "black"),
    # panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.title.y = element_text(size = 25),
    # axis.ticks.x = element_blank()
  ) +
  
  labs(
    title = "",
    y = "DeCAF probability",
    x = "",
    fill = NULL
  ) +
  
  scale_y_continuous(
    breaks = c(0.00, 0.5, 1.0),
    labels = c("0", "0.5", "1.0"),
    limits = c(0, 1.2)
  ) +
  guides(
    color = guide_legend(title = NULL),
    fill = guide_legend(title = NULL)
  )
p1

pdf("/work/users/c/h/changfei/CAPER_Paper/00_DeCAF_paper_summary/Fig3/Fig3_I/figures/DeCAFpaper_boxplot_datasets_calls_all.pdf",width = 6,height = 6)
p1
dev.off()


