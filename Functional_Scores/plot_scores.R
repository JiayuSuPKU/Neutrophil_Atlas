## This file contains vidualization R codes for several functional scores

setwd("path")
load("./Integrated.RData")

library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)

# Secretory vesicle scores for all cells
p0 = ggplot(cell.info, aes(x=cluster_condition, y=secretory_vesicle.score, fill=condition)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(title="Secretory Vesicle") + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                              brewer.pal(9, "Set1")[5])) + 
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none",
    plot.title = element_text(size=7, hjust = 0.5)
  )

## Visualizations of functional scores for G1~G5 cells
cell.info.G1_5 = subset(cell.info, cell.info$cluster_union != "G0")

# Primary granule scores
p1 = ggplot(cell.info.G1_5, aes(x=cluster_union, y=Pri.graules.score, fill=condition)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                               brewer.pal(9, "Set1")[5])) + 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=7,angle = 90,colour = "black"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size =7),
        plot.title = element_blank(),
        legend.title = element_text(colour="black", size=7),
        legend.text=element_text(size=7),
        legend.background = element_blank())+
  labs(title = "Azurophil", x = "Cluster", y = "Azurophil score") + 
  labs(fill="Condition") +
  stat_compare_means(aes(group = cell.info.G1_5$condition),
                     label = "p.signif",
                     size = 2,
                     method = "t.test")

# Second granule scores
p2 = ggplot(cell.info.G1_5, aes(x=cluster_union, y=second_granules.score, fill=condition)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                               brewer.pal(9, "Set1")[5])) + 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=7,angle = 90,colour = "black"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size =7),
        plot.title = element_blank(),
        legend.title = element_text(colour="black", size=7),
        legend.text=element_text(size=7),
        legend.background = element_blank())+
  labs(title = "Specific", x = "Cluster", y = "Specific score") + 
  labs(fill="Condition") +
  stat_compare_means(aes(group = cell.info.G1_5$condition),
                     label = "p.signif",
                     size = 2,
                     method = "t.test")

# Tertiary granule scores
p3 = ggplot(cell.info.G1_5, aes(x=cluster_union, y=tertiary_granules.score, fill=condition)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                               brewer.pal(9, "Set1")[5])) + 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=7,angle = 90,colour = "black"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size =7),
        plot.title = element_blank(),
        legend.title = element_text(colour="black", size=7),
        legend.text=element_text(size=7),
        legend.background = element_blank())+
  labs(title = "Gelatinase", x = "Cluster", y = "Gelatinase score") + 
  labs(fill="Condition") +
  stat_compare_means(aes(group = cell.info.G1_5$condition),
                     label = "p.signif",
                     size = 2,
                     method = "t.test")

# Secretory vesicle scores
p4 = ggplot(cell.info.G1_5, aes(x=cluster_union, y=secretory_vesicle.score, fill=condition)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                               brewer.pal(9, "Set1")[5])) + 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=7,angle = 90,colour = "black"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size =7),
        plot.title = element_blank(),
        legend.title = element_text(colour="black", size=7),
        legend.text=element_text(size=7),
        legend.background = element_blank())+
  labs(title = "Secretory", x = "Cluster", y = "Secretory score") + 
  labs(fill="Condition") +
  stat_compare_means(aes(group = cell.info.G1_5$condition),
                     label = "p.signif",
                     size = 2,
                     method = "t.test")

# NADPH scores
p5 = ggplot(cell.info.G1_5, aes(x=cluster_union, y=NADPH.score, fill=condition)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                               brewer.pal(9, "Set1")[5])) + 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=7,angle = 90,colour = "black"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size =7),
        plot.title = element_blank(),
        legend.title = element_text(colour="black", size=7),
        legend.text=element_text(size=7),
        legend.background = element_blank())+
  labs(title = "NADPH", x = "Cluster", y = "NADPH complex score") + 
  labs(fill="Condition") +
  stat_compare_means(aes(group = cell.info.G1_5$condition),
                     label = "p.signif",
                     size = 2,
                     method = "t.test")

# Chemotaxis scores
p6 = ggplot(cell.info.G1_5, aes(x=cluster_union, y=chemotaxis.score, fill=condition)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                               brewer.pal(9, "Set1")[5])) + 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=7,angle = 90,colour = "black"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size =7),
        plot.title = element_blank(),
        legend.title = element_text(colour="black", size=7),
        legend.text=element_text(size=7),
        legend.background = element_blank())+
  labs(title = "Chemotaxis", x = "Cluster", y = "Chemotaxis score") + 
  labs(fill="Condition") +
  stat_compare_means(aes(group = cell.info.G1_5$condition),
                     label = "p.signif",
                     size = 2,
                     method = "t.test")

# Phagocytosis scores
p7 = ggplot(cell.info.G1_5, aes(x=cluster_union, y=phagocytosis.score, fill=condition)) +
  geom_boxplot(notch=TRUE,outlier.shape=NA) + 
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[9], 
                               brewer.pal(9, "Set1")[5])) + 
  theme_classic()+
  theme(legend.position="none",
        axis.text=element_text(size=7,angle = 90,colour = "black"),
        axis.title.x =element_blank(),
        axis.title.y = element_text(size =7),
        plot.title = element_blank(),
        legend.title = element_text(colour="black", size=7),
        legend.text=element_text(size=7),
        legend.background = element_blank())+
  labs(title = "Phagocytosis", x = "Cluster", y = "Phagocytosis score") + 
  labs(fill="Condition") +
  stat_compare_means(aes(group = cell.info.G1_5$condition),
                     label = "p.signif",
                     size = 2,
                     method = "t.test")

## Save figures
pdf("path/name.pdf", width = 6.5,height = 2.7)
cowplot::plot_grid(p1,p2,p3,p4,p5,p,p6,p7,nrow = 2, align = "hv")
dev.off()

## Pie plot for cell proportions
# Pie chart
p1 = ggplot(data.frame(table(as.character(cell.info[cell.info$condition == "Control",]$tissue))), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  labs(fill="") +
  scale_fill_manual(values=brewer.pal(12,"Paired")[c(1,2,4,8)]) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    legend.text = element_text(size = 14, color = "black"),
    legend.position = "left",
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
)

p2 = ggplot(data.frame(table(as.character(cell.info[cell.info$condition == "Ecoli",]$tissue))), aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  labs(fill="") +
  scale_fill_manual(values=brewer.pal(12,"Paired")[c(1,2,10,4,12,8)]) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    legend.text = element_text(size = 14, color = "black"),
    legend.position = "left",
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
)
pdf("path/name.pdf", width = 65*4/72, height = 142*4/72)
cowplot::plot_grid(p1,p2,nrow = 2, align = "hv")
dev.off()
