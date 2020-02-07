setwd("path")
load("./Data.RData")

library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ComplexHeatmap)

# Color annotation (for ComplexHeatmap)
ann.color <- list(Cluster = c(G0=pal_npg()(10)[1], 
                              G1=pal_npg()(10)[2], 
                              G2=pal_npg()(10)[3], 
                              G3a=pal_npg()(10)[4],
                              G4=pal_npg()(10)[5], 
                              G5a=pal_npg()(10)[6], 
                              G5b=pal_npg()(10)[7], 
                              G5c=pal_npg()(10)[8], 
                              Mono=pal_npg()(10)[9],
                              G3b=pal_npg()(10)[10],
                              DC="gray9",
                              B="#7A57D1",
                              HSPC="#0000A1",
                              "T"="#1F6ED4"))

## ---------------------------------------------------
## Part 1: Aging Score
## ---------------------------------------------------

# Define age-related genes and their relative weights. 
# Genes were summarized from previous literatures.
aging.genes = c("Sell","Itgam","Itga4","Cxcr4","Cxcr2",
                "Cd47","Cd24a","Tlr4","Icam1","Itgax")
weights = c(-1,1,1,1,-1,-1,1,1,1,1)/10
# Remove undetected genes
aging.genes.refined = c("Cxcr4", "Sell", "Itga4", "Itgam", 
                        "Cxcr2", "Cd24a", "Icam1", "Cd47")
weights.refined = c(1,-1,1,1,-1,1,1,-1)/10

# Centered-scale expression matrix 
genes.used = rowSums(ref.merge.neu_G5@raw.data>0)>10
z_matrix = t(scale(t(ref.merge.neu_G5@data[genes.used,]),
                   center = T, scale = T))

# Calculate aging score as the weighted average of age-related gene expressions
aging.score = t(z_matrix[aging.genes.refined,]) %*% weights.refined 
# Define aged subgroups by Gaussian mixture model
require(mixtools)
out = normalmixEM(aging.score, k=2)
plot(out,2) # check the order
aging.group = apply(out$posterior, 1, which.max)

# Optional: Multi-variable GMM
#cov = cor(as.matrix(t(ref.merge.neu_G5@data[aging.genes,])))
#out = mvnormalmixEM(t(ref.merge.neu_G5@data[aging.genes,]), sigma = list(cov, cov))

# Optional: Find de novo age-related genes 
#corr = apply(ref.merge.neu_G5@data[genes.used,], 1, function(x){cor(x,aging.score)})
#corr = corr[!is.na(corr)]
#corr = cor(as.matrix(t(ref.merge.neu_G5@data[genes.used,])), 
#           as.matrix(t(ref.merge.neu_G5@data[aging.genes,])))


# Set data ready for visualization
df.aging = data.frame(aging.score, ref.merge.neu_G5@meta.data)
# Check the aged group's composition
table(aging.group, df.aging$cluster)
tapply(df.aging$aging.score, df.aging$cluster, summary)
# Check the proportion
t(t(table(aging.group, df.aging$cluster))/colSums(
  table(aging.group,df.aging$cluster)))
ggplot(df.aging, aes(x=cluster,y=aging.score, fill=cluster)) + 
  geom_violin() + theme_bw()

# Visualization
require(ggpubr)
p1 = ggplot(df.aging, aes(x=cluster,y=aging.score, fill=cluster)) + 
  geom_violin() + geom_boxplot(width=0.2,outlier.shape = NA) + 
  ylab("Aging Score")+
  stat_compare_means(comparisons = list(c("G5a", "G5b"), c("G5a", "G5c")), 
                     label = "p.signif", label.y = c(1.7, 2),
                     size=3, method = "t.test")+
  scale_fill_manual(values = ann.color$Cluster)+
  theme_minimal() + 
  theme(
    axis.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

df.aging.group = data.frame(Group=c("G5a","G5b","G5c"), 
                            prop = c(0.04654443,0.06561361,0.15179114))
p2 = ggplot(df.aging.group, aes(x=Group, y=prop, fill=Group)) +
  geom_bar(width = 0.75, stat = "identity", color="black") +
  ylab("Proportion of aged cells")+
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )
pdf("Aging_score_prop.pdf", width = 2.4, height = 6)
cowplot::plot_grid(p1,p2,nrow = 2, align = "hv")
dev.off()

## ---------------------------------------------------
## Part 2: Apoptosis Score
## ---------------------------------------------------

# Apoptosis score is pre-calculated in a similar way as scores in Part3
# Load pre-calculated apoptosis score
df.apoptosis = data.frame(
  score = ref.merge.neu_all@meta.data[ref.merge.neu_all@meta.data$cluster 
                                      %in% c("G3a","G4","G5a","G5b","G5c"),]$apoptosis.score,
  cluster = ref.merge.neu_all@meta.data[ref.merge.neu_all@meta.data$cluster 
                                        %in% c("G3a","G4","G5a","G5b","G5c"),]$cluster)

# Apply Gaussian mixture model
require(mixtools)
out = normalmixEM(df.apoptosis$score, k=2, maxit = 5000)
plot(out,2)

# Check compositions
df.apoptosis$apop_group = apply(out$posterior, 1, which.max)
tb = table(df.apoptosis$apop_group, df.apoptosis$cluster)
t(t(tb)/colSums(tb))

# Visualizations
p1 = ggplot(df.apoptosis, aes(x=cluster, y=score, fill=cluster)) +
  geom_violin(scale="width") +
  #scale_x_discrete(limits = rev(levels(factor(df.apoptosis$cluster)))) + 
  #coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Apoptosis score") + 
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none",
    plot.title = element_text(size=7, hjust = 0.5)
  )


df.apop.group = data.frame(Group=c("G3a","G4","G5a","G5b","G5c"),
                           prop = c(0.02764423,0.03647251,0.08152327,
                                    0.16767922,0.30297511))
p2 = ggplot(df.apop.group, aes(x=Group, y=prop, fill=Group)) +
  geom_bar(width = 0.75, stat = "identity", color="black") +
  labs(y="Proportion of apoptotic cells")+
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7,color = "Black"),
    plot.title = element_text(size=7,color = "Black", hjust=0.5),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

pdf("E:/Scores/Apop.pdf", width=189/72,height=136/72)
cowplot::plot_grid(p1,p2,align = "h")
dev.off()

## ---------------------------------------------------
## Part 3: Pyroptosis and necroptic Scores
## ---------------------------------------------------

## Calculate pyroptosis scores
pyroptosis.gene = read.csv("pyroptosis.csv",header = T)
colnames(pyroptosis.gene) = "genes"
ref.merge.neu_all@meta.data$pyroptosis.score = 
  Matrix::colMeans(ref.merge.neu_all@data[pyroptosis.gene$genes,])


## Calculate necroptic scores
necroptic.genes = read.csv("necroptic_process.csv", header = T)
colnames(necroptic.genes) = "genes"
ref.merge.neu_all@meta.data$necroptic.score = 
  colMeans(ref.merge.neu_all@data[necroptic.genes$genes,])

df.G3_G5 = data.frame(
  pyroptosis = ref.merge.neu_all@meta.data[ref.merge.neu_all@meta.data$cluster 
                                           %in% c("G3a","G4","G5a","G5b","G5c"),]$pyroptosis.score,
  necroptic = ref.merge.neu_all@meta.data[ref.merge.neu_all@meta.data$cluster 
                                          %in% c("G3a","G4","G5a","G5b","G5c"),]$necroptic.score,
  cluster = ref.merge.neu_all@meta.data[ref.merge.neu_all@meta.data$cluster 
                                        %in% c("G3a","G4","G5a","G5b","G5c"),]$cluster)

# Apply GMM
require(mixtools)
out = normalmixEM(df.G3_G5[df.G3_G5$pyroptosis>0,]$pyroptosis, k=2)
plot(out,2)

# Set data ready for visualization
df.G3_G5$pyro_group = 1
df.G3_G5[df.G3_G5$pyroptosis>0,]$pyro_group = apply(out$posterior, 1, which.max)
tb = table(df.G3_G5$pyro_group, df.G3_G5$cluster)
t(t(tb)/colSums(tb))
df.pyro.group = data.frame(Group=c("G3a","G4","G5a","G5b","G5c"),
                           prop = c(0.03245192,0.05498095,0.06854725,
                                    0.10206561,0.08378871))

df.G3_G5$nec_group = apply(out$posterior, 1, which.max)
tb = table(df.G3_G5$nec_group, df.G3_G5$cluster)
t(t(tb)/colSums(tb))
df.nec.group = data.frame(Group=c("G3a","G4","G5a","G5b","G5c"),
                          prop = c(0.1718750,0.2090365,0.2578279,
                                   0.2831106,0.2659381))

ggplot(df.pyro.group, aes(x=Group, y=prop, fill=Group)) +
  geom_bar(width = 0.75, stat = "identity", color="black") +
  ylab("Proportion of pyroptotic cells")+
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

ggplot(df.nec.group, aes(x=Group, y=prop, fill=Group)) +
  geom_bar(width = 0.75, stat = "identity", color="black") +
  ylab("Proportion of necroptic cells")+
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

## ---------------------------------------------------
## Part 4: Granule Scores
## ---------------------------------------------------

# Visualization for granule gene expressions
genes.used = c("Mpo","Elane","Ctsg","Prtn3","Prss57","Ctsc",
               "Camp","Ltf","Cybb","Cyba","Lcn2","B2m",
               "Mmp8","Mmp9","Hp","Slpi","Itgam",
               "Cd63","Fcgr3","Cd177","Mmp25","Fpr1","Scamp1","Vamp2","Stxbp4")
data.granule = as.matrix(ref.merge.neu_all@data[genes.used,])
hm.1 = Heatmap(data.granule, 
               name = "Exp",
               cluster_rows = F, cluster_columns = T,
               show_row_names = T,
               show_column_names = F,
               show_column_dend = F,
               col = circlize::colorRamp2(seq(-2,2,0.5), 
                                          rev(brewer.pal(n=9, "RdBu"))) ,
               column_names_gp = gpar(fontsize=7),
               row_names_gp = gpar(fontsize=7),
               top_annotation = HeatmapAnnotation(cluster=ref.merge.neu_all@meta.data$cluster,
                                                  col=ann.color, which="column",
                                                  show_annotation_name = F),
               #width = unit(2.5, "cm"), height = unit(5,"cm"),
               heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                           labels_gp = gpar(fontsize=7),
                                           grid_width = unit(2.5, "mm")))


df.granule = data.frame(ref.merge.neu_all@dr$umap@cell.embeddings,
                        ref.merge.neu_all@meta.data,
                        t(ref.merge.neu_all@data[c("Elane","Prtn3","Ltf","Camp","Mmp8","Mmp9"),]))
p1 = ggplot(df.granule, aes(x=UMAP1, y=UMAP2, color=Elane)) + 
  geom_point(alpha=0.8, size=0.1) + 
  scale_color_material("indigo") +
  ggtitle("Elane") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid=element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p2 = ggplot(df.granule, aes(x=UMAP1, y=UMAP2, color=Prtn3)) + 
  geom_point(alpha=0.8, size=0.1) + 
  scale_color_material("indigo") +
  ggtitle("Prtn3") + 
  theme_minimal() + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid=element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p3 = ggplot(df.granule, aes(x=UMAP1, y=UMAP2, color=Ltf)) + 
  geom_point(alpha=0.8, size=0.1) + 
  scale_color_material("indigo") +
  ggtitle("Ltf") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid=element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p4 = ggplot(df.granule, aes(x=UMAP1, y=UMAP2, color=Camp)) + 
  geom_point(alpha=0.8, size=0.1) + 
  scale_color_material("indigo") +
  ggtitle("Camp") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid=element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p5 = ggplot(df.granule, aes(x=UMAP1, y=UMAP2, color=Mmp8)) + 
  geom_point(alpha=0.8, size=0.1) + 
  scale_color_material("indigo") +
  ggtitle("Mmp8") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid=element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p6 = ggplot(df.granule, aes(x=UMAP1, y=UMAP2, color=Mmp9)) + 
  geom_point(alpha=0.8, size=0.1) + 
  scale_color_material("indigo") +
  ggtitle("Mmp9") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid=element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
pdf("F2B.pdf", width=2.44, height=4.25)
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,nrow = 3, align = "hv")
dev.off()

p1 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=Pri.graules.score, fill=cluster)) +
  geom_violin(scale="width") +
  scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Azurophil") + 
  scale_fill_manual(values=ann.color$Cluster) +
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
p2 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=second_granules.score, fill=cluster)) +
  geom_violin(scale="width") +
  scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Specific") + 
  scale_fill_manual(values=ann.color$Cluster) +
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
p3 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=tertiary_granules.score, fill=cluster)) +
  geom_violin(scale="width") +
  scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Gelatinase") + 
  scale_fill_manual(values=ann.color$Cluster) +
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
p4 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=Secretory_Vesicle.score, fill=cluster)) +
  geom_violin(scale="width") +
  scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Secretory Vesicle") + 
  scale_fill_manual(values=ann.color$Cluster) +
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
pdf("E:/Figures/Figure 2/F2C.pdf", width = 3, height = 2.67)
cowplot::plot_grid(p1,p2,p3,p4,nrow = 2,align = "hv")
dev.off()

## ---------------------------------------------------
## Part 5: Other figures in the paper
## ---------------------------------------------------

### Fig.S3
p1 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=nUMI, fill=cluster)) +
  geom_violin(scale="width") +
  #scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  #coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Number of UMI") + 
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none",
    plot.title = element_text(size=7, hjust = 0.5)
  )
p2 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=nGene, fill=cluster)) +
  geom_violin(scale="width") +
  #scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  #coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Number of genes") + 
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none",
    plot.title = element_text(size=7, hjust = 0.5)
  )
p3 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=nUMI.nGene.ratio, fill=cluster)) +
  geom_violin(scale="width") +
  #scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  #coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="UMI per gene") + 
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none",
    plot.title = element_text(size=7, hjust = 0.5)
  )
p4 = ggplot(ref.merge.neu_all@meta.data, aes(x=cluster, y=percent.mito, fill=cluster)) +
  geom_violin(scale="width") +
  #scale_x_discrete(limits = rev(levels(factor(ref.merge.neu_all@meta.data$cluster)))) + 
  #coord_flip() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  labs(title="Mitochondrial count percetage") + 
  scale_fill_manual(values=ann.color$Cluster) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none",
    plot.title = element_text(size=7, hjust = 0.5)
  )
pdf("E:/Figures/FS3B.pdf", width=240/72, height = 160/72)
cowplot::plot_grid(p1,p2,p3,p4, nrow = 2, align = "hv")
dev.off()
