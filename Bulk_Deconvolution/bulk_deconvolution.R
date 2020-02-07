setwd("path")
load("./BulkClassifier.RData")

library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

## ---------------------------------------------------
## Part 1: Load and preprocess bulk RNA-seq data
## ---------------------------------------------------
library(edgeR)

# Load bulk RNA-seq data
files = dir(path="./Bulk Data/", pattern="*\\.tab$")
files = paste("./Bulk Data/", files, sep = "")
dge = readDGE(files, columns=c(1,3), header=FALSE)

# Load 10x gene annotation
feature.anno <- read.table("path/features.tsv.gz")
rownames(feature.anno) = feature.anno$V1

# Keep shared genes only
gene.intersect = intersect(rownames(dge), feature.anno$V1)
dge = dge[gene.intersect,]
rownames(dge) = feature.anno[rownames(dge),]$V2
  
targets<-read.table('./Bulk Data/targets.2.txt', header=F)
colnames(targets)<-c("Sample", "Replicate", "Genotype", "Celltype", "Location", "Stage")
colnames(dge)=targets$Sample

## ---------------------------------------------------
## Part 2: Map bulk samples to sc clusters
## ---------------------------------------------------

# Quantile normalization function
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

# Correlation matrix between bulk samples after normalization
bulk.norm = FQnorm(dge$counts)
bulk.corr = cor(bulk.norm)
bulk.corr = cor(dge$counts)
Heatmap(bulk.corr)

## Using pca transfering to map each bulk sample back to single cell clusters
bulk.pca = prcomp(t(bulk.norm))
bulk.pca.transfer = t(bulk.norm[ref.merge.neu_all@var.genes,]) %*% ref.merge.neu_all@dr$pca@gene.loadings
df.bulk = data.frame(bulk.pca$x[,1:2], celltype=targets$Celltype)
ggplot(df.bulk, aes(x=PC1, y=PC2, color=celltype)) + geom_point()

# Combine bulk replicates
bulk.norm.merge = cbind(rowMeans(bulk.norm[,1:3]),rowMeans(bulk.norm[,4:6]),rowMeans(bulk.norm[,7:9]),
                        rowMeans(bulk.norm[,10:12]),rowMeans(bulk.norm[,13:15]))
colnames(bulk.norm.merge) = c("Neu", "MB", "MC", "MM", "PM")
bulk.norm.merge = bulk.norm.merge[,c(2,5,3,4,1)]

## ---------------------------------------------------
## Part 3: Bulk deconvolution by correlation and 
##   by nnls regression
## ---------------------------------------------------

# Load sc signatures genes 
DEGenes = read.csv("bm_deg.csv", row.names = 1)
# nnls deconvolution is robust to the choice of num_markers
num_markers = 20
top <- DEGenes %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
genes.used = as.character(top$gene)[top$gene %in% rownames(bulk.norm.merge)[rowSums(bulk.norm.merge)>0]]
# Add Cd34 to the previous list
genes.used = c("Cd34", genes.used)

# # Other marker genes lists used in the paper
# previous.genes = c("Itgb2l","Wfdc21","Anxa1","Cd177","Lcn2",
#                     "Ly6g","AA467197","Ifitm6","Retnlg","Mmp8",
#                     "H2afz","Birc5","Fcnb","H2afx",
#                     "Mpo","Ctsg","Elane","Prtn3","Ms4a3",
#                     "Cd34","Rpl12","Rpl32","Rplp1","Rpl26","Eef1a1")
# previous.genes.index = match(previous.genes, genes.used)

# highlight.genes = top %>% group_by(cluster) %>% top_n(n=5, wt = avg_logFC)
# highlight.genes = intersect(genes.used, highlight.genes$gene)
# highlight.genes.index = match(highlight.genes, genes.used)

# Load pooled sc expression profiles
scdata.merge = read.csv("mouse_ave_exp_cluster_all.csv")
scdata.merge = scdata.merge[!duplicated(scdata.merge$Gene.symbol),1:6]
rownames(scdata.merge) = scdata.merge$Gene.symbol
scdata.merge = scdata.merge[,2:6]
# Log-normalization
scdata.merge.norm = log1p(scdata.merge)

# Visualization
hm1 = Heatmap(bulk.norm.merge.scaled, 
        name = "Exp",
        cluster_rows = T, cluster_columns = F,
        col = rev(brewer.pal(n=9, "RdBu")),
        column_names_gp = gpar(fontsize=7),
        show_row_names = F,
        width = unit(2.5, "cm"), height = unit(10,"cm"),
        right_annotation = rowAnnotation(
          hightlight = anno_mark(at = previous.genes.index, labels = previous.genes,
                                 labels_gp = gpar(fontsize = 7),
                                 link_width = unit(7, "mm"))),
        heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                    labels_gp = gpar(fontsize=7),
                                    grid_width = unit(2.5, "mm")))

hm2 = Heatmap(scdata.merge.scaled, 
        name = "Exp",
        cluster_rows = T, cluster_columns = F,
        col = rev(brewer.pal(n=9, "RdBu")),
        column_names_gp = gpar(fontsize=7),
        show_row_names = F,
        show_row_dend = F,
        width = unit(2.5, "cm"), height = unit(10,"cm"),
        right_annotation = rowAnnotation(
          hightlight = anno_mark(at = previous.genes.index, labels = previous.genes,
                                 labels_gp = gpar(fontsize = 7),
                                 link_width = unit(7, "mm"))),
        heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                    labels_gp = gpar(fontsize=7),
                                    grid_width = unit(2.5, "mm")))

## 1) Deconvolution by correlation
# Strongly not recommended because of the correlations within bulk samples
genes.used = unique(intersect(rownames(bulk.norm.merge), as.character(DEGenes$gene)))
colnames(scdata.merge) = c("G0","G1","G2","G3","G4")
Heatmap(t(cor(x=bulk.norm.merge[genes.used,], 
              y=scdata.merge[genes.used,], method = "spearman")), 
        cluster_rows = F,
        cluster_columns = F,
        col = circlize::colorRamp2(c(-0.2,0,0.2,0.4,0.6), rev(brewer.pal(5, "RdBu"))),
        heatmap_legend_param = list(title = "Spearman correlation"))
Heatmap(cor(scdata.merge[genes.used,],method = "spearman"), 
        cluster_rows = F,
        cluster_columns = F,
        col = circlize::colorRamp2(c(-0.2,0,0.2,0.4,0.6), rev(brewer.pal(5, "RdBu"))),
        heatmap_legend_param = list(title = "Spearman correlation"))
Heatmap(t(c), col = brewer.pal(n=9, "Reds"), name = "Proportion", 
        cluster_rows = F, cluster_columns = F)

## 2) Deconvolution by nnls regression
library(nnls)

# Scale each gene to have the same mean and variance
scdata.merge.scaled = t(scale(t(scdata.merge[genes.used,])))
bulk.norm.scaled = t(scale(t(bulk.norm[genes.used,])))
bulk.norm.merge.scaled = t(scale(t(bulk.norm.merge[genes.used,])))
deconv.result = c()

# Apply nnls regression
for (i in 1:dim(bulk.norm.scaled)[2]/3){
  out = nnls(
    rbind(scdata.merge.scaled,
          scdata.merge.scaled,
          scdata.merge.scaled),
    c(bulk.norm.scaled[,3*i-2],
      bulk.norm.scaled[,3*i-1],
      bulk.norm.scaled[,3*i]))
  deconv.result = c(deconv.result, out)
}

# deconv.result
cmatrix = rbind(c(1.047437,0.6470013,0,0.4004393,0.3953729),
                c(0.3011903,0.8216298,0.1412882,0,0.4173733),
                c(0.1861662,0.3115887,0.7696528,0.1423949,0),
                c(0.1994302,0,0.8283148,0.491128,0.3154364),
                c(0.04599633,0,0.04096413,0.7462578,0.6520374))
rownames(cmatrix) = c("MB", "PM", "MC", "MM", "Neu")
colnames(cmatrix) = c("G0", "G1", "G2", "G3", "G4")

# Optional: Scale the composition matrix by group sizes
# Idea inspired by Cibersort
scale.prop = matrix(table(df.scClass.non5$cluster)/dim(df.scClass.non5)[1])[,1]
cmatrix = t(cmatrix)/scale.prop
c = t(cmatrix)/colSums(cmatrix)

# Visualization
hm1 = Heatmap(t(c), col = brewer.pal(n=9, "OrRd"), name = "Proportion", 
              cluster_rows = F, cluster_columns = F,
              width = unit(2.65, "cm"), height = unit(6.5, "cm"),
              row_names_side = "left",
              heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                          labels_gp = gpar(fontsize=7),
                                          grid_height = unit(0.2, "cm"),
                                          legend_direction = "horizontal"),
              column_names_rot = 90,
              row_names_gp = gpar(fontsize=7),
              column_names_gp = gpar(fontsize=7))
draw(hm1, heatmap_legend_side="top")
decorate_heatmap_body("Proportion", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.5))
})

## 3) Bulk deconvolution using nonnegative lasso
library(dplyr)
library(glmnet)

Deconv.nnls <- function(x, y){
  require(glmnet)
  deconv = c()
  for (i in 1:ncol(y)){
    # Under this setting it's equivalent to nnls regression
    mod = glmnet(x, y[,i], lambda = 0, lower.limits = 0, intercept = TRUE)
    deconv = cbind(deconv, coef(mod))
  }
  colnames(deconv) = colnames(y)
  avg = ave()
  return(as.matrix(deconv))
}

Deconv.nnls(x=scdata.merge.scaled, y=bulk.norm.merge.scaled)

## ---------------------------------------------------
## Part 4: Bayesian Classifier for sc annotation using
##   bulk as reference (Zilionis et al. 2019, Immunity)
## ---------------------------------------------------

# Bayesian Clasifier function
BayesianClassifier <- function(reference, query, prior=NULL){
  require(Matrix)
  genes = intersect(rownames(reference), rownames(query))
  reference = reference[genes,]
  query = query[genes,]
  
  P_ij = t(reference)/colSums(reference)
  llh = P_ij %*% query
  
  if (is.null(prior)){
    return(t(llh)/colSums(llh))
  }else {
    llh = llh * prior
    return(t(llh)/colSums(llh))
  }
} 

# Map each single cells to bulk groups
llh.sc = BayesianClassifier(reference = bulk.norm.merge, 
                            query = ref.merge.neu_all@raw.data)
llh.sc = llh.sc[rownames(ref.merge.neu_all@meta.data),]
bay.group = apply(llh.sc, 1, function(x){colnames(llh.sc)[which.max(x)]})

# Confusion matrix
df.scClass = data.frame(llh.sc, ref.merge.neu_all@meta.data[rownames(llh.sc),], 
                        bay.group=bay.group)
cmatrix = as.matrix(table(df.scClass$bay.group, df.scClass$cluster))
c = matrix(cmatrix, nrow=dim(cmatrix)[1], dimnames = list(rownames(cmatrix), colnames(cmatrix)))
c = t(c)/colSums(c) 

# Visualizations 
Heatmap(c, col = brewer.pal(n=9, "OrRd"), name = "Percentage", 
        cluster_rows = F, cluster_columns = F)

df.scClass.non5 = df.scClass[!df.scClass$cluster %in% c("G5a", "G5b", "G5c", "G3b"),]

p1 = ggplot(df.scClass.non5, aes(x = cluster, y=MB, fill=cluster)) + 
  geom_violin(scale="width") + theme_bw() + scale_fill_manual(values=ann.color$Cluster) + 
  ggtitle("MB posterior") +
  ylim(0,0.4) + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=7, hjust = 0.5),
        axis.text = element_text(size=7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p2 = ggplot(df.scClass.non5, aes(x = cluster, y=PM, fill=cluster)) + 
  geom_violin(scale="width") + theme_bw() + scale_fill_manual(values=ann.color$Cluster) +
  ggtitle("PM posterior") +
  ylim(0,0.4) + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=7, hjust = 0.5),
        axis.text = element_text(size=7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p3 = ggplot(df.scClass.non5, aes(x = cluster, y=MC, fill=cluster)) + 
  geom_violin(scale="width") + theme_bw() + scale_fill_manual(values=ann.color$Cluster) +
  ggtitle("MC posterior") +
  ylim(0,0.4) + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=7, hjust = 0.5),
        axis.text = element_text(size=7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p4 = ggplot(df.scClass.non5, aes(x = cluster, y=MM, fill=cluster)) + 
  geom_violin(scale="width") + theme_bw() + scale_fill_manual(values=ann.color$Cluster) +
  ggtitle("MM posterior") +
  ylim(0,0.4) + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=7, hjust = 0.5),
        axis.text = element_text(size=7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')
p5 = ggplot(df.scClass.non5, aes(x = cluster, y=Neu, fill=cluster)) + 
  geom_violin(scale="width") + theme_bw() + scale_fill_manual(values=ann.color$Cluster) +
  ggtitle("Neu posterior") +
  ylim(0,0.4) + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=7, hjust = 0.5),
        axis.text = element_text(size=7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')

pdf("E:/Figures/Figure 2/FS2B.pdf", width = 7.28, height = 1.34)
cowplot::plot_grid(p1, p2, p3, p4, p5, nrow = 1, align = "hv")
dev.off()


bulk.mb.cells = rank(df.scClass.non5$MB) > 5477
table(df.scClass.non5[bulk.mb.cells,]$cluster)
bulk.pm.cells = rank(df.scClass.non5$PM) > 5477
table(df.scClass.non5[bulk.pm.cells,]$cluster)
bulk.mc.cells = rank(df.scClass.non5$MC) > 5477
table(df.scClass.non5[bulk.mc.cells,]$cluster)
bulk.mm.cells = rank(df.scClass.non5$MM) > 5477
table(df.scClass.non5[bulk.mm.cells,]$cluster)
bulk.neu.cells = rank(df.scClass.non5$Neu) > 5477
table(df.scClass.non5[bulk.neu.cells,]$cluster)

cmatrix = rbind(c(394,96,0,1,9),
                c(303,185,0,2,10),
                c(0,1,242,257,0),
                c(0,0,182,318,0),
                #c(0,0,0,26,28,1427),
                c(0,0,0,12,488))
rownames(cmatrix) = c("MB", "PM", "MC", "MM", "Neu")
colnames(cmatrix) = c("G0", "G1", "G2", "G3a", "G4")

scale.prop = matrix(table(df.scClass.non5$cluster)/dim(df.scClass.non5)[1])[,1]
cmatrix = t(cmatrix)/scale.prop
c = t(cmatrix)/colSums(cmatrix)
hm1 = Heatmap(c, col = brewer.pal(n=9, "OrRd"), name = "Percentage", 
        cluster_rows = F, cluster_columns = F,
        width = unit(1.4, "inch"), height = unit(1.2, "inch"),
        bottom_annotation = HeatmapAnnotation(Cluster = c("G0", "G1", "G2", "G3a", "G4"),
                                              which = "col", col=ann.color, 
                                              simple_anno_size = unit(0.1, "inch"),
                                              show_legend = F,
                                              show_annotation_name = F),
        heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                    labels_gp = gpar(fontsize=7),
                                    grid_height = unit(0.1, "inch"),
                                    legend_direction = "horizontal"),
        column_names_rot = 0,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=7))
draw(hm1, heatmap_legend_side = "top")

tapply(df.scClass.non5$MB, df.scClass.non5$cluster, mean)
tapply(df.scClass.non5$PM, df.scClass.non5$cluster, mean)
tapply(df.scClass.non5$MC, df.scClass.non5$cluster, mean)
tapply(df.scClass.non5$MM, df.scClass.non5$cluster, mean)



## ---------------------------------------------------
## Part 5: Other related figures in the paper
## ---------------------------------------------------

# Marker gene expressions
markers.b = c("Mpo","Elane","Ctsg","Prtn3","Prss57",
              "Ctsc","Camp","Ltf","Cybb","Cyba",
              "Lcn2","B2m","Mmp8","Mmp9","Hp",
              "Slpi","Cd63","Itgam","Fcgr3",
              "Cd177","Mmp25","Fpr1","Scamp1","Vamp2")
markers.c = c("Mki67","Ranbp1","Spc24","Anp32b","H2afx",
              "Ccne2","Cdc20","Smc2","Pmf1","Cks2",
              "Top2a","Esco2","Cks1b","Cenpf","Cdkn3",
              "Tuba1c","Kif22", "Stmn1","Smc4")

hm.b = Heatmap(t(scale(t(bulk.norm.merge[markers.b,]))), 
              name = "Exp",
              cluster_rows = F, cluster_columns = F,
              show_row_names = T,
              col = rev(brewer.pal(n=9, "RdBu")),
              column_names_gp = gpar(fontsize=7),
              row_names_gp = gpar(fontsize=7),
              width = unit(2.5, "cm"), height = unit(5,"cm"),
              heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                          labels_gp = gpar(fontsize=7),
                                          grid_width = unit(2.5, "mm")))

hm.c = Heatmap(t(scale(t(bulk.norm.merge[markers.c,]))), 
               name = "Exp",
               cluster_rows = T, cluster_columns = F,
               show_row_names = T,
               show_row_dend = F,
               col = rev(brewer.pal(n=9, "RdBu")),
               column_names_gp = gpar(fontsize=7),
               row_names_gp = gpar(fontsize=7),
               width = unit(2.5, "cm"), height = unit(5,"cm"),
               heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                           labels_gp = gpar(fontsize=7),
                                           grid_width = unit(2.5, "mm")))

# Heatmap for all DEGenes 
all.degenes = as.character(DEGenes$gene[
  DEGenes$gene %in% rownames(bulk.norm.merge)[rowSums(bulk.norm.merge)>0]])

bulk.alldeg.scaled = t(scale(t(bulk.norm.merge[all.degenes,])))
scdata.alldeg.scaled = t(scale(t(scdata.merge[all.degenes,])))
nonna.rows = rowSums(is.na(scdata.alldeg.scaled))==0 & 
  rowSums(is.na(bulk.alldeg.scaled)) == 0
bulk.alldeg.scaled = bulk.alldeg.scaled[nonna.rows,]
scdata.alldeg.scaled = scdata.alldeg.scaled[nonna.rows,]

hm.sc.alldeg = Heatmap(scdata.alldeg.scaled, 
               name = "Exp",
               cluster_rows = T, cluster_columns = F,
               show_row_names = F,
               show_row_dend = F,
               col = rev(brewer.pal(n=9, "RdBu")),
               column_names_gp = gpar(fontsize=7),
               row_names_gp = gpar(fontsize=7),
               column_title = "Single-cell\ngroups",
               column_title_gp = gpar(fontsize=7),
               width = unit(1.5, "cm"), height = unit(6.5,"cm"),
               heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                           labels_gp = gpar(fontsize=7),
                                           grid_width = unit(2.5, "mm")))
hm.bulk.alldeg = Heatmap(bulk.alldeg.scaled, 
               name = "Exp",
               cluster_rows = F, cluster_columns = F,
               show_row_names = F,
               show_row_dend = F,
               col = rev(brewer.pal(n=9, "RdBu")),
               column_names_gp = gpar(fontsize=7),
               row_names_gp = gpar(fontsize=7),
               column_title = "Morphological\ngroups",
               column_title_gp = gpar(fontsize=7),
               width = unit(1.5, "cm"), height = unit(6.5,"cm"),
               heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                           labels_gp = gpar(fontsize=7),
                                           grid_width = unit(2.5, "mm")))
hm.sc.alldeg + hm.bulk.alldeg