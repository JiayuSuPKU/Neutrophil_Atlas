## This file is for downstream network analysis & visualizations
## Input: 
##   1) regulons inferred by pySCENIC
##   2) AUCell matrix calculated by pySCENIC

setwd("paths")
load("./AUCell.RData")

## ---------------------------------------------------
## Part 1: Load and preprocess data
## ---------------------------------------------------

# Load and preprocess regulons
regulons = read.table("regulons_weights.new.txt", header = F, row.names = 1, 
                      sep = "\t", stringsAsFactors = F)
colnames(regulons) = c("TF", "Motif", "Score", "nTarget_Gene","Target_Gene")
write.csv(regulons, "regulons_weights.csv")

regulons$auc_names = gsub("[:(+):]", ".", rownames(regulons))
regulons[c("1810024B03Rik(+)","2010315B03Rik(+)","9130023H24Rik(+)"),]$auc_names = c("X1810024B03Rik...","X2010315B03Rik...","X9130023H24Rik...")
table(regulons$auc_names %in% colnames(auc))
regulons$Formal_Name = paste(regulons$TF, sprintf("(%d genes)",regulons$nTarget_Gene))
rownames(regulons) = regulons$auc_names

# Create a mapping dataframe for different versions of regulon names
collab = regulons[colnames(auc),]$Formal_Name
names(collab) = colnames(auc)

# Load AUCell matrix
auc.wt = read.csv("./AUC_mtx_wt.csv", header = T, row.names=1)
auc.ecoli = read.csv("./AUC_mtx_ecoli.csv", header = T, row.names=1)
auc.nneu = read.csv("./AUC_mtx_nneu.csv", header = T, row.names=1)

# # Cell names might be different in different versions of metadata 
# rownames(auc.nneu) = gsub("ctl", "ctrl", rownames(auc.nneu))

auc = rbind(auc.wt, auc.ecoli, auc.nneu)
cells = rownames(auc)

## ---------------------------------------------------
## Part 2: Dimension reduction of the AUCell matrix
## ---------------------------------------------------
library(umap)
library("RColorBrewer")
library(ggplot2)

## UMAP for wt neutrophils and non-neutrophils
auc.wt_nneu = auc[rownames(df.umap)[!is.na(df.umap$cluster_all)],]
cells.wt_nneu = rownames(auc.wt_nneu)

# Dimension reduction
auc.pca.wt_nneu = prcomp(auc.wt_nneu)
auc.umap.wt_nneu = umap(auc.pca.wt_nneu$x[,1:20])
df.umap.wt_nneu = data.frame(UMAP=auc.umap.wt_nneu$layout[cells.wt_nneu,], 
                             auc.wt_nneu, cell.info[cells.wt_nneu,])

# Add annotations for later visualizations
df.umap.wt_nneu$cluster_all = as.character(df.umap.wt_nneu$cluster)
df.umap.wt_nneu[rownames(auc.nneu),]$cluster_all = 
  as.character(meta.data[gsub("ctrl", "ctl", rownames(auc.nneu)),]$cell_type)

# K-means clustering
cell.kmc = kmeans(auc.pca.wt_nneu$x[,1:20], centers=7, iter.max=2000)
t = table(cell.kmc$cluster, df.umap.wt_nneu$cluster_all)
cmatrix = matrix(as.numeric(t), 7)
colnames(cmatrix) = colnames(t)
cmatrix = cmatrix[,c(1,2,13,12,11,3,4,5,6,7,8,9,10)]

# Visualizations
c = apply(cmatrix, 2, function(x) x/sum(x))
c = c[c(4,7,1,3,6,5,2),]
rownames(c) = c("R1","R2","R3","R4","R5","R6","R7")
ca = HeatmapAnnotation(Cluster=rownames(c), col=list(Cluster = c(
  R1=pal_jama()(7)[1], 
  R2=pal_jama()(7)[2],
  R3=pal_jama()(7)[3],
  R4=pal_jama()(7)[4], 
  R5=pal_jama()(7)[5],
  R6=pal_jama()(7)[6],
  R7=pal_jama()(7)[7])),
  show_annotation_name = F, show_legend = F,
  simple_anno_size = unit(0.2, "cm"))
hm.c = Heatmap(t(c), name = "Percentage", col = brewer.pal(n = 9, name = "Reds"),
        cluster_rows = F, cluster_columns = F,
        row_names_side = "left",
        column_names_side = "top",
        column_names_rot = 0,
        top_annotation = ca,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=7),
        width = unit(4, "cm"), height = unit(4, "cm"),
        #column_title = "K-means clusters using regulon activity",
        #column_title_gp = gpar(fontsize=7),
        #row_title = "Seurat clusters using gene expression",
        #row_title_gp = gpar(fontsize=7),
        heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                    labels_gp = gpar(fontsize=7),
                                    grid_width = unit(0.2, "cm"),
                                    legend_direction="vertical"))
draw(hm.c, heatmap_legend_side="right")
decorate_heatmap_body("Percentage", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.5))
})

## Rename regulon-based clusters
df.umap.wt_nneu$R_cluster = cell.kmc$cluster
df.umap.wt_nneu$R_cluster[cell.kmc$cluster == 4] = "R1"
df.umap.wt_nneu$R_cluster[cell.kmc$cluster == 7] = "R2"
df.umap.wt_nneu$R_cluster[cell.kmc$cluster == 1] = "R3"
df.umap.wt_nneu$R_cluster[cell.kmc$cluster == 3] = "R4"
df.umap.wt_nneu$R_cluster[cell.kmc$cluster == 6] = "R5"
df.umap.wt_nneu$R_cluster[cell.kmc$cluster == 5] = "R6"
df.umap.wt_nneu$R_cluster[cell.kmc$cluster == 2] = "R7"

pdf("E:/Figures/Figure 5/F5B_v2.pdf", width = 4, height = 4.92)
ggplot(df.umap.wt_nneu, aes(x=UMAP.1, y=UMAP.2, color=R_cluster)) + 
  #scale_color_manual(values=ann.color$Cluster) +
  scale_color_jama()+
  geom_point(alpha=1, size=0.3) + theme_bw() +
  guides(colour=guide_legend(override.aes = list(size=5))) +  
  xlab("UMAP 1") + ylab("UMAP 2") + labs(color="Cluster") +
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14), 
        legend.position = "top",
        legend.key = element_blank())
dev.off()

## UMAP for all neutrophils and non-neutrophils
auc.pca = prcomp(auc)
auc.umap = umap(auc.pca$x[,1:20])

# Visualizations
df.umap = data.frame(UMAP=auc.umap$layout[cells,], auc, cell.info[cells,])
# `cluster_union_all` contains information of control + ecoli groups
df.umap$cluster_union_all = df.umap$cluster_union
df.umap[rownames(auc.nneu),]$cluster_union_all = 
  as.character(meta.data[gsub("ctrl", "ctl", rownames(auc.nneu)),]$cell_type)
# `cluster_all` contains information of control groups
df.umap$cluster_all = as.character(df.umap$cluster)
df.umap[rownames(auc.nneu),]$cluster_all = 
  as.character(meta.data[gsub("ctrl", "ctl", rownames(auc.nneu)),]$cell_type)

# Remove low-quality cells
na.cells = rownames(df.umap)[is.na(df.umap$cluster_union_all)] 
df.umap = df.umap[!is.na(df.umap$cluster_union_all),]

ggplot(df.umap, aes(x=UMAP.1, y=UMAP.2, color=cell_type)) + 
  geom_point(alpha=0.8, size=0.5) + theme_bw()
ggplot(df.umap, aes(x=UMAP.1, y=UMAP.2, color=condition)) + 
  geom_point(alpha=0.8, size=0.5) + theme_bw()
ggplot(df.umap, aes(x=UMAP.1, y=UMAP.2, color=cluster_union_all)) + 
  geom_point(alpha=0.8, size=0.5) + theme_bw()

ggplot(df.umap, aes(x=UMAP.1, y=UMAP.2, color=Nfil3...)) + 
  geom_point(alpha=1, size=0.3) + 
  scale_color_material("indigo") +
  #ggtitle("Nfil3 Activity (119 target genes)") + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=15, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none')

## ---------------------------------------------------
## Part 3: Differentially Expression Analysis using 
##   generalized linear regression model
## ---------------------------------------------------
library(ComplexHeatmap)

## 1) Compare non-neu and neutrophils
df.umap$nneu_compare = df.umap$cluster_all
df.umap[rownames(auc.nneu),]$nneu_compare = "_NonNeu"
model = glm(as.formula("Cebpe...~nneu_compare"), data=df.umap)

# Apply glm model
nneu_compare_t_list = list()
for (i in colnames(auc)){
  f = as.formula(paste(i, "nneu_compare", sep = "~"))
  model = glm(f, data = df.umap)
  nneu_compare_t_list[[i]] = summary(model)$coefficients[,c(3)]
}
nneu_compare_t_mat = matrix(unlist(nneu_compare_t_list), 
                            nrow = dim(auc)[2], byrow = T)
rownames(nneu_compare_t_mat) = collab
#rownames(nneu_compare_t_mat) = colnames(auc)
colnames(nneu_compare_t_mat) = c("Intercept.nneu", "G0", "G1", "G2", "G3",
                                 "G4", "G5a", "G5b", "G5c")

# Save results
neu_up = c("E2f2...","Rad21...","Mafg...", "Cebpe...","Klf5...",
           "Ltf...","Klf7...","Mlx...", "Xbp1...","Nfil3...",
           "Cebpb...","Spi1...","Max...","Nfe2...", "Mxd1...")
neu_down = c("Myc...", "Cebpz...","Bclaf1...", "Trp53...")
write.csv(regulons[c(neu_up, neu_down),1:5], "Neu_specific_regulons.csv", row.names = F)

# Visualizations
color.fun = circlize::colorRamp2(
  c(-100,-75,-50,-25,0,25,50,75,100), 
  rev(c("#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7",
        "#D1E5F0","#92C5DE","#4393C3","#2166AC")))
color.fun2 = circlize::colorRamp2(
  c(-50,-25,0,25,50), 
  rev(c("#B2182B","#D6604D","#F7F7F7","#4393C3","#2166AC")))
hm1 = Heatmap(nneu_compare_t_mat[,-c(1)], 
              row_names_gp = gpar(fontsize=5), cluster_columns = F,
              heatmap_legend_param=list(title="Regulon activity\n (t-value group\n versus Non-Neu)",
                                        legend_direction = "horizontal"))
hm1 = Heatmap(nneu_compare_t_mat[rownames(nneu_compare_t_mat) %in% collab[c(neu_up, neu_down)],-c(1)], 
        row_names_gp = gpar(fontsize=7), cluster_columns = F,
        column_names_gp = gpar(fontsize=7),
        show_row_dend = F,
        #col = rev(brewer.pal(n=9,name="RdBu")),
        col = color.fun,
        width = unit(2.5, "cm"), height = unit(5.7, "cm"),
        heatmap_legend_param=list(title="Regulon activity\n (t-value group\n versus Non-Neu)",
                                  legend_direction = "vertical", 
                                  title_gp = gpar(fontsize=7),
                                  labels_gp = gpar(fontsize=7),
                                  grid_width = unit(0.2, "cm")))
        #left_annotation = rowAnnotation(foo2 = factor(c(rep(1, 15), rep(0,4)))))
draw(hm1, heatmap_legend_side = "left")

## 2) Compare wt and challenged neutrophils
challenge_compare_t_mat = matrix(0, nrow = dim(auc)[2], ncol = 8)
colnames(challenge_compare_t_mat) = c("G0", "G1", "G2", "G3",
                                      "G4", "G5a", "G5b", "G5c")
rownames(challenge_compare_t_mat) = collab

# Apply glm model
for (i in c("G0", "G1", "G2", "G3",
            "G4", "G5a", "G5b", "G5c")){
  for (j in colnames(auc)){
    f = as.formula(paste(j, "condition", sep = "~"))
    model = glm(f, data = df.umap[df.umap$cluster_union==i,])
    challenge_compare_t_mat[j,i] = summary(model)$coefficient[2,3]
  }
}
# Regulons were assigned to different categories by t-value thresholds (skipped here)

# Load selected regulons for visualization
challenged_up = scan("Challenged_regulon_up.txt", character())
challenged_down = scan("Challenged_regulon_down.txt", character())
challenged_up_down = scan("Challenged_regulon_up_down.txt", character())
challenged_cluster = data.frame(names=c(challenged_up, challenged_up_down, challenged_down))
rownames(challenged_cluster) = challenged_cluster$names
challenged_cluster$type = "Up"
challenged_cluster[challenged_down,]$type = "Down"
challenged_cluster[challenged_up_down,]$type = "Up_Down"
challenged_cluster$type = factor(challenged_cluster$type)

write.csv(regulons[as.character(challenged_cluster$names),1:5],"Challenged_regulons.csv", row.names = F)  

mouse_DEGs_clusters_Ecoli = scan("mouse_DEGs_clusters_Ecoli.txt", character())

# Visualizations
hm2 = Heatmap(challenge_compare_t_mat, cluster_columns = F,
              row_names_gp = gpar(fontsize=5), 
              heatmap_legend_param=list(title="Regulon activity\n (t-value challenged group\n versus control)",
                                        legend_direction = "horizontal"))
row_anno = HeatmapAnnotation(Type=challenged_cluster$type, which = "row", 
                             col = list(Type=c("Up"="darkorchid3",
                                               "Down"="chartreuse4", 
                                               "Up_Down"="gold")),
                             simple_anno_size = unit(0.2, "cm"),
                             show_legend = F, show_annotation_name = F)
hm2 = Heatmap(challenge_compare_t_mat[collab[rownames(challenged_cluster)],], 
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize=7),
              column_names_gp = gpar(fontsize=7),
              width = unit(4.4, "cm"), height = unit(11, "cm"),
              row_names_side = "left",
              left_annotation = row_anno,
              col = color.fun2,
              heatmap_legend_param=list(title="Regulon activity\n (t-value challenged group\n versus control)",
                                        legend_direction = "horizontal",
                                        title_gp = gpar(fontsize=7),
                                        labels_gp = gpar(fontsize=7),
                                        grid_height = unit(0.2, "cm")))
draw(hm2, heatmap_legend_side = "top")

hm2 = Heatmap(dev_wt_compare_t_mat[collab[rownames(challenged_cluster)],1:7], 
              cluster_columns = F,
              cluster_rows = F, 
              show_row_names = F,
              width = unit(3.85, "cm"), height = unit(11, "cm"),
              left_annotation = row_anno,
              col = color.fun2,
              column_names_gp = gpar(fontsize=7), 
              column_names_rot = 45,
              heatmap_legend_param=list(title="Regulon activity\n (t-value group\n versus previous group)",
                                        legend_direction = "horizontal",
                                        title_gp = gpar(fontsize=7),
                                        labels_gp = gpar(fontsize=7),
                                        grid_height = unit(0.2, "cm")))
draw(hm2, heatmap_legend_side = "top")

## 3) Compare regulon activity before and after wt group transition events
GroupCompared_wt = rbind(c("G0_Control", "G1_Control"), 
                         c("G1_Control", "G2_Control"),
                         c("G2_Control", "G3a_Control"),
                         c("G3a_Control", "G4_Control"),
                         c("G4_Control", "G5a_Control"),
                         c("G4_Control", "G5b_Control"),
                         c("G4_Control", "G5c_Control"),
                         c("G5a_Control", "G5b_Control"),
                         c("G5a_Control", "G5c_Control"),
                         c("G5b_Control", "G5c_Control"))

# Apply glm model
dev_wt_compare_t_mat = matrix(0, nrow = dim(auc)[2], ncol = 10)
rownames(dev_wt_compare_t_mat) = colnames(auc)
for (i in 1:10){
  for (j in colnames(auc)){
    f = as.formula(paste(j, "cluster_condition", sep = "~"))
    model = glm(f, data = df.umap[df.umap$cluster_condition %in% GroupCompared_wt[i,],])
    dev_wt_compare_t_mat[j,i] = summary(model)$coefficient[2,3]
  }
}
colnames(dev_wt_compare_t_mat) = c("G1_vs_G0","G2_vs_G1","G3a_vs_G2","G4_vs_G3a",
                                   "G5a_vs_G4","G5b_vs_G4","G5c_vs_G4","G5b_vs_G5a",
                                   "G5c_vs_G5a","Gc_vs_G5b")
rownames(dev_wt_compare_t_mat) = collab

# Filter regulons based on t-values
keep.top.wt = rowSums(abs(dev_wt_compare_t_mat) > 40, na.rm = T)>0
table(keep.top.wt)
write.csv(regulons[colnames(auc)[keep.top.wt], 1:5], 
          "Dev_switch_regulons.csv", row.names = F)

# Visualizations
hm3 = Heatmap(dev_wt_compare_t_mat, cluster_columns = F,
              row_names_gp = gpar(fontsize=5), 
              heatmap_legend_param=list(title="Regulon activity\n (t-value group\n versus previous group)",
                                        legend_direction = "horizontal"))
hm3 = Heatmap(dev_wt_compare_t_mat[keep.top.wt,1:7], cluster_columns = F,
              width = unit(2.5, "cm"), height = unit(12.54, "cm"),
              show_row_dend = F,
              right_annotation = rowAnnotation(pattern = anno_block(
                gp = gpar(fill = brewer.pal(5, "Spectral"))), 
                width = unit(0.2, "cm")),
              row_km = 5, 
              row_names_gp = gpar(fontsize=7), 
              column_names_gp = gpar(fontsize=7),
              col = color.fun2,
              heatmap_legend_param=list(title="Regulon activity\n (t-value group\n versus previous group)",
                                        legend_direction = "horizontal",
                                        title_gp = gpar(fontsize=7),
                                        labels_gp = gpar(fontsize=7),
                                        grid_height = unit(0.2, "cm")))
draw(hm3, heatmap_legend_side = "top")


GroupCompared_ecoli = rbind(c("G0_Ecoli", "G1_Ecoli"), 
                            c("G1_Ecoli", "G2_Ecoli"), 
                            c("G2_Ecoli", "G3a_Ecoli"), 
                            c("G3a_Ecoli", "G4_Ecoli"), 
                            c("G4_Ecoli", "G5a_Ecoli"), 
                            c("G4_Ecoli", "G5b_Ecoli"), 
                            c("G4_Ecoli", "G5c_Ecoli"), 
                            c("G5a_Ecoli", "G5b_Ecoli"), 
                            c("G5a_Ecoli", "G5c_Ecoli"), 
                            c("G5b_Ecoli", "G5c_Ecoli"))

# ## 4) Compare transition events in ecoli sample (results not shown in the paper)
# dev_ecoli_compare_t_mat = matrix(0, nrow = dim(auc)[2], ncol = 10)
# rownames(dev_ecoli_compare_t_mat) = colnames(auc)
# for (i in 1:10){
#   for (j in colnames(auc)){
#     f = as.formula(paste(j, "cluster_condition", sep = "~"))
#     model = glm(f, data = df.umap[df.umap$cluster_condition %in% GroupCompared_ecoli[i,],])
#     dev_ecoli_compare_t_mat[j,i] = summary(model)$coefficient[2,3]
#   }
# }
# colnames(dev_ecoli_compare_t_mat) = c("G1_vs_G0","G2_vs_G1","G3a_vs_G2","G4_vs_G3a",
#                                       "G5a_vs_G4","G5b_vs_G4","G5c_vs_G4","G5b_vs_G5a",
#                                       "G5c_vs_G5a","Gc_vs_G5b")
# rownames(dev_ecoli_compare_t_mat) = collab
# hm4 = Heatmap(dev_ecoli_compare_t_mat, cluster_columns = F,
#               row_names_gp = gpar(fontsize=5), 
#               heatmap_legend_param=list(title="Regulon activity\n (t-value group\n versus previous group)",
#                                         legend_direction = "horizontal"))
# draw(hm4, heatmap_legend_side = "top")

# keep.top.ecoli = rowSums(abs(dev_ecoli_compare_t_mat) > 40, na.rm = T)>0
# table(keep.top.ecoli)
# table(keep.top.wt & keep.top.wt)

# hm4 = Heatmap(dev_ecoli_compare_t_mat[keep.top.wt,], 
#               cluster_columns = F,
#               row_names_gp = gpar(fontsize=12), 
#               heatmap_legend_param=list(title="Regulon activity\n (t-value group\n versus previous group)",
#                                         legend_direction = "horizontal"))
# draw(hm4, heatmap_legend_side = "top")


## ---------------------------------------------------
## Part 4: Calculate consistancy between wt 
##   development and E.coli challenged development
## ---------------------------------------------------

# AUC scoring helper function
recovery_curve <- function(reference, query, top=50){
  recovery_rate = c()
  for(i in 1:top){
    recovery_rate = c(recovery_rate, mean(reference %in% query[1:i]))
  }
  return(list(x=1:top, y=recovery_rate, auc=mean(recovery_rate)))
}

# Calculate consistency score defined by the AUC score,
# a metric that measures the overlaps between DE regulons/genes
# under different conditions
calc_c_score <- function(wt.mat, ecoli.mat, top=dim(wt.mat)[1]){
  up.score = c()
  for (i in 1:10){
    ref.up = colnames(auc)[wt.mat[,i] > 5]
    rank_ecoli_up = rownames(ecoli.mat)[order(ecoli.mat[,i], decreasing = T)]
    up.score = c(up.score, recovery_curve(reference = ref.up,
                                          query = rank_ecoli_up,
                                          top = top)$auc)
  }
  down.score = c()
  for (i in 1:10){
    ref.down = colnames(auc)[wt.mat[,i] < -5]
    rank_ecoli_down = rownames(ecoli.mat)[order(ecoli.mat[,i], decreasing = F)]
    down.score = c(down.score, recovery_curve(reference = ref.down,
                                              query = rank_ecoli_down,
                                              top = top)$auc)
  }
  
  return(list(up.score=up.score, down.score=down.score))
}

# Calculate Pearson correlation
calc_cor <- function(wt.mat, ecoli.mat){
  na.regulons = is.na(wt.mat) | is.na(ecoli.mat)
  corr = c()
  for (i in 1:dim(wt.mat)[2]){
    corr = c(corr, cor(wt.mat[,i][!na.regulons[,i]],
                       ecoli.mat[,i][!na.regulons[,i]]))
  }
  return(corr)
}

## calculate consistancy scores and correlation between two conditions
calc_c_score(wt.mat = dev_wt_compare_t_mat,
             ecoli.mat = dev_ecoli_compare_t_mat)
calc_cor(wt.mat = dev_wt_compare_t_mat,
         ecoli.mat = dev_ecoli_compare_t_mat)

# Visualizations
corr = cor(x=dev_wt_compare_t_mat[rowSums(na.regulons)==0,1:7],
           y=dev_ecoli_compare_t_mat[rowSums(na.regulons)==0,1:7])

Heatmap(t(corr),col = rev(brewer.pal(n = 9, name = "RdBu")),
        cluster_rows = F, cluster_columns = F,
        width = unit(1.15, "inch"), height = unit(1.15, "inch"),
        column_names_gp = gpar(fontsize=7),
        row_names_gp = gpar(fontsize=7),
        column_names_rot = 45,
        row_title = "Group switches\nafter E. coli challenge",
        column_title = "Group switches\nunder normal condition",
        column_title_gp = gpar(fontsize=7),
        row_title_gp = gpar(fontsize=7),
        heatmap_legend_param = list(title="Pearson\ncorrelation",
                                    title_gp = gpar(fontsize=7),
                                    labels_gp = gpar(fontsize=7),
                                    grid_width = unit(0.2, "cm")))

corr = cor(dev_wt_compare_t_mat[rowSums(na.regulons)==0,1:7])
Heatmap(corr,col = rev(brewer.pal(n = 9, name = "RdBu")),
        cluster_rows = F, cluster_columns = F,
        show_heatmap_legend = F,
        width = unit(1.15, "inch"), height = unit(1.15, "inch"),
        column_names_gp = gpar(fontsize=7),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize=7),
        column_title = "Group switches\nunder normal condition",
        row_title = "Group switches\nunder normal condition",
        column_title_gp = gpar(fontsize=7),
        row_title_gp = gpar(fontsize=7),
        heatmap_legend_param = list(title="Pearson\ncorrelation",
                                    title_gp = gpar(fontsize=7),
                                    labels_gp = gpar(fontsize=7),
                                    grid_width = unit(0.2, "cm")))