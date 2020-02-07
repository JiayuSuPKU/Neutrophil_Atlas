## This file contains code for figures related to GRN analysis

setwd("paths")
load("./AUCell.RData")

library(ggplot2)
library(ComplexHeatmap)

ann.color <- list(Cluster = c(G0=pal_npg()(10)[1], 
                              G1=pal_npg()(10)[2], 
                              G2=pal_npg()(10)[3], 
                              G3a=pal_npg()(10)[4], 
                              G3b=pal_npg()(10)[5], 
                              G4=pal_npg()(10)[6], 
                              G5a=pal_npg()(10)[7], 
                              G5b=pal_npg()(10)[8], 
                              G5c=pal_npg()(10)[9], 
                              Mono=pal_npg()(10)[10],
                              DC="gray9",
                              B="#7A57D1",
                              HSPC="#0000A1",
                              "T"="#1F6ED4"))

p1 = ggplot(df.umap, aes(x=UMAP.1, y=UMAP.2, color=cluster_union_all)) + 
  scale_color_manual(values=ann.color$Cluster) +
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
        legend.key = element_blank())

p2 = ggplot(df.umap, aes(x=UMAP.1, y=UMAP.2, color=condition)) + 
  geom_point(alpha=1, size=0.3) + theme_bw() +
  guides(colour=guide_legend(override.aes = list(size=5))) +
  xlab("UMAP 1") + ylab("UMAP 2") + labs(color="Condition") +
  scale_color_manual(values = c(brewer.pal(9, "Set1")[9], 
                                brewer.pal(9, "Set1")[5])) + 
  theme(panel.background = element_rect(fill=NA, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14), 
        legend.key = element_blank())

cowplot::plot_grid(p1, p2, align = "hv")

load("E:/Score.RData")
p1 = ggplot(data.frame(merge_ref_wt_eco@dr$umap@cell.embeddings), 
       aes(x=UMAP1, y=UMAP2, color=merge_ref_wt_eco@data["Xbp1",])) + 
  geom_point(alpha=0.7, size=0.5) +
  labs(color="Il1r1") + theme_bw() + 
  scale_color_gradient(low="light grey", high = "red") 
p2 = ggplot(data.frame(merge_ref_wt_eco@dr$umap@cell.embeddings), 
            aes(x=UMAP1, y=UMAP2, color=merge_ref_wt_eco@data["Il1r2",])) + 
  geom_point(alpha=0.7, size=0.5) +
  labs(color="Il1r2") + theme_bw() +
  scale_color_gradient(low="light grey", high = "red")
cowplot::plot_grid(p1, p2, align = "hv")

IL1R.regulons = c("Nfic...", "Trp53...", "Bcl3...", "Cebpb...", 
                  "Ets1...", "Fosl2...", "Irf7...", "Irf9...",
                  "Junb...", "Ltf...", "Nfe2...", "Nfil3...",
                  "Spi1...")
hm2 = Heatmap(challenge_compare_t_mat[rownames(challenge_compare_t_mat) %in% 
                                        collab[IL1R.regulons],], 
              cluster_columns = F,
              col = circlize::colorRamp2(c(-50,0,50), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize=12), 
              heatmap_legend_param=list(title="Regulon activity\n (t-value challenged group\n versus control)",
                                        legend_direction = "horizontal"))
draw(hm2, heatmap_legend_side = "top")


ggplot(df.umap.wt_nneu, aes(x=UMAP.1, y=UMAP.2, color=Xbp1...)) + 
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
