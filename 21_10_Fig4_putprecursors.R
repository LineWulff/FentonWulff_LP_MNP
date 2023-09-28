### Plotting for Figure 4 - putative DC precursors
## Author: Line Wulff
###############################################################
#######    F4 DC traj fig updates september 2021       ########
###############################################################
rm(list=ls())
setwd()

library(PairedData)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(tidyverse)
library(tidyr)
library(Seurat)
#library(Seurat, lib.loc = "/Library/Frameworks/R.framework/Versions/3.5/Resources/other_libs") #3.1.1
library(velocyto.R)
library(pagoda)
library(ggrastr)

#### Variables ####
##date in format year_month_day 
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
pseud_col <- c('magenta','darkorange','gold','black')
pseud_col <- myColorRamp(pseud_col,seq(1:50))

#### Read in necessary data ####
# reclustered DC object, with tSpace PCs and high resolution clustering
DCtraj <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata_v2/DCtraj_v5.rds")

DimPlot(DCtraj, group.by = "tSP_clus_comp", reduction = "flat_3D", label = T, label.size =6, pt.size = 2)

# divide subsets to overall populations from Fig 2
DC1A_lev <- as.character(c(43, 30, 25, 41, 8, 39, 26))
DC2_lev <- as.character(c(33, 0, 18, 46, 11, 52, 48, 34, 24, 20, 28, 22, 49, 3))
DC3_lev <-as.character(c(13, 27, 44, 37, 23, 17, 6, 42, 14, 38, 5, 21, 31))
ambig_lev <-as.character(c(15, 12, 10, 2, 4, 9, 45, 7, 16, 40, 1, 47))
HLA_lev <-as.character(c(53, 19, 36, 32, 50, 51, 35, 29))

DCtraj@meta.data$tSP_clus_comp_v2 <- NA
for (clus in unique(DCtraj@meta.data$integrated_snn_res.3)){
  if (clus %in% DC1A_lev){
    DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3==clus,]$tSP_clus_comp_v2 <-"cDC1"}
  else if (clus %in% DC2_lev){
    DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3==clus,]$tSP_clus_comp_v2 <-"cDC2"}
  else if (clus %in% DC3_lev){
    DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3==clus,]$tSP_clus_comp_v2 <-"cDC3"}
  else if (clus %in% ambig_lev){
    DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3==clus,]$tSP_clus_comp_v2 <-"amb"}
  else if (clus %in% HLA_lev){
    DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3==clus,]$tSP_clus_comp_v2 <-"HLA low"}
  else {print(paste(clus,"wasn't in levs."))}
}

DCtraj@meta.data$tSP_clus_comp_v2 <- factor(DCtraj@meta.data$tSP_clus_comp_v2,levels=c("cDC1","cDC2","amb","cDC3","HLA low"))
pdf(paste(dato,"tSP_DCtraj_compclust_v2_wlabels.pdf"), width=4, height=3.5)
DimPlot(DCtraj, group.by = "tSP_clus_comp_v2", reduction = "flat_3D", label = T, label.size =6, pt.size = 2)+NoLegend()
dev.off()
pdf(paste(dato,"tSP_DCtraj_compclust_v2_wlabels.pdf"), width=4, height=3.5)
DimPlot(DCtraj, group.by = "tSP_clus_comp_v2", reduction = "flat_3D")
dev.off()


#### A) w HLA low as separate ####
DCtraj@meta.data$tSP_clus_comp_v3 <- as.character(DCtraj@meta.data$tSP_clus_comp_v2)
DCtraj@meta.data[DCtraj@meta.data$tSP_clus_comp_v2=="HLA low",]$tSP_clus_comp_v3 <- as.character(DCtraj@meta.data[DCtraj@meta.data$tSP_clus_comp_v2=="HLA low",]$tSP_clus_comp)

#save colors as in tSP_clus_comp_v3
clus_v3_col <- c("cDC1"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
                 "53"="#F8766D","19"="#FF61CC","32"="#7CAE00","36"="#00BE67","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","51"="#CD9600")

pdf(paste(dato,"tSP_DCtraj_compclust_v3_wolabels.pdf"), width=4, height=3.5)
DimPlot(DCtraj, group.by = "tSP_clus_comp_v3", reduction = "flat_3D",pt.size = 1)+scale_color_manual(values=clus_v3_col)
dev.off()

pdf(paste(dato,"tSP_DCtraj_compclust_v3_wolabels.pdf"), width=4, height=3.5)
ggplot(cbind(DCtraj@meta.data,DCtraj@reductions$flat_3D@cell.embeddings), aes(x=UMAP_1,y=UMAP_2,colour=tSP_clus_comp_v3))+
  geom_point_rast()+theme_classic()+scale_color_manual(values = clus_v3_col)+
  xlab("Flat tUMAP1")+ylab("Flat tUMAP2")+
  theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()


pdf(paste(dato,"tSP_DCtraj_compclust_v3_wlabels.pdf"), width=4, height=3.5)
DimPlot(DCtraj, group.by = "tSP_clus_comp_v3", reduction = "flat_3D", label=T, label.size=6, pt.size = 1)+scale_color_manual(values=clus_v3_col)+NoLegend()
dev.off()


#### C) CD11c expression ####
mycols_blue <-  c("#bdbdbd","#d9d9d9","skyblue3","steelblue","steelblue4","#1B3346")

pdf(paste(dato,"tSP_DCtraj_CD11cgeneIntExp.pdf"), width=4, height=3.5)
ggplot(as.data.frame(cbind(DCtraj@reductions$flat_3D@cell.embeddings, t(as.matrix(DCtraj@assays$integrated@data)))), 
       aes(x=UMAP_1, UMAP_2, colour=ITGAX))+
  geom_point_rast()+theme_classic()+scale_color_gradientn(colours=mycols_blue)+
  xlab("Flat tUMAP1")+ylab("Flat tUMAP2")+
  theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()

#### B) Proliferation gene module ####
ccyc <- list(c('KIAA0101', 'TUBA1B', 'MKI67', 'HIST1H4C', 'TUBB', 'UBE2C', 'STMN1', 'H2AFZ'))
DCtraj <- AddModuleScore(DCtraj, features = ccyc, ctrl = 5, name = 'cc.score')

pdf(paste(dato,"tSP_DCtraj_cellcyclemodulescore.pdf"), width=4, height=3.5)
ggplot(as.data.frame(cbind(DCtraj@reductions$flat_3D@cell.embeddings, DCtraj@meta.data)), 
       aes(x=UMAP_1, UMAP_2, colour=cc.score1))+
  geom_point_rast()+theme_classic()+scale_color_gradientn(colours=mycols)+
  xlab("Flat tUMAP1")+ylab("Flat tUMAP2")+
  theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()

## same HLA score, using all HLA genes included in the integrated data set
HLA_genes <- list(rownames(DCtraj@assays$integrated@data)[startsWith(rownames(DCtraj@assays$integrated@data),"HLA")])
DCtraj <- AddModuleScore(DCtraj, features = HLA_genes, ctrl = 5, name = 'HLA.score')

pdf(paste(dato,"tSP_DCtraj_HLAmodulescore.pdf"), width=4, height=3.5)
ggplot(as.data.frame(cbind(DCtraj@reductions$flat_3D@cell.embeddings, DCtraj@meta.data)), 
       aes(x=UMAP_1, UMAP_2, colour=HLA.score1))+
  geom_point()+theme_classic()+scale_color_gradientn(colours=mycols)+
  xlab("Flat tUMAP1")+ylab("Flat tUMAP2")+
  theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()


#### E) Velocity calculation and embedding on Flat 3D tSP UMAP ####
# before this loom files were concantenated, renamed to fit object and split into spliced (emat) and unspliced (nmat) matrices 
nmat <- readRDS(file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R1_AllCells_nmat.rds")
emat <- readRDS(file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R1_AllCells_emat.rds")

nmat_all <- nmat
emat_all <- emat

## Velo with extra precursor clusters
clus_v3.5_col <- c("cDC1"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
                   "53"="#F8766D","19"="#FF61CC","43"="darksalmon","30"="rosybrown","32"="#7CAE00","36"="#00BE67","33"="palegreen3","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","51"="#CD9600")
DCtraj@meta.data$tSP_clus_comp_v3.5 <- as.character(DCtraj@meta.data$tSP_clus_comp_v3)
DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3 %in% c(43,30,33),]$tSP_clus_comp_v3.5 <- as.character(DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3 %in% c(43,30,33),]$integrated_snn_res.3)

# check colouring and included clusters fit
DimPlot(DCtraj, group.by = "tSP_clus_comp_v3.5", reduction = "flat_3D", label=T, label.size=6, pt.size = 1)+scale_color_manual(values=clus_v3.5_col)+NoLegend()

visu_obj <- DCtraj

Idents(visu_obj) <- visu_obj@meta.data$tSP_clus_comp_v3.5
Idents(visu_obj) <- factor(Idents(visu_obj), levels = names(clus_v3.5_col))
cluster.label <- Idents(visu_obj)
#specific colours from clustering
#add for each cell w. pagoda
cell.colors <- pagoda2:::fac2col(cluster.label,level.colors = clus_v3.5_col)

#and include only cells of interest
incl_cells <- intersect(colnames(visu_obj),colnames(emat))
emat <- emat[,incl_cells]
nmat <- nmat[,incl_cells]

### Calculation and embedding in visual spaces
#calculate cell-to-cell-distance matrix
cell.dist <- as.dist(1-armaCor(t(visu_obj@reductions$pca@cell.embeddings)))
#tSNE/UMAP or PCA embeddings from seurat
emb <- visu_obj@reductions$flat_3D@cell.embeddings[,c(1,2)]

#filtering of genes based on minimum average expression magnitude
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat))) #5262 genes after filtering if 0.2 

#genes for velocity calculation should overlap with integrated features
genes_velo <- intersect(rownames(emat),rownames(emat))
var_overlap <- intersect(genes_velo, rownames(visu_obj@assays$integrated@data))
emat.var <- emat[var_overlap,]
dim(emat.var)

genes_velo <- intersect(rownames(nmat),rownames(nmat))
var_overlap <- intersect(genes_velo, rownames(visu_obj@assays$integrated@data))
nmat.var <- nmat[var_overlap,]

# Manual gate, can't be reproduced exactly, instead use saved csv with list of cells to get exact same outputs
#DC1s <- CellSelector(DimPlot(visu_obj, group.by = "tSP_clus_comp_v3", reduction="flat_3D")+scale_color_manual(values=clus_v3_col))
#DC2s <- CellSelector(DimPlot(visu_obj, group.by = "tSP_clus_comp_v3", reduction="flat_3D")+scale_color_manual(values=clus_v3_col))
#DC3s <- CellSelector(DimPlot(visu_obj, group.by = "tSP_clus_comp_v3", reduction="flat_3D")+scale_color_manual(values=clus_v3_col))
#interest <- c(DC1s,DC2s,DC3s)
#write.csv(interest,"velocyto_cells_v3.csv")
interest <- read.csv("velocyto_cells_v3.csv")
interest <- as.character(interest$x)
length(interest) #4255 cells
nmat.int <- nmat.var[,interest]
emat.int <- emat.var[,interest]

# Estimation of RNA vel. with k=10% cells kNN pooling and top/bottom, 2% quantiles for gamma fit
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat.int,nmat.int,
                                            deltaT=1,kCells=440, 
                                            cell.dist=cell.dist,fit.quantile=fit.quantile, 
                                            zero.offset = T, n.cores = 4)

# Embed velocity estimation on existing flat 3d tSP umap embedding
pdf(paste(dato,"DCtraj_velocytoCellInterest_PearsonEmb_tSPclus3.5_v1.pdf",sep="_"), height = 5, width = 5)
show.velocity.on.embedding.cor(emb, rvel.cd, n = 400, scale = 'sqrt', 
                               cell.colors = ac(cell.colors, alpha = 0.70), cex = 1.5, n.cores = 4,
                               arrow.scale = 2, show.grid.flow = TRUE, min.grid.cell.mass = 0.1,
                               grid.n = 25, arrow.lwd = 0.8, do.par = T, cell.border.alpha = 0,
                               axes = FALSE)
dev.off()

# Without velo. for legend purposes
pdf(paste(dato,"DCtraj_legend_clus3.5_v1.pdf",sep="_"), height = 5, width = 5)
ggplot(cbind(DCtraj@meta.data,DCtraj@reductions$flat_3D@cell.embeddings), aes(x=UMAP_1,y=UMAP_2,colour=tSP_clus_comp_v3.5))+
  geom_point_rast()+theme_classic()+scale_color_manual(values = clus_v3.5_col)+
  xlab("Flat tUMAP1")+ylab("Flat tUMAP2")+
  theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()


#### F) PCA on comitted precursors ####
DEGs <- FindAllMarkers(subset(visu_obj, idents = c("19","29","36")), #c("19","29","36") c("cDC1A","cDC3","cDC2") c("36","53","35")
                       only.pos = T, #only including UPregulated DEGs per cluster
                       test.use = "wilcox", #test to use
                       min.pct = 0.25, # genes included in results should be expressed in atleast 25% of the cluster
                       logfc.threshold = 0.25)

DEGs <- DEGs[DEGs$p_val_adj<0.05,]
dim(DEGs)
top50 <- DEGs %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
# top 50 genes from each cluster based on avg_logFC
prec_genes <- top50$gene

# subset cells from HLA low clusters to include in PCA
Idents(DCtraj) <- 'tSP_clus_comp_v3'
precursors <- subset(DCtraj, idents = c("19","29","36","53","35","51","32","50")) #,"19","29","36","53","35"
precursors <- ScaleData(precursors, vars.to.regress = c("nCount_RNA", "percent.mt","G2M.Score","S.Score","cc.exp"), do.scale = TRUE) 
precursor_mat <- t(precursors@assays$integrated@scale.data[prec_genes,])

prec_pca <- prcomp(precursor_mat)
prec_df <- as.data.frame(prec_pca$x)
prec_pca_2 <- PCA(precursor_mat)
# add coluring labels
prec_df$cluster <- as.character(precursors@meta.data[rownames(prec_df),]$tSP_clus_comp_v3)

# check PC variance pr PC
plot(1:50, prec_pca_2$var[1:50]/prec_pca_2$totalvar, cex=2, ylab = "% Variance explained", xlab = "PCs")
lines(1:50, prec_pca_2$var[1:50]/prec_pca_2$totalvar)

# Plot only cells that showed commitment from F4D-E
clus_v3_alp <- c("53"= 1,"19"= 1,"32"= 0,"36"= 1,"35"=0,"29"= 1,"50"= 0,"51"= 0,"cDC1"=0,"cDC2"=0)
prec_df$alpha <- 0
for (clus in unique(prec_df$cluster)) {
  prec_df[prec_df$cluster==clus,]$alpha <- clus_v3_alp[clus]
}

pdf(paste(dato,"4F_precommColoured_w35.pdf"), width=4, height=3.5)
ggbiplot(prec_pca, var.axes = F, groups = prec_df$cluster, alpha=0)+
  geom_point_rast(aes(colour = prec_df$cluster), alpha=prec_df$alpha)+
  theme_classic()+
  scale_color_manual(values=clus_v3_col)+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

clus_v3_alp <- c("53"= 0.5,"19"= 0.5,"32"= 1,"36"= 0.5,"35"=1,"29"= 0.5,"50"= 1,"51"= 1)
prec_df$alpha <- 0
for (clus in unique(prec_df$cluster)) {
  prec_df[prec_df$cluster==clus,]$alpha <- clus_v3_alp[clus]
}
clus_v3_col_2 <- c("53"="grey","19"="grey","32"="#7CAE00","36"="grey","35"="#00BFC4","29"="grey","50"="#C77CFF","51"="#CD9600")

pdf(paste(dato,"4Gcomb_UncommColoured_w35.pdf"), width=4, height=3.5)
ggbiplot(prec_pca, var.axes = F, groups = prec_df$cluster, alpha=0)+
  geom_point_rast(aes(colour = prec_df$cluster), alpha=prec_df$alpha)+
  theme_classic()+
  scale_color_manual(values=clus_v3_col_2)+
  theme(axis.ticks = element_blank(),axis.text = element_blank())
dev.off()

### Add PC1&2 to precursors to be able to gate
precursors@reductions$prec_pca <- precursors@reductions$pca
precursors@reductions$prec_pca@cell.embeddings <- prec_pca$x[rownames(precursors@reductions$prec_pca@cell.embeddings),c(1,2)]
colnames(precursors@reductions$prec_pca@cell.embeddings) <- c("PC_1","PC_2")

DimPlot(precursors, reduction = "prec_pca",pt.size = 2)+scale_color_manual(values=clus_v3_col_2)#+scale_alpha_manual(values=clus_v3_alp)
#DC1_prec <- CellSelector(DimPlot(precursors, reduction = "prec_pca",pt.size = 2)+scale_color_manual(values=clus_v3_col_2))
#DC3_prec <- CellSelector(DimPlot(precursors, reduction = "prec_pca",pt.size = 2)+scale_color_manual(values=clus_v3_col_2))
#uncomm_prec <- CellSelector(DimPlot(precursors, reduction = "prec_pca",pt.size = 2)+scale_color_manual(values=clus_v3_col_2))

#DC1_prec <- DC1_prec[DC1_prec %in% rownames(precursors@meta.data[precursors@meta.data$tSP_clus_comp_v3=="50" | precursors@meta.data$tSP_clus_comp_v3=="35",])]
#DC3_prec <- DC3_prec[DC3_prec %in% rownames(precursors@meta.data[precursors@meta.data$tSP_clus_comp_v3=="50",])]

DC1_prec <- as.character(read.csv("DC1prec_cl5035.csv")$x)
DC1_prec <- as.character(read.csv("DC3prec_cl50.csv")$x)

unique(DCtraj@meta.data$tSP_clus_comp_v3)
DCtraj@meta.data$tSP_clus_comp_v5 <- as.character(DCtraj@meta.data$tSP_clus_comp_v3)
DCtraj@meta.data[DC1_prec,]$tSP_clus_comp_v5 <- "pre-cDC1"
DCtraj@meta.data[DC3_prec,]$tSP_clus_comp_v5 <- "pre-cDC3"

DimPlot(DCtraj, group.by = "tSP_clus_comp_v5", reduction = "flat_3D", pt.size = 2)

precursors@meta.data$tSP_clus_comp_v5 <- as.character(precursors@meta.data$tSP_clus_comp_v3)
precursors@meta.data[DC1_prec,]$tSP_clus_comp_v5 <- "pre-cDC1"
precursors@meta.data[DC3_prec,]$tSP_clus_comp_v5 <- "pre-cDC3"

DimPlot(precursors, group.by = "tSP_clus_comp_v5", reduction = "flat_3D", pt.size = 2)
ggplot(cbind(precursors@reductions$flat_3D@cell.embeddings, precursors@meta.data), 
       aes(x=UMAP_1,y=UMAP_2,colour=tSP_clus_comp_v5))+geom_point()+
  theme_classic()+scale_color_manual(values=c('#6a3d9a','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6'))

ggplot(cbind(precursors@reductions$flat_3D@cell.embeddings, precursors@meta.data), 
       aes(x=UMAP_1,y=UMAP_2,colour=tSP_clus_comp_v5))+geom_point()+
  theme_classic()+scale_color_manual(values=c('#6a3d9a','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6'))

cl_35 <- rownames(DCtraj@meta.data[DCtraj@meta.data$tSP_clus_comp_v3=="35",])
ggplot(cbind(precursors@reductions$flat_3D@cell.embeddings[cl_35,], precursors@meta.data[cl_35,]), 
       aes(x=UMAP_1,y=UMAP_2))+geom_point(colour="grey", alpha=0.5, size=2)+
  geom_point(data = cbind(precursors@reductions$flat_3D@cell.embeddings[DC1_prec,], precursors@meta.data[DC1_prec,]), colour="forestgreen")+
  theme_classic()

ggplot(cbind(DCtraj@reductions$flat_3D@cell.embeddings, DCtraj@meta.data), 
       aes(x=UMAP_1,y=UMAP_2))+geom_point(colour="grey", alpha=0.2, size=2)+
  geom_point(data = cbind(precursors@reductions$flat_3D@cell.embeddings[cl_35,], precursors@meta.data[cl_35,]), colour='#cab2d6')+
  geom_point(data = cbind(precursors@reductions$flat_3D@cell.embeddings[DC1_prec,], precursors@meta.data[DC1_prec,]), colour="forestgreen")+
  theme_classic()




#### G) HM with committed branches and sep. 50,51, 32 ####
HM_E_genes <- c("STMN1", "KIAA0101", "TUBA1B","MKI67", "HIST1H4C","LILRA4","GZMB","E2F1","CD38",
                "RUNX3","THBD","CLEC4C", "BATF3","BTLA","CADM1","XCR1","IRF8","ID2","CLEC9A",
                "CCR6","NR4A1","ATF3","SOX4","S100B","PRDM16","CXXC5","NGFR","CDH1","KIT","LTB",
                "IL22RA2","IRF4","PLAC8","CD207","AXL","CD1C","IL3RA","C1QA","CD14","CD163", "FCGR1A",
                "S100A9", "FCN1", "VCAN","HLA-DQA1","HLA-DRA")

# Av_visu <- AverageExpression(visu_obj, return.seurat = T, assay= "integrated")
# Idents(Av_visu) <- factor(Idents(Av_visu), levels = c("51","50","32","53","19","cDC1","36","cDC2","35","29","cDC3","amb"))
# 
# col_ord <- factor(names(clus_v3_col),levels = c("51","50","32","53","19","cDC1","36","cDC2","35","29","cDC3","amb"))
# 
# pdf(paste(dato,"G_HM_selcetedgenes_allgroups.pdf",sep="_"), width = 4, height = 8)
# DoHeatmap(Av_visu, features = HM_E_genes, group.colors = clus_v3_col[order(col_ord)], draw.lines = F)+
#   scale_fill_gradientn(colours = mycols, na.value = "white")
# dev.off()
# 
# Av_visu_sub <- subset(Av_visu, idents = c("53","19","cDC1","36","cDC2","35","29","cDC3"))
# col_ord_2 <- factor(names(clus_v3_col),levels = c("53","19","cDC1","36","cDC2","35","29","cDC3"))
# 
# pdf(paste(dato,"G_HM_selcetedgenes_comittedandmat.pdf",sep="_"), width = 4, height = 8)
# DoHeatmap(Av_visu_sub, features = HM_E_genes, group.colors = clus_v3_col[order(col_ord_2)], draw.lines = F)+
#   scale_fill_gradientn(colours = mycols, na.value = "white")
# dev.off()
DCtraj@meta.data$tSP_clus_comp_v4 <- as.character(DCtraj@meta.data$tSP_clus_comp_v3)
cl50 <- rownames(DCtraj@meta.data[DCtraj@meta.data$tSP_clus_comp_v3=="50",])
cl50_DC3 <- cl50[which(DCtraj@assays$integrated@data["FCER1A",cl50]>1)]
cl50_DC1 <- cl50[which(DCtraj@assays$integrated@data["CLEC9A",cl50]>1)]

DCtraj@meta.data[cl50_DC1,]$tSP_clus_comp_v4 <- "50_DC1"
DCtraj@meta.data[cl50_DC3,]$tSP_clus_comp_v4 <- "50_DC3"

Idents(DCtraj) <- 'tSP_clus_comp_v4'
clus_v4_col <- c("cDC1"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
                 "53"="#F8766D","19"="#FF61CC","32"="#7CAE00","36"="#00BE67","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","50_DC1"="#C77CFF","50_DC3"="#C77CFF","51"="#CD9600")
clus_ord <- rev(c("51","50_DC1","53","19","cDC1","32","36","cDC2","50_DC3","35","29","cDC3"))
Idents(DCtraj) <- factor(Idents(DCtraj), levels=rev(clus_ord))

Av_visu <- AverageExpression(subset(DCtraj, idents = names(clus_v4_col)[names(clus_v4_col) %in% clus_ord]), return.seurat = T, assay= "integrated")
Idents(Av_visu) <- factor(Idents(Av_visu), levels = clus_ord)

col_ord <- factor(names(clus_v4_col),levels =rev(clus_ord)) #clus_v3_col[order(col_ord)]

pdf(paste(dato,"G_HM_selcetedgenes_incl51.pdf",sep="_"), width = 4, height = 8)
DoHeatmap(Av_visu, features = HM_E_genes, group.colors = clus_v4_col[order(col_ord)], draw.lines = F)+
  scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()

Av_visu_sub <- subset(Av_visu, idents = c("50_DC1","53","19","cDC1","32","36","cDC2","50_DC3","35","29","cDC3"))
col_ord_2 <- factor(names(clus_v4_col),levels = c("50_DC1","53","19","cDC1","32","36","cDC2","50_DC3","35","29","cDC3"))

pdf(paste(dato,"G_HM_selcetedgenes_comittedandmat.pdf",sep="_"), width = 4, height = 8)
DoHeatmap(Av_visu_sub, features = HM_E_genes, group.colors = clus_v4_col[order(col_ord_2)], draw.lines = F)+
  scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()


#### H) TF expression across branches ####
# first update clustering to split for the above found commitments
DCtraj@meta.data$tSP_clus_comp_v4 <- as.character(DCtraj@meta.data$tSP_clus_comp_v3)
cl50 <- rownames(visu_obj@meta.data[visu_obj@meta.data$tSP_clus_comp_v3=="50",])
cl50_DC3 <- cl50[which(DCtraj@assays$integrated@data["FCER1A",cl50]>1)]
cl50_DC1 <- cl50[which(DCtraj@assays$integrated@data["CLEC9A",cl50]>1)]

DCtraj@meta.data[cl50_DC1,]$tSP_clus_comp_v4 <- "50_DC1"
DCtraj@meta.data[cl50_DC3,]$tSP_clus_comp_v4 <- "50_DC3"

#save colors as in tSP_clus_comp_v2
clus_v4_col <- c("cDC1"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
                 "53"="#F8766D","19"="#FF61CC","32"="#7CAE00","36"="#00BE67","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","50_DC1"="#C77CFF","50_DC3"="#C77CFF","51"="#CD9600")

pdf(paste(dato,"tSP_DCtraj_compclust_v3_wlabels.pdf"), width=4, height=3.5)
DimPlot(DCtraj, group.by = "tSP_clus_comp_v4", reduction = "flat_3D")+scale_color_manual(values=clus_v4_col)
dev.off()

Idents(DCtraj) <- 'tSP_clus_comp_v5'
DC_subsets <- list(DC1 = c("51","50_DC1","53","19","cDC1"), DC2 = c("51","32","36","cDC2"), DC3 = c("51","50_DC3","35","29","cDC3"))
cols_incl <- c("tSP_clus_comp_v5","orig.ident","segment")

Tom_TF_list <- read.csv("/Volumes/SUND/Public/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata_v2/DC_graj_TFs_Tom.csv")
Tom_TF_list <- Tom_TF_list[1:8,11:14]

Tom_lims <- list(GROUP1=c(-1,1.2),GROUP2=c(-1.2,1.5), GROUP3=c(-1.2,1.6),GROUP4=c(-1.2,3.2),GROUP5=c(-5,3.5))
#Idents(DCtraj) <- DCtraj@meta.data$tSP_clus_comp_v4

### new list og genes to study without looking at Tom's
Tom_TF_list["GROUP5"] <- c("KIT","ITGAX","HLA-DRA","HLA-DRB1","CD14","CD207","","")
Tom_TF_list <- Tom_TF_list[5]

for (group in names(Tom_TF_list)){
  for (DC_sub in rev(names(DC_subsets))){
    print(DC_sub); print(DC_subsets[[DC_sub]])
    DCtraj_sub <- subset(DCtraj, idents = DC_subsets[[DC_sub]])
    # RESCALE
    #DCtraj_sub <- ScaleData(DCtraj_sub, vars.to.regress = c("nCount_RNA", "percent.mt","G2M.Score","S.Score","cc.exp"), do.scale = TRUE) 
    TF_tot_sub <- as.character(Tom_TF_list[[group]][Tom_TF_list[[group]]!=""])
    print(TF_tot_sub)
    # data frame for plotting
    DC_df <- cbind(DCtraj_sub@meta.data[,cols_incl],t(as.matrix(DCtraj_sub@assays$integrated@scale.data[TF_tot_sub,])))
    DC_df_av <- data.frame()
    clusters <- c()
    # avergaes of data frame
    for (clus in unique(DC_df$tSP_clus_comp_v5)){
       DC_df_av <- rbind(DC_df_av, c(colMeans(DC_df[DC_df$tSP_clus_comp_v5==clus,4:dim(DC_df)[2]])))
       colnames(DC_df_av) <- c(TF_tot_sub)
       clusters <- append(clusters,clus)
    }
    DC_df_av <- cbind(DC_df_av, tSP_clus_comp_v5=clusters)
    #colnames(DC_df_av) <- colnames(DC_df)[c(4:dim(DC_df)[2],1)]
    DC_df_av <- gather(DC_df_av, gene, expression, TF_tot_sub[1]:TF_tot_sub[length(TF_tot_sub)], factor_key=TRUE)
    DC_df_av$tSP_clus_comp_v5 <- factor(DC_df_av$tSP_clus_comp_v5, levels=DC_subsets[[DC_sub]])
    plot1 <- ggplot(DC_df_av, aes(x=tSP_clus_comp_v5, y=as.numeric(expression), color=gene))+geom_line(aes(group = gene), size=1.5)+
       theme_classic()+xlab("")+ylab("scaled gene expressioon")+ggtitle(paste(DC_sub,group))+
       ylim(Tom_lims[[group]])+
       theme(plot.title = element_text(hjust=0.5))
    pdf(paste(dato,"NotRescaled_DCtraj_TFexp_branch",DC_sub,group,".pdf",sep="_"), height = 5, width = 5)
    print(plot1)
    dev.off()
  }}





#### H) Gene expression (well known) across branches ####
HM_E_genes <- c("STMN1", "KIAA0101", "TUBA1B","MKI67", "HIST1H4C","LILRA4","GZMB","E2F1","CD38",
                "RUNX3","THBD","CLEC4C", "BATF3","BTLA","CADM1","XCR1","IRF8","ID2","CLEC9A",
                "CCR6","NR4A1","ATF3","SOX4","S100B","PRDM16","CXXC5","NGFR","CDH1","KIT","LTB",
                "IL22RA2","IRF4","PLAC8","CD207","AXL","CD1C","IL3RA","C1QA","CD14","CD163", "FCGR1A",
                "S100A9", "FCN1", "VCAN","HLA-DQA1","HLA-DRA")

Sign_list <- list(GROUPDC1=c("BATF3","IRF8","XCR1","CLEC9A","CADM1"),GROUPDC2=c("LTB","CD207","CD1C","IRF4"),GROUPDC3=c("CD163","CD14","S100A9","C1QA","MERTK"))

Sign_lims <- list(GROUPDC1=c(-1.1,2.5),GROUPDC2=c(-1.5,2.7), GROUPDC3=c(-0.85,1))
DC_subsets <- list(DC1 = c("51","50_DC1","35_DC1","53","19","cDC1"))#, DC2 = c("51","32","36","cDC2"), DC3 = c("51","50_DC3","35","29","cDC3"))

for (group in names(Sign_list)){
  for (DC_sub in rev(names(DC_subsets))){
    print(DC_sub); print(DC_subsets[[DC_sub]])
    DCtraj_sub <- subset(DCtraj, idents = DC_subsets[[DC_sub]])
    # RESCALE
    #DCtraj_sub <- ScaleData(DCtraj_sub, vars.to.regress = c("nCount_RNA", "percent.mt","G2M.Score","S.Score","cc.exp"), do.scale = TRUE) 
    TF_tot_sub <- as.character(Sign_list[[group]][Sign_list[[group]]!=""])
    print(TF_tot_sub)
    # data frame for plotting
    DC_df <- cbind(DCtraj_sub@meta.data[,cols_incl],t(as.matrix(DCtraj_sub@assays$integrated@scale.data[TF_tot_sub,])))
    DC_df_av <- data.frame()
    clusters <- c()
    # avergaes of data frame
    for (clus in unique(DC_df$tSP_clus_comp_v5)){
      DC_df_av <- rbind(DC_df_av, c(colMeans(DC_df[DC_df$tSP_clus_comp_v5==clus,4:dim(DC_df)[2]])))
      colnames(DC_df_av) <- c(TF_tot_sub)
      clusters <- append(clusters,clus)
    }
    DC_df_av <- cbind(DC_df_av, tSP_clus_comp_v5=clusters)
    #colnames(DC_df_av) <- colnames(DC_df)[c(4:dim(DC_df)[2],1)]
    DC_df_av <- gather(DC_df_av, gene, expression, TF_tot_sub[1]:TF_tot_sub[length(TF_tot_sub)], factor_key=TRUE)
    DC_df_av$tSP_clus_comp_v5 <- factor(DC_df_av$tSP_clus_comp_v5, levels=DC_subsets[[DC_sub]])
    plot1 <- ggplot(DC_df_av, aes(x=tSP_clus_comp_v5, y=as.numeric(expression), color=gene))+geom_line(aes(group = gene), size=1.5)+
      theme_classic()+xlab("")+ylab("scaled gene expressioon")+ggtitle(paste(DC_sub,group))+
      ylim(Sign_lims[[group]])+
      theme(plot.title = element_text(hjust=0.5))
    pdf(paste(dato,"v2w35DC1_NotRescaled_DCtraj_TFexp_branch",DC_sub,group,".pdf",sep="_"), height = 5, width = 5)
    print(plot1)
    dev.off()
  }}

### specific signalling pathways
KEGG <- read.csv("/Volumes/sund$/WA group/Tom and Line data/cLP_SILP_merged/R6//Rdata_v2/KEGG_2021_Human.csv", na.strings = "")

Idents(DCtraj) <- 'tSP_clus_comp_v4'
clus_v4_col <- c("cDC1"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
                 "53"="#F8766D","19"="#FF61CC","32"="#7CAE00","36"="#00BE67","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","50_DC1"="#C77CFF","50_DC3"="#C77CFF","51"="#CD9600")
clus_ord <- rev(c("51","50_DC1","53","19","cDC1","32","36","cDC2","50_DC3","35","29","cDC3"))
Idents(DCtraj) <- factor(Idents(DCtraj), levels=rev(clus_ord))

Av_visu <- AverageExpression(subset(DCtraj, idents = names(clus_v4_col)[names(clus_v4_col) %in% clus_ord]), return.seurat = T)
Idents(Av_visu) <- factor(Idents(Av_visu), levels = rev(clus_ord))

col_ord <- factor(names(clus_v4_col),levels =rev(clus_ord)) #clus_v3_col[order(col_ord)]

## LTB - NF-kappa signalling
NFkappa <- KEGG[ "NF.kappa.B.signaling.pathway"][!is.na(KEGG[ "NF.kappa.B.signaling.pathway"])]

pdf(paste(dato,"G_HM_selcetedgenes_incl51.pdf",sep="_"), width = 4, height = 8)
DoHeatmap(Av_visu, features = NFkappa, group.colors = clus_v4_col[order(col_ord)], draw.lines = F, assay = "RNA", slot = "data")+
  scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()


## NOTCH2 signalling
NOTCH <- KEGG[ "Notch.signaling.pathway" ][!is.na(KEGG[ "Notch.signaling.pathway" ])]

pdf(paste(dato,"G_HM_selcetedgenes_incl51.pdf",sep="_"), width = 4, height = 8)
DoHeatmap(Av_visu, features = NOTCH, group.colors = clus_v4_col[order(col_ord)], draw.lines = F, assay = "RNA", slot = "data")+
  scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()

#### New H) - AUG 2023 - KIT, ITGAX and HLA score exp exp ####
# first colouring for DCs
clus_v3_col <- c("cDC1"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
                 "53"="#F8766D","19"="#FF61CC","32"="#7CAE00","36"="#00BE67",
                 "35"="#00BFC4","29"="#00A9FF",
                 "50"="#C77CFF","51"="#CD9600")
Idents(DCtraj) <- factor(DCtraj@meta.data$tSP_clus_comp_v3, levels = 
                           c("51","50","53","19","32","36","35","29","cDC1","cDC2","amb","cDC3"))
#check
DimPlot(DCtraj, group.by = "tSP_clus_comp_v3", reduction = "flat_3D")+scale_color_manual(values = clus_v3_col)

# HLA module
HLA_genes <- list(rownames(DCtraj@assays$integrated@data)[startsWith(rownames(DCtraj@assays$integrated@data),"HLA")])
DCtraj <- AddModuleScore(DCtraj, features = HLA_genes, ctrl = 5, name = 'HLA.score')

# Plot and save
pdf(paste(dato,"F4_Dctraj_Vln","CD11cgene.pdf",sep="_"), width=4, height=4)
VlnPlot(subset(DCtraj, idents = c("51","50","53","19","32","36","35","29")), 
        features = c("ITGAX"), pt.size = 0)+NoLegend()+
  scale_fill_manual(values=clus_v3_col)
dev.off()
pdf(paste(dato,"F4_Dctraj_Vln","KITgene.pdf",sep="_"), width=4, height=4)
VlnPlot(subset(DCtraj, idents = c("51","50","53","19","32","36","35","29","cDC2")), 
        features = c("KIT"), pt.size = 0)+NoLegend()+
  scale_fill_manual(values=clus_v3_col)
dev.off()
pdf(paste(dato,"F4_Dctraj_Vln","HLAscore.pdf",sep="_"), width=4, height=4)
VlnPlot(subset(DCtraj, idents = c("51","50","53","19","32","36","35","29")), 
        features = c("HLA.score1"), pt.size = 0)+NoLegend()+
  scale_fill_manual(values=clus_v3_col)
dev.off()

head()

#### I) CITE-seq expression in precursors ####
rownames(DCtraj@assays$CITE@data)[4:13]
clus_ord <- rev(c("51","50","50_DC1","53","19","cDC1","32","36","cDC2","50_DC3","35","29","cDC3"))

DCtraj_sub <- subset(DCtraj, cells=rownames(DCtraj@meta.data[DCtraj@meta.data$cite2=="YES" & DCtraj@meta.data$segment=="cLP" & DCtraj@meta.data$tSP_clus_comp_v4!="amb",]))
Idents(DCtraj_sub) <- factor(Idents(DCtraj_sub), levels=clus_ord)
RidgePlot(DCtraj_sub, features="cite_CD14")+scale_fill_manual(values=clus_v4_col)

## In individual branches
DC_subsets <- list(DC1 = c("50_DC1","53","19","cDC1"), DC2 = c("32","36","cDC2"), DC3 = c("35","29","cDC3"))
genes <- rownames(DCtraj@assays$CITE@data)[c(4,6:12)]
genes_lim <- list("CD14"=c(-2,5),"CD1c"=c(-1,5),"CD55"=c(-2,5),"CD5"=c(-1,12),"CD207"=c(-1,12),"CD206"=c(-2,5),"CD11a"=c(-2,4),"CD103"=c(-1,10))
#gene_name="cite_CD14"
for (gene_name in genes){
  gennam <- gene_name
  gene_name <- paste("cite",gene_name,sep="_")
for (DC_sub in rev(names(DC_subsets))){
    print(DC_sub); print(DC_subsets[[DC_sub]])
    cell_sub <- rownames(DCtraj@meta.data[DCtraj@meta.data$cite=="YES" & DCtraj@meta.data$segment=="cLP" & DCtraj@meta.data$tSP_clus_comp_v4 %in% DC_subsets[[DC_sub]],])
    DCtraj_sub <- subset(DCtraj, cells = cell_sub)
    Idents(DCtraj_sub) <- factor(Idents(DCtraj_sub), levels=DC_subsets[[DC_sub]])
    plot1 <- RidgePlot(DCtraj_sub, features = gene_name)+scale_fill_manual(values=clus_v4_col)+NoLegend()+xlim(genes_lim[[gennam]])
    pdf(paste(dato,"RidgePlot",DC_sub,gene_name,".pdf",sep="_"), height = 3, width = 3)
    print(plot1)
    dev.off()
}}

genes <- rownames(DCtraj@assays$CITE@data)[c(5,13)]
genes_lim <- list("CD209"=c(-2.5,5),"CD11c"=c(-2,3.5))
#gene_name="cite_CD14"
for (gene_name in genes){
  gennam <- gene_name
  gene_name <- paste("cite",gene_name,sep="_")
  for (DC_sub in rev(names(DC_subsets))){
    print(DC_sub); print(DC_subsets[[DC_sub]])
    cell_sub <- rownames(DCtraj@meta.data[DCtraj@meta.data$cite2=="YES" & DCtraj@meta.data$segment=="cLP" & DCtraj@meta.data$tSP_clus_comp_v4 %in% DC_subsets[[DC_sub]],])
    DCtraj_sub <- subset(DCtraj, cells = cell_sub)
    Idents(DCtraj_sub) <- factor(Idents(DCtraj_sub), levels=DC_subsets[[DC_sub]])
    plot1 <- RidgePlot(DCtraj_sub, features = gene_name)+scale_fill_manual(values=clus_v4_col)+NoLegend()+xlim(genes_lim[[gennam]])
    pdf(paste(dato,"RidgePlot",DC_sub,gene_name,".pdf",sep="_"), height = 3, width = 3)
    print(plot1)
    dev.off()
  }}






#### Supplementary Figure ####
#### A) Velo justifications ####
Idents(DCtraj) <- 'tSP_clus_comp'
DEGs <- FindAllMarkers(subset(DCtraj, idents = c("cDC1A","cDC2","cDC3")), #c("19","29","36") c("cDC1A","cDC3","cDC2") c("36","53","35")
                       only.pos = T, #only including UPregulated DEGs per cluster
                       test.use = "wilcox", #test to use
                       min.pct = 0.25, # genes included in results should be expressed in atleast 25% of the cluster
                       logfc.threshold = 0.25)

DEGs <- DEGs[DEGs$p_val_adj<0.05,]
dim(DEGs)
top50 <- DEGs %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

mat_genes <- top50$gene

Idents(DCtraj) <- 'tSP_clus_comp'
clus_col <- c("cDC1A"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
              "53"="#F8766D","19"="#FF61CC","32"="#7CAE00","36"="#00BE67","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","51"="#CD9600")
clus_ord <- c("51","50","32","53","19","cDC1A","36","cDC2","35","29","cDC3")
Idents(DCtraj) <- factor(Idents(DCtraj), levels=rev(clus_ord))

Av_visu <- AverageExpression(subset(DCtraj, idents = names(clus_col)[names(clus_col) %in% clus_ord]), return.seurat = T, assay= "integrated")
Idents(Av_visu) <- factor(Idents(Av_visu), levels = clus_ord)

pdf(paste(dato,"SF4_velojustictation_top50matDEGs.pdf",sep="_"), width = 4, height = 6)
DoHeatmap(Av_visu, features = mat_genes, group.colors = clus_col[clus_ord], draw.lines = F)+
  scale_fill_gradientn(colours = mycols, na.value = "white")+theme(axis.text.y = element_blank())
dev.off()

pdf(paste(dato,"SF4_velojustictation_top50matDEGs_withgenenames.pdf",sep="_"), width = 4, height = 10)
DoHeatmap(Av_visu, features = mat_genes, group.colors = clus_col[clus_ord], draw.lines = F)+
  scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()

#### B) Velo split in SI and LI ####

interest_SI <- c()

nmat.int_SI <- nmat.var[,interest_SI]
emat.int_SI <- emat.var[,interest_SI]

#Estimation of RNA velwith k=20 cells kNN pooling and top/bottom, 2% quantiles for gamma fit
fit.quantile <- 0.02
#try other values for deltaT, try kGenes insetad of kCells and kCells=5 
rvel.cd <- gene.relative.velocity.estimates(emat.int,nmat.int,deltaT=1,kCells=440, cell.dist=cell.dist,fit.quantile=fit.quantile, zero.offset = T, n.cores = 4)

#visualization of RNA vel 
#play around with this
#e.g. suggested kCells=5
#worked well but slow
pdf(paste(dato,"DCtraj_velocytoCellInterest_PearsonEmb_tSPv1001.pdf",sep="_"), height = 5, width = 5)
show.velocity.on.embedding.cor(emb,rvel.cd,n=400,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.70),cex=1.5,n.cores = 4,
                               arrow.scale=2,show.grid.flow=TRUE,min.grid.cell.mass=0.1,grid.n=25,arrow.lwd=0.8,do.par=T,cell.border.alpha = 0)
dev.off()

## Velo with extra precursor clusters
clus_v3.5_col <- c("cDC1"="grey","cDC2"="grey","amb"="grey","cDC3"="grey",
                   "53"="#F8766D","19"="#FF61CC","43"="darksalmon","30"="rosybrown","32"="#7CAE00","36"="#00BE67","33"="palegreen3","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","51"="#CD9600")
DCtraj@meta.data$tSP_clus_comp_v3.5 <- as.character(DCtraj@meta.data$tSP_clus_comp_v3)
DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3 %in% c(43,30,33),]$tSP_clus_comp_v3.5 <- as.character(DCtraj@meta.data[DCtraj@meta.data$integrated_snn_res.3 %in% c(43,30,33),]$integrated_snn_res.3)

DimPlot(DCtraj, group.by = "tSP_clus_comp_v3.5", reduction = "flat_3D", label=T, label.size=6, pt.size = 1)+scale_color_manual(values=clus_v3.5_col)+NoLegend()

visu_obj <- DCtraj

Idents(visu_obj) <- visu_obj@meta.data$tSP_clus_comp_v3.5
#DC_traj_ord <- c(4,3,7,1,5,2,0,6)
Idents(visu_obj) <- factor(Idents(visu_obj), levels = names(clus_v3.5_col))
cluster.label <- Idents(visu_obj) #numbers now 0-5 - 6 clusters - going to change with new clustering
#specific colours from clustering
#add for each cell w. pagoda
cell.colors <- pagoda2:::fac2col(cluster.label,level.colors = clus_v3.5_col)

pdf(paste(dato,"DCtraj_velocytoCellInterest_PearsonEmb_tSPclus3.5_v1.pdf",sep="_"), height = 5, width = 5)
show.velocity.on.embedding.cor(emb,rvel.cd,n=400,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.70),cex=1.5,n.cores = 4,
                               arrow.scale=2,show.grid.flow=TRUE,min.grid.cell.mass=0.1,grid.n=25,arrow.lwd=0.8,do.par=T,cell.border.alpha = 0,
                               axes = FALSE)
dev.off()


#### top 5 DEGs before and ####
#### HM D) - OBS v3.1.1 required for HM ####
DEGs <- FindAllMarkers(subset(DCtraj, idents = c("cDC1","cDC2","cDC3")), assay = "RNA",
                       base = exp(1), only.pos = T, logfc.threshold = 0.2)
DEGs_top20 <- DEGs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DEGs_top20 <- as.character(DEGs_top20$gene[c(1:20,41:60,21:40)])

Idents(DCtraj) <- "tSP_clus_comp_v3"
Av_visu <- AverageExpression(DCtraj, return.seurat = T)
Idents(Av_visu) <- factor(Idents(Av_visu), levels = c("51","50","32","53","19","cDC1","36","cDC2","35","29","cDC3","amb"))
Av_visu@meta.data$orig.ident <- rownames(Av_visu@meta.data)

#save colors as in tSP_clus_comp_v2
Av_visu <- subset(Av_visu, cells = colnames(Av_visu)[colnames(Av_visu)!="amb"])
clus_v3_col <- c("cDC1"="grey","cDC2"="grey","cDC3"="grey",
                 "53"="#F8766D","19"="#FF61CC","32"="#7CAE00","36"="#00BE67","35"="#00BFC4","29"="#00A9FF","50"="#C77CFF","51"="#CD9600")
col_ord <- factor(names(clus_v3_col),levels = c("51","50","32","53","19","cDC1","36","cDC2","35","29","cDC3"))

pdf(paste(dato,"D_HM_selcetedgenes_allgroups.pdf",sep="_"), width = 4, height = 8)
DoHeatmap(Av_visu, features = DEGs_top20, assay="integrated", label = F, draw.lines = F, group.colors = clus_v3_col[order(col_ord)])+
  scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()

#### D) - lit. preDC exp. across branch - 8/8-23 #### 
#Idents(DCtraj) <- 'tSP_clus_comp_v5'
DC_subsets <- list(DC1 = c("51","50_DC1","35_DC1","53","19","cDC1"), DC2 = c("51","32","36","cDC2"), DC3 = c("51","50_DC3","35","29","cDC3"))
cols_incl <- c("tSP_clus_comp_v5","orig.ident","segment")

Lit_list_org <- read.csv("/Users/lwc/Desktop/temp/2307_preDCmarkers_Literature.csv")
Lit_list_org <- Lit_list_org[,c("gene","group")]

#transform to proper format list for below plotting function
Lit_list <- list()
for (group in unique(Lit_list_org$group)){
  Lit_sub_list <- Lit_list_org[Lit_list_org$group==group,]
  Lit_list[[group]] <- Lit_sub_list$gene
}

#check all in integrated / scaled data
which(duplicated(unlist(Lit_list))) #no duplicates
unlist(Lit_list)[!unlist(Lit_list) %in% rownames(DCtraj@assays$integrated@scale.data)] #yes some only in RNA others not there at all

Lit_lims <- list(group1a=c(-1,1.5),
                 group1b=c(-1,3.5),
                 group2=c(-1.5,2.5), 
                 group3=c(-2,1.5),
                 group4=c(-2.5,2),
                 group5=c(-2,1),
                 group6=c(-5,5))

### new list og genes to study without looking at Tom's
TF_tot_sub_lack <- c() #empty for those not in integrated/scaled data
for (group in names(Lit_list)){
  for (DC_sub in rev(names(DC_subsets))){
    print(DC_sub); print(DC_subsets[[DC_sub]])
    DCtraj_sub <- subset(DCtraj, idents = DC_subsets[[DC_sub]])
    # RESCALE
    #DCtraj_sub <- ScaleData(DCtraj_sub, vars.to.regress = c("nCount_RNA", "percent.mt","G2M.Score","S.Score","cc.exp"), do.scale = TRUE) 
    TF_tot_sub <- as.character(Lit_list[[group]][Lit_list[[group]]!=""])
    TF_tot_sub_lack <- c( TF_tot_sub_lack,TF_tot_sub[!TF_tot_sub %in% rownames(DCtraj_sub@assays$integrated@scale.data)])
    TF_tot_sub <- TF_tot_sub[TF_tot_sub %in% rownames(DCtraj_sub@assays$integrated@scale.data)]
    print(TF_tot_sub)
    # data frame for plotting
    DC_df <- cbind(DCtraj_sub@meta.data[,cols_incl],t(as.matrix(DCtraj_sub@assays$integrated@scale.data[TF_tot_sub,])))
    DC_df_av <- data.frame()
    clusters <- c()
    # avergaes of data frame
    for (clus in unique(DC_df$tSP_clus_comp_v5)){
      DC_df_av <- rbind(DC_df_av, c(colMeans(DC_df[DC_df$tSP_clus_comp_v5==clus,4:dim(DC_df)[2]])))
      colnames(DC_df_av) <- c(TF_tot_sub)
      clusters <- append(clusters,clus)
    }
    DC_df_av <- cbind(DC_df_av, tSP_clus_comp_v5=clusters)
    #colnames(DC_df_av) <- colnames(DC_df)[c(4:dim(DC_df)[2],1)]
    DC_df_av <- gather(DC_df_av, gene, expression, TF_tot_sub[1]:TF_tot_sub[length(TF_tot_sub)], factor_key=TRUE)
    DC_df_av$tSP_clus_comp_v5 <- factor(DC_df_av$tSP_clus_comp_v5, levels=DC_subsets[[DC_sub]])
    plot1 <- ggplot(DC_df_av, aes(x=tSP_clus_comp_v5, y=as.numeric(expression), color=gene))+geom_line(aes(group = gene), size=1.5)+
      theme_classic()+xlab("")+ylab("scaled gene expressioon")+ggtitle(paste(DC_sub,group))+
      ylim(Lit_lims[[group]])+
      theme(plot.title = element_text(hjust=0.5))
    pdf(paste(dato,"NotRescaled_DCtraj_Litexp_branch",DC_sub,group,".pdf",sep="_"), height = 5, width = 5)
    print(plot1)
    dev.off()
  }}

# Now Vln plot of those with only RNA info
TF_tot_sub_lack <- unique(TF_tot_sub_lack)

Idents(DCtraj) <- factor(DCtraj@meta.data$tSP_clus_comp_v5, 
                         levels = c("51","50",   "50_DC1","35_DC1","53","19","cDC1",
                                    "32","36","cDC2",
                                    "50_DC3", "35", "29","cDC3",
                                    "amb"))
DCtraj_sub

for (gene in TF_tot_sub_lack[TF_tot_sub_lack %in% rownames(DCtraj@assays$RNA@data)]){
plot1 <- VlnPlot(DCtraj_sub, features = gene, assay = "RNA", pt.size = 0)
pdf(paste(dato,"DCtraj_Litexp_VlnBranch_RNAavail",gene,".pdf",sep="_"), height = 4, width = 6)
print(plot1)
dev.off()
}

