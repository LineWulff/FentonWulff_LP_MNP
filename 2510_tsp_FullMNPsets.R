###### Libraries ######
library(gplots)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(viridis)
library(scales)
library(Seurat)
library(ggdendro)
library(plotly)
library(ggrastr)
library(plot3D)
library(rgl)
library(umap)
library(htmlwidgets)

### Clear and read in data ###
rm(list=ls())

### variables
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
# project
project <- "MNP_LP"

outdir <- "/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/tspace_plots"

#### data read in and variables ####
MNP <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata_v2/R3_CD14CD1Cprol_MNPLP.rds")
ggplot(MNP@reductions$umap@cell.embeddings, aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=MNP@meta.data$integrated_snn_res.0.5))
MNP <- UpdateSeuratObject(MNP)
subs1 <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/tspace_output/R3_CD14CD1Cprol_MNPLP_tspacefile.rds")

head(subs1$ts_file)

MNP@meta.data$integrated_snn_res.2.8

ord <- c(12,24,22,4,33,16,19,18,35,17,28,0,10,15,2,5,29,21,3,13,37,34,36,38,26,7,8,30,32,27,11,20,9,1,14,6,31,25,23)
ord <- as.character(ord)
names(ord) <- c(1:length(ord))
MNP@meta.data$res2.8_ord <- factor(names(ord)[MNP@meta.data$integrated_snn_res.2.8], levels = c(1:39))


Idents(MNP) <- 'res2.8_ord'

MNP@meta.data$res2.8_ord

#### Prepare for visualization
visu <- subs1$ts_file
visu$res2.8_ord <- MNP@meta.data[rownames(visu),]$res2.8_ord
visu$lineage <- NA
visu[visu$res2.8_ord %in% as.character(seq(1,20)),]$lineage <- "mono/mac"
visu[visu$res2.8_ord %in% as.character(seq(21,39)),]$lineage <- "cDC2/3"

pdf(paste0(outdir,"/",dato,"CD1CCD14tsp_2D_Fig1colouring.pdf"), height = 5, width = 6)
ggplot(visu, aes(x=umap1, y=umap2, color=res2.8_ord))+
  geom_point_rast(size=0.5)+
  xlab("tUMAP1")+ylab("tUMAP2")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

pdf(paste0(outdir,"/",dato,"CD1CCD14tsp_2D_LineageColouring.pdf"), height = 5, width = 6)
ggplot(visu, aes(x=umap1, y=umap2, color=lineage))+
  geom_point_rast(size=0.5)+
  xlab("tUMAP1")+ylab("tUMAP2")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()


visu$Lincol <-visu$lineage
lin_cols <- hue_pal()(length(unique(visu$Lincol)))
names(lin_cols) <- unique(visu$lineage)[order(unique(visu$lineage))]
visu$Lincol <- lin_cols[visu$lineage]

visu$LinCol <- MNP@meta.data[rownames(visu),]$res2.8_ord
res2.8_ord_cols <- hue_pal()(length(unique(visu$res2.8_ord)))
names(res2.8_ord_cols) <- unique(visu$res2.8_ord)[order(unique(visu$res2.8_ord))]
visu$res2.8_ord_col <- res2.8_ord_cols[visu$res2.8_ord]


seed <- 1111
config_tspace <- umap::umap.defaults;config_tspace$n_neighbors <- 7;config_tspace$min_dist <- 0.2;config_tspace$metric <- "manhattan"
# 2D vs 3D umap outputs
config_tspace$n_components <- 3
set.seed(seed)
umap_tspace <- umap::umap(as.matrix(visu[,c(4:18)]), config = config_tspace)
umap_out <- as.data.frame(umap_tspace$layout)
colnames(umap_out) <- paste0("umap", seq(1, ncol(umap_tspace$layout),1))
head(umap_out)

plot3d(x=umap_out[,"umap1"],y=umap_out[,"umap2"],z=umap_out[,"umap3"],
                col=visu$res2.8_ord_col,
                xlab = "tUMAP1",ylab = "tUMAP2",zlab = "tUMAP3")
writeWebGL(dir = outdir, filename = paste0(dato,"CD1CD14tsp_3D_Fig1Colouring.html"))

plot3d(x=umap_out[,"umap1"],y=umap_out[,"umap2"],z=umap_out[,"umap3"],
       col=visu$Lincol,
       xlab = "tUMAP1",ylab = "tUMAP2",zlab = "tUMAP3")
writeWebGL(dir = outdir, filename = paste0(dato,"CD1CD14tsp_3D_LineageColouring.html"))

#### All MNP version ####
load("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata/R2_MNPs_MNP_LP.Rdata")
MNP <- UpdateSeuratObject(MNP_LP)
CD14CD1Cmeta <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata_v2/R3_CD14CD1Cprol_MNPLP_METAONLY_AUG22.rds")
monomac <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata_v2/monomac_traj_versionApril2022.rds")
monomac <- monomac@meta.data

MNP@meta.data$F2F5 <- NA
MNP@meta.data[rownames(CD14CD1Cmeta),]$F2F5 <- CD14CD1Cmeta$semisupDC_v5
MNP@meta.data[MNP@meta.data$integrated_snn_res.0.1==3,]$F2F5 <- "cDC1"
MNP@meta.data[MNP@meta.data$integrated_snn_res.0.1==6,]$F2F5 <- "CD16+ mono"
MNP@meta.data[MNP@meta.data$integrated_snn_res.0.1==5,]$F2F5 <- "pDC"
MNP@meta.data[rownames(monomac),]$F2F5 <- as.character(monomac$tSP_clustering_F4)

DimPlot(MNP, group.by = "F2F5", label = T)

saveRDS(MNP,"R2_AllMNPs.rds")

subs1 <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/tspace_output/R2_AllMNPs_2Dtspacefile.rds")
visu <- subs1$ts_file

visu$F2F5 <- MNP@meta.data[rownames(visu),]$F2F5

pdf(paste0(outdir,"/",dato,"FullMNPtsp_2D_indfigcolouring.pdf"), height = 5, width = 6)
ggplot(visu, aes(x=umap1, y=umap2, color=F2F5))+
  geom_point_rast(size=0.5)+
  xlab("tUMAP1")+ylab("tUMAP2")+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

MNP@reductions$tumap <- MNP@reductions$umap
MNP@reductions$tumap@cell.embeddings <- as.matrix(visu[Cells(MNP),c("umap1","umap2")])
MNP@meta.data <- cbind(MNP@meta.data, visu[Cells(MNP),c("umap1","umap2")])

subs1 <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/tspace_output/R2_AllMNPs_3Dtspacefile.rds")
visu <- subs1$ts_file

visu$F2F5 <- MNP@meta.data[rownames(visu),]$F2F5
F2F5_cols <- hue_pal()(length(unique(visu$F2F5)))
names(F2F5_cols) <- unique(visu$F2F5)[order(unique(visu$F2F5))]
visu$F2F5_col <- F2F5_cols[visu$F2F5]

seed <- 1111
config_tspace <- umap::umap.defaults;config_tspace$n_neighbors <- 7;config_tspace$min_dist <- 0.2;config_tspace$metric <- "manhattan"
# 2D vs 3D umap outputs
config_tspace$n_components <- 3
set.seed(seed)
umap_tspace <- umap::umap(as.matrix(visu[,c(4:18)]), config = config_tspace)
umap_out <- as.data.frame(umap_tspace$layout)
colnames(umap_out) <- paste0("umap", seq(1, ncol(umap_tspace$layout),1))
head(umap_out)

plot3d(x=umap_out[,"umap1"],y=umap_out[,"umap2"],z=umap_out[,"umap3"],
       col=visu$F2F5_col,
       xlab = "tUMAP1",ylab = "tUMAP2",zlab = "tUMAP3")
writeWebGL(dir = outdir, filename = paste0(dato,"FullMNPtsp_3D_indfigcolouring.html"))

## remove CD16+ and pDC
dim(visu) #28758
visu <- visu[!visu$F2F5 %in% c("CD16+ mono","pDC","mono/mac"),]
dim(visu) #27918

seed <- 1111
config_tspace <- umap::umap.defaults;config_tspace$n_neighbors <- 7;config_tspace$min_dist <- 0.2;config_tspace$metric <- "manhattan"
# 2D vs 3D umap outputs
config_tspace$n_components <- 3
set.seed(seed)
umap_tspace <- umap::umap(as.matrix(visu[,c(4:18)]), config = config_tspace)
umap_out <- as.data.frame(umap_tspace$layout)
colnames(umap_out) <- paste0("umap", seq(1, ncol(umap_tspace$layout),1))
head(umap_out)

## extract 3D UMAP and cluster labels
plot.data <- umap_out
plot.data$label <- paste(rownames(plot.data))
plot.data$F2F5 <- visu$F2F5

## color palette for clusters
cluster_colors <- unique(visu$F2F5_col)
names(cluster_colors) <- unique(visu$F2F5)

## 2Ds
pdf(paste0(outdir,"/",dato,"CD14CD1CcDC1tsp_2D_tUMAP1vstUMAP2_indfigcolouring.pdf"), height = 5, width = 6)
ggplot(visu, aes(x=umap_out$umap1, y=umap_out$umap2, color=F2F5))+
  geom_point_rast(size=0.5)+
  xlab("tUMAP1")+ylab("tUMAP2")+
  scale_color_manual(values = cluster_colors)+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

pdf(paste0(outdir,"/",dato,"CD14CD1CcDC1tsp_2D_tUMAP1vstUMAP3_indfigcolouring.pdf"), height = 5, width = 6)
ggplot(visu, aes(x=umap_out$umap1, y=umap_out$umap3, color=F2F5))+
  geom_point_rast(size=0.5)+
  xlab("tUMAP1")+ylab("tUMAP3")+
  scale_color_manual(values = cluster_colors)+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

pdf(paste0(outdir,"/",dato,"CD14CD1CcDC1tsp_2D_tUMAP2vstUMAP3_indfigcolouring.pdf"), height = 5, width = 6)
ggplot(visu, aes(x=umap_out$umap2, y=umap_out$umap3, color=F2F5))+
  geom_point_rast(size=0.5)+
  xlab("tUMAP2")+ylab("tUMAP3")+
  scale_color_manual(values = cluster_colors)+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

## create interactive 3D scatter plot
fig <- plot_ly(data = plot.data, 
               x = ~umap1, y = ~umap2, z = ~umap3,
               color = ~F2F5,
               jitter = T,
               colors = cluster_colors,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 4, width = 1), 
               text = ~label, 
               hoverinfo = "text") %>%
  layout(scene = list(camera = list(
    eye = list(
      x = 1.25,
      y = 1.25,
      z = 1.25
    ),
    center = list(x = 0,
                  y = 0,
                  z = 0)
  ))) %>%
  onRender("
      function(el, x){
  var id = el.getAttribute('id');
  var gd = document.getElementById(id);
  Plotly.update(id).then(attach);
  function attach() {
    var cnt = 0;
    
    function run() {
      rotate('scene', Math.PI / 1000);
      requestAnimationFrame(run);
    } 
    run();
    
    function rotate(id, angle) {
      var eye0 = gd.layout[id].camera.eye
      var rtz = xyz2rtz(eye0);
      rtz.t += angle;
      
      var eye1 = rtz2xyz(rtz);
      Plotly.relayout(gd, id + '.camera.eye', eye1)
    }
    
    function xyz2rtz(xyz) {
      return {
        r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
        t: Math.atan2(xyz.y, xyz.x),
        z: xyz.z
      };
    }
    
    function rtz2xyz(rtz) {
      return {
        x: rtz.r * Math.cos(rtz.t),
        y: rtz.r * Math.sin(rtz.t),
        z: rtz.z
      };
    }
  };
}
    ") %>% layout(legend = list(itemsizing = "constant")) # this one sizes the legend markers according to markers in plot

## save as standalone HTML
saveWidget(partial_bundle(fig), file = paste0(outdir,dato,"CD14CD1CcDC1tsp_3D.HTML"), selfcontained = TRUE)
  
### 2D
seed <- 1111
config_tspace <- umap::umap.defaults;config_tspace$n_neighbors <- 7;config_tspace$min_dist <- 0.2;config_tspace$metric <- "manhattan"
# 2D vs 3D umap outputs
config_tspace$n_components <- 2
set.seed(seed)
umap_tspace <- umap::umap(as.matrix(visu[,c(4:18)]), config = config_tspace)
umap_out <- as.data.frame(umap_tspace$layout)
colnames(umap_out) <- paste0("umap", seq(1, ncol(umap_tspace$layout),1))
head(umap_out)

pdf(paste0(outdir,"/",dato,"CD14CD1CcDC1tsp_2D_indfigcolouring.pdf"), height = 5, width = 6)
ggplot(visu, aes(x=umap_out$umap1, y=umap_out$umap2, color=F2F5))+
  geom_point_rast(size=0.5)+
  xlab("tUMAP1")+ylab("tUMAP2")+
  scale_color_manual(values = cluster_colors)+
  theme_classic()+
  theme(axis.ticks = element_blank(),axis.text = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

