###############################################################################################
############################ Analysis of CITE-seq data on MNP LP ##############################
###############################################################################################
setwd("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/clustering/CD14CD1C_prol")

rm(list=ls())

#### Read in libraries ####
library(gplots)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(Rmagic)
library(ggplot2)
library(viridis)
library(clustree)
library(ccRemover)
library(rgl)
library(dsb)
library(Seurat)

#### Variables to use ####
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
# project
project <- "MNP_LP_cite-dsb-norm"

#### Read in data ####
load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R1_AllCells_MNP_LP.Rdata")
#MNP_LP <- readRDS("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R1_AllCells_MNP_LP.rds")
DimPlot(MNP_LP, group.by = "patient")
#DimPlot(MNP_LP, group.by = "integrated_snn_res.2.8", label=T)
unique(MNP_LP@meta.data$orig.ident)

load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_5/analysis/Rdata/Pat5AllCells_CITEseq_dsbnorm.Rdata")
normalized_matrix_pat5cLP <- normalized_matrix
colnames(normalized_matrix_pat5cLP)[1:5]
rownames(MNP_LP@meta.data[MNP_LP@meta.data$orig.ident=="cLP_pat5",])[1:5]
colnames(normalized_matrix_pat5cLP) <- paste0("P5_cLP_",colnames(normalized_matrix_pat5cLP))

load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_6/analysis/Rdata/Pat6AllCells_CITEseq_dsbnorm.Rdata")
normalized_matrix_pat6cLP <- normalized_matrix
colnames(normalized_matrix_pat6cLP)[1:5]
rownames(MNP_LP@meta.data[MNP_LP@meta.data$orig.ident=="cLP_pat6",])[1:5]
colnames(normalized_matrix_pat6cLP) <- paste0("P6_cLP_",colnames(normalized_matrix_pat6cLP))

load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/18/analysis/Rdata/Pat4AllCells_CITEseq_dsbnorm.Rdata")
normalized_matrix_pat4cLP <- normalized_matrix
colnames(normalized_matrix_pat4cLP)[1:5]
rownames(MNP_LP@meta.data[MNP_LP@meta.data$orig.ident=="cLP_pat4",])[1:5]
colnames(normalized_matrix_pat4cLP) <- paste0("P4_cLP_",colnames(normalized_matrix_pat4cLP))

load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/analysis/Rdata/Pat4AllCells_CITEseq_dsbnorm.Rdata")
normalized_matrix_pat4SILP <- normalized_matrix
colnames(normalized_matrix_pat4SILP)[1:5]
rownames(MNP_LP@meta.data[MNP_LP@meta.data$orig.ident=="SILP_pat4",])[1:5]
colnames(normalized_matrix_pat4SILP) <- paste0("P4_SILP_",colnames(normalized_matrix_pat4SILP))

citenames <- c(colnames(normalized_matrix_pat5cLP),colnames(normalized_matrix_pat6cLP),colnames(normalized_matrix_pat4cLP),colnames(normalized_matrix_pat4SILP))

#rownames(MNP_LP@meta.data[MNP_LP@meta.data$patient=="P4",])
#colnames(normalized_matrix)
#colnames(normalized_matrix) <- paste0("P4_cLP_",colnames(normalized_matrix))

Inc_cite <- rownames(MNP_LP@meta.data)
Not_NA <- citenames[citenames %in% rownames(MNP_LP@meta.data)]
Tobe_NA <- Inc_cite[!Inc_cite %in% Not_NA]

## Samples without cite-seq data
notCITE_mat <- rep(NA, length(rownames(normalized_matrix_pat4cLP))*length(Tobe_NA))
dim(notCITE_mat) <- c(length(rownames(normalized_matrix)),length(Tobe_NA))
colnames(notCITE_mat) <- Tobe_NA
rownames(notCITE_mat) <- rownames(normalized_matrix_pat4cLP)

## the patient without isotypes and CD209 and CD11c, pat5 cLP
notCITE_mat_pat5 <- rep(NA, 5*length(colnames(normalized_matrix_pat5cLP)))
dim(notCITE_mat_pat5) <- c(5,length(colnames(normalized_matrix_pat5cLP)))
colnames(notCITE_mat_pat5) <- colnames(normalized_matrix_pat5cLP)
rownames(notCITE_mat_pat5) <- rownames(normalized_matrix_pat4cLP)[!rownames(normalized_matrix_pat4cLP) %in% rownames(normalized_matrix_pat5cLP)]
normalized_matrix_pat5cLP <- rbind(normalized_matrix_pat5cLP, notCITE_mat_pat5)
normalized_matrix_pat5cLP <- normalized_matrix_pat5cLP[rownames(normalized_matrix_pat4cLP),]


## Bind all citeseq data into one matrix
normalized_matrix <- normalized_matrix_pat4cLP %>% cbind(normalized_matrix_pat4SILP) %>%
  cbind(normalized_matrix_pat6cLP) %>% cbind(normalized_matrix_pat5cLP)
dim(normalized_matrix) #18136 cells


## Bind both citeseq and NA data together and choose only cells existing in MNP_LP object
pat04_CITE_wNA <- cbind(normalized_matrix,notCITE_mat)
dim(pat04_CITE_wNA) #42506 cells
saveRDS(pat04_CITE_wNA, file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/Citeseq_DSBNormMat_AllCells.rds")
pat04_CITE_wNA <- pat04_CITE_wNA[,Inc_cite]
#check that there are the same number of cells in the objects
dim(pat04_CITE_wNA)
MNP_LP


MNP_LP[["CITE"]] <- CreateAssayObject(data = pat04_CITE_wNA)
MNP_LP <- ScaleData(MNP_LP, assay = "CITE")

MNP_LP@meta.data$cite <- "NO"
MNP_LP@meta.data[Not_NA,]$cite <- "YES"
MNP_LP@meta.data$cite2 <- "NO"
MNP_LP@meta.data[Not_NA[!Not_NA %in% colnames(normalized_matrix_pat5cLP)],]$cite2 <- "YES"


Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.0.5
DoHeatmap(subset(MNP_LP, cells = Not_NA), features = rownames(MNP_LP@assays$CITE@data)[1:13], assay = "CITE", slot = "data", 
          raster = F, lines.width = 20, disp.max = 2)

FeaturePlot(MNP_LP, features = c("cite_CD14","cite_CD207","cite_CD55","cite_CD1c"), cols=mycols_b, cells = rownames(MNP_LP@meta.data[MNP_LP@meta.data$cite=="YES",]))
FeaturePlot(MNP_LP, features = c("cite_CD14","cite_CD207","cite_CD11a","cite_CD209"), cols=mycols_b, cells = rownames(MNP_LP@meta.data[MNP_LP@meta.data$cite2=="YES",]))
DimPlot(MNP_LP, group.by = "cite")

#### Save Rdata objects  ####
saveRDS(MNP_LP, file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R1_AllCells_MNP_LP_wCITEdsbnorm.rds")
#saveRDS(MNP_LP,file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1Cprol_MNP_LP_wCITEdsbnorm.rds")

AllCells <- MNP_LP

##### ADD to other subsets of data ####
rm(normalized_matrix_pat4cLP,normalized_matrix_pat4SILP,normalized_matrix_pat5cLP,normalized_matrix_pat6cLP)
rm(CD14CD1C)
rm(MNP_LP)

#### R2 MNP LPs #####
MNP_LP <- readRDS("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R2_MNPs_MNP_LP.rds")
MNP_LP[["CITE"]] <- CreateAssayObject(data = pat04_CITE_wNA[,colnames(MNP_LP)])

MNP_LP@meta.data$cite <- AllCells@meta.data[colnames(MNP_LP),]$cite
MNP_LP@meta.data$cite2 <- AllCells@meta.data[colnames(MNP_LP),]$cite2
saveRDS(MNP_LP, file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R2_MNPs_MNP_LP_wCITEdsbnorm.rds")

#### R5 DC2 like cells ####
MNP_LP <- readRDS("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R5_DC2like.rds")
MNP_LP[["CITE"]] <- CreateAssayObject(data = pat04_CITE_wNA[,colnames(MNP_LP)])

MNP_LP@meta.data$cite <- AllCells@meta.data[colnames(MNP_LP),]$cite
MNP_LP@meta.data$cite2 <- AllCells@meta.data[colnames(MNP_LP),]$cite2
saveRDS(MNP_LP, file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R5_DC2like_wCITEdsbnorm.rds")

#### R3 CD14CD1C monomac_traj_noprol ####
MNP_LP <- readRDS("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1C_monomac_traj_noprol.rds")
MNP_LP[["CITE"]] <- CreateAssayObject(data = pat04_CITE_wNA[,colnames(MNP_LP)])

MNP_LP@meta.data$cite <- AllCells@meta.data[colnames(MNP_LP),]$cite
MNP_LP@meta.data$cite2 <- AllCells@meta.data[colnames(MNP_LP),]$cite2
saveRDS(MNP_LP, file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1C_monomac_traj_noprol_wCITEdsbnorm.rds")







