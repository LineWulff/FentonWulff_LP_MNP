###################################################################################
###### ----------------- normalising CITE-seq with DSB ----------------- ##########
###################################################################################
setwd("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/citeseq")

rm(list=ls())

#### Read in libraries ####
library(gplots)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(magicBatch)
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
project <- "citeseq_MNP"


#### Read in raw data to identify empty droplets/debris ####
pat04_GEX <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/18/filtered_feature_bc_matrix/")
pat04_CITE <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/18/citeseq_data_full/umi_count/",gene.column=1)
rownames(pat04_CITE) <- c("Iso_IgG1","Iso_IgG2a","Iso_IgG2b","CD14","CD209","CD1c","CD55","CD5","CD207","CD206","CD11a","CD103","CD11c","unmapped")
colnames(pat04_CITE) <- paste0(colnames(pat04_CITE),"-1")
dim(pat04_GEX)
dim(pat04_CITE)

pat04_SILP_GEX <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/filtered_feature_bc_matrix/")
pat04_SILP_CITE <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/citeseq_data/umi_count/",gene.column=1)
rownames(pat04_SILP_CITE) <- c("Iso_IgG1","Iso_IgG2a","Iso_IgG2b","CD14","CD209","CD1c","CD55","CD5","CD207","CD206","CD11a","CD103","CD11c","unmapped")
colnames(pat04_SILP_CITE) <- paste0(colnames(pat04_SILP_CITE),"-1")
dim(pat04_SILP_GEX)
dim(pat04_SILP_CITE)

pat05_GEX <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_5/filtered_feature_bc_matrix/")
pat05_CITE <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_5/citeseq_data/umi_count/",gene.column=1)
rownames(pat05_CITE) <- c("Iso_IgG1","Iso_IgG2a","Iso_IgG2b","CD14","CD209","CD1c","CD55","CD5","CD207","CD206","CD11a","CD103","CD11c","unmapped")
colnames(pat05_CITE) <- paste0(colnames(pat05_CITE),"-1")
dim(pat05_GEX)
dim(pat05_CITE)

pat06_GEX <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_6/filtered_feature_bc_matrix/")
pat06_CITE <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_6/citeseq_data/umi_count/",gene.column=1)
rownames(pat06_CITE) <- c("Iso_IgG1","Iso_IgG2a","Iso_IgG2b","CD14","CD209","CD1c","CD55","CD5","CD207","CD206","CD11a","CD103","CD11c","unmapped")
colnames(pat06_CITE) <- paste0(colnames(pat06_CITE),"-1")
dim(pat05_GEX)
dim(pat05_CITE)

###### As Seurat objects ######
#### pat4_cLP ####
pat4_cLP <- CreateSeuratObject(counts = pat04_GEX, project = project, min.cells = 3, min.features = 100)
pat4_cLP[["percent.mt"]] <- PercentageFeatureSet(pat4_cLP, pattern = "^MT-")

GEX_cells <- rownames(pat4_cLP@meta.data)
CITE_cells <- colnames(pat04_CITE)
GEX_cells[!GEX_cells %in% CITE_cells] #how many cells are we excluding from the cite-seq data
CITE_cells[!CITE_cells %in% GEX_cells] #how many from the gene expression?

pat4_cLP <- subset(pat4_cLP, cells = GEX_cells[GEX_cells %in% CITE_cells])
pat4_cLP[["CITE"]] <- CreateAssayObject(counts = pat04_CITE)

#from QC cut offs
neg_object <- subset(pat4_cLP, subset = nFeature_RNA < 800)
pat4_cLP <- subset(pat4_cLP, subset = nFeature_RNA > 800 & nFeature_RNA < 5500 & percent.mt < 10)

# non sparse CITEseq data actually store better in a regular materix so the as.matrix() call is not memory intensive.
neg_adt_matrix = GetAssayData(neg_object, assay = "CITE", slot = 'counts') %>% as.matrix()
positive_adt_matrix = GetAssayData(pat4_cLP, assay = "CITE", slot = 'counts') %>% as.matrix()

isotypes <- c("Iso-IgG1","Iso-IgG2a","Iso-IgG2b")

# normalize the data with dsb
# make sure you've run devtools::install_github(repo = 'MattPM/dsb')

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,
                                        empty_drop_matrix = neg_adt_matrix,
                                        use.isotype.control = TRUE,
                                        isotype.control.name.vec = isotypes)

#colnames(normalized_matrix) <- str_replace(colnames(normalized_matrix),"-1","")
#save the normalized matrix to be able to upload to other objects
save(normalized_matrix,file= "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/18/analysis/Rdata/Pat4AllCells_CITEseq_dsbnorm.Rdata")

#rename cells in norm matrix according to object you want to fit to
load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/analysis/Rdata/R2_pat4_cLP.Rdata")
rownames(pat4_cLP@meta.data)
#colnames(normalized_matrix) <- paste0("cLP_",colnames(normalized_matrix))

#load in the data you want to plot with
normalized_matrix <- normalized_matrix[,rownames(pat4_cLP@meta.data)] #the data you want to pull 

# now add the normalized dat back to the object (the singlets defined above as "object")
#pat4_cLP = SetAssayData(object = pat4_cLP, assay = "CITE",slot = "data", new.data = normalized_matrix)

pat4_cLP[["CITE"]] <- CreateAssayObject(data = normalized_matrix)

save(pat4_cLP,file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/18/analysis/Rdata/R2_pat4_wCITE_dsbnorm.Rdata")

#Idents(pat4_cLP) <- pat4_cLP@meta.data$RNA_snn_res.0.4
#DoHeatmap(pat4_cLP, features = rownames(pat4_cLP@assays$CITE@data), assay = "CITE", slot = "data", 
#           raster = F, lines.width = 20, disp.max = 2)


#### pat4_SILP ####
pat4_SILP <- CreateSeuratObject(counts = pat04_SILP_GEX, project = project, min.cells = 3, min.features = 100)
pat4_SILP[["percent.mt"]] <- PercentageFeatureSet(pat4_SILP, pattern = "^MT-")

GEX_cells <- rownames(pat4_SILP@meta.data)
CITE_cells <- colnames(pat04_SILP_CITE)
dim(pat04_SILP_GEX)
#dim(pat04_SILP_CITE)
length(GEX_cells[GEX_cells %in% CITE_cells]) #how many cells are we excluding from the cite-seq data
length(CITE_cells[CITE_cells %in% GEX_cells]) #how many from the gene expression?


pat4_SILP <- subset(pat4_SILP, cells = GEX_cells[GEX_cells %in% CITE_cells])
pat4_SILP[["CITE"]] <- CreateAssayObject(counts = pat04_SILP_CITE[,GEX_cells[GEX_cells %in% CITE_cells]])

#from QC cut offs
neg_object <- subset(pat4_SILP, subset = nFeature_RNA < 700)
pat4_SILP <- subset(pat4_SILP, subset = nFeature_RNA > 700 & nFeature_RNA < 5900 & percent.mt < 10)

# non sparse CITEseq data actually store better in a regular materix so the as.matrix() call is not memory intensive.
neg_adt_matrix = GetAssayData(neg_object, assay = "CITE", slot = 'counts') %>% as.matrix()
positive_adt_matrix = GetAssayData(pat4_SILP, assay = "CITE", slot = 'counts') %>% as.matrix()

isotypes <- c("Iso-IgG1","Iso-IgG2a","Iso-IgG2b")

# normalize the data with dsb
# make sure you've run devtools::install_github(repo = 'MattPM/dsb')

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,
                                        empty_drop_matrix = neg_adt_matrix,
                                        use.isotype.control = TRUE,
                                        isotype.control.name.vec = isotypes)

colnames(normalized_matrix) <- str_replace(colnames(normalized_matrix),"-1","")
#save the normalized matrix to be able to upload to other objects
save(normalized_matrix,file= "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/analysis/Rdata/Pat4AllCells_CITEseq_dsbnorm.Rdata")

#rename cells in norm matrix according to object you want to fit to
load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/analysis/Rdata/R1_pat4_SILP.Rdata")
rownames(pat4_SILP@meta.data)
colnames(normalized_matrix) <- paste0("SILP_",colnames(normalized_matrix))

#load in the data you want to plot with
length(colnames(pat4_SILP))
length(colnames(pat4_SILP)[colnames(pat4_SILP) %in% colnames(normalized_matrix)])
cells <- colnames(pat4_SILP)[colnames(pat4_SILP) %in% colnames(normalized_matrix)]
normalized_matrix <- normalized_matrix[,cells] #the data you want to pull 

pat4_SILP <- subset(pat4_SILP, cells = cells)

# now add the normalized dat back to the object (the singlets defined above as "object")
#pat4_SILP = SetAssayData(object = pat4_SILP, assay = "CITE",slot = "data", new.data = normalized_matrix)

pat4_SILP[["CITE"]] <- CreateAssayObject(data = normalized_matrix)

save(pat4_SILP,file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/analysis/Rdata/R1_pat4_wCITE_dsbnorm.Rdata")

Idents(pat4_SILP) <- pat4_SILP@meta.data$RNA_snn_res.0.4
DoHeatmap(pat4_SILP, features = rownames(pat4_SILP@assays$CITE@data), assay = "CITE", slot = "data", 
raster = F, lines.width = 20, disp.max = 2)


#### pat5 cLP ####
pat5_cLP <- CreateSeuratObject(counts = pat05_GEX, project = project, min.cells = 3, min.features = 100)
pat5_cLP[["percent.mt"]] <- PercentageFeatureSet(pat5_cLP, pattern = "^MT-")

## pat5 is the one with no CD11c and CD209 and Isotypes!!
pat05_CITE <- pat05_CITE[c("CD14","CD1c","CD55","CD5","CD207","CD206","CD11a","CD103","unmapped"),]

GEX_cells <- rownames(pat5_cLP@meta.data)
CITE_cells <- colnames(pat05_CITE)
dim(pat05_GEX)
#dim(pat05_CITE)
length(GEX_cells[GEX_cells %in% CITE_cells]) #how many cells are we excluding from the cite-seq data
length(CITE_cells[CITE_cells %in% GEX_cells]) #how many from the gene expression?


pat5_cLP <- subset(pat5_cLP, cells = GEX_cells[GEX_cells %in% CITE_cells])
pat5_cLP[["CITE"]] <- CreateAssayObject(counts = pat05_CITE[,GEX_cells[GEX_cells %in% CITE_cells]])

#from QC cut offs
neg_object <- subset(pat5_cLP, subset = nFeature_RNA < 700)
pat5_cLP <- subset(pat5_cLP, subset = nFeature_RNA > 700 & nFeature_RNA < 6000 & percent.mt < 10)

# non sparse CITEseq data actually store better in a regular materix so the as.matrix() call is not memory intensive.
neg_adt_matrix = GetAssayData(neg_object, assay = "CITE", slot = 'counts') %>% as.matrix()
positive_adt_matrix = GetAssayData(pat5_cLP, assay = "CITE", slot = 'counts') %>% as.matrix()

# no isotype controls in this sample
#isotypes <- c("Iso-IgG1","Iso-IgG2a","Iso-IgG2b")

# normalize the data with dsb
# make sure you've run devtools::install_github(repo = 'MattPM/dsb')

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,
                                        empty_drop_matrix = neg_adt_matrix,
                                        use.isotype.control = FALSE)

colnames(normalized_matrix) <- str_replace(colnames(normalized_matrix),"-1","")
#save the normalized matrix to be able to upload to other objects
save(normalized_matrix,file= "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_5/analysis/Rdata/Pat5AllCells_CITEseq_dsbnorm.Rdata")

#rename cells in norm matrix according to object you want to fit to
load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_5/analysis/Rdata/R1_pat5_cLP.Rdata")
rownames(pat5_cLP@meta.data)
colnames(normalized_matrix) <- paste0("cLP_",colnames(normalized_matrix))

#load in the data you want to plot with
length(colnames(pat5_cLP))
length(colnames(pat5_cLP)[colnames(pat5_cLP) %in% colnames(normalized_matrix)])
cells <- colnames(pat5_cLP)[colnames(pat5_cLP) %in% colnames(normalized_matrix)]
normalized_matrix <- normalized_matrix[,cells] #the data you want to pull 

pat5_cLP <- subset(pat5_cLP, cells = cells)

# now add the normalized dat back to the object (the singlets defined above as "object")
#pat5_cLP = SetAssayData(object = pat5_cLP, assay = "CITE",slot = "data", new.data = normalized_matrix)

pat5_cLP[["CITE"]] <- CreateAssayObject(data = normalized_matrix)

save(pat5_cLP,file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_5/analysis/Rdata/R1_pat5_wCITE_dsbnorm.Rdata")

Idents(pat5_cLP) <- pat5_cLP@meta.data$RNA_snn_res.0.4
DoHeatmap(pat5_cLP, features = rownames(pat5_cLP@assays$CITE@data), assay = "CITE", slot = "data", 
          raster = F, lines.width = 20, disp.max = 2)

#### Pat6 cLP ####
pat6_cLP <- CreateSeuratObject(counts = pat06_GEX, project = project, min.cells = 3, min.features = 100)
pat6_cLP[["percent.mt"]] <- PercentageFeatureSet(pat6_cLP, pattern = "^MT-")

GEX_cells <- rownames(pat6_cLP@meta.data)
length(GEX_cells)
CITE_cells <- colnames(pat06_CITE)
length(GEX_cells[GEX_cells %in% CITE_cells]) #how many cells are we excluding from the cite-seq data
length(CITE_cells[CITE_cells %in% GEX_cells]) #how many from the gene expression?

pat6_cLP <- subset(pat6_cLP, cells = GEX_cells[GEX_cells %in% CITE_cells])
pat6_cLP[["CITE"]] <- CreateAssayObject(counts = pat06_CITE[,GEX_cells[GEX_cells %in% CITE_cells]])

#from QC cut offs
neg_object <- subset(pat6_cLP, subset = nFeature_RNA < 700)
pat6_cLP <- subset(pat6_cLP, subset = nFeature_RNA > 700 & nFeature_RNA < 5700 & percent.mt < 10)

# non sparse CITEseq data actually store better in a regular materix so the as.matrix() call is not memory intensive.
neg_adt_matrix = GetAssayData(neg_object, assay = "CITE", slot = 'counts') %>% as.matrix()
positive_adt_matrix = GetAssayData(pat6_cLP, assay = "CITE", slot = 'counts') %>% as.matrix()

isotypes <- c("Iso-IgG1","Iso-IgG2a","Iso-IgG2b")

# normalize the data with dsb
# make sure you've run devtools::install_github(repo = 'MattPM/dsb')

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,
                                        empty_drop_matrix = neg_adt_matrix,
                                        use.isotype.control = TRUE,
                                        isotype.control.name.vec = isotypes)

colnames(normalized_matrix) <- str_replace(colnames(normalized_matrix),"-1","")
#save the normalized matrix to be able to upload to other objects
save(normalized_matrix,file= "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_6/analysis/Rdata/Pat6AllCells_CITEseq_dsbnorm.Rdata")

#rename cells in norm matrix according to object you want to fit to
load("/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_6/analysis/Rdata/R1_pat6_cLP.Rdata")
rownames(pat6_cLP@meta.data)
colnames(normalized_matrix) <- paste0("cLP_",colnames(normalized_matrix))

#load in the data you want to plot with
cells = rownames(pat6_cLP@meta.data)[rownames(pat6_cLP@meta.data) %in% colnames(normalized_matrix)]
normalized_matrix <- normalized_matrix[,cells] #the data you want to pull 

# now add the normalized dat back to the object (the singlets defined above as "object")
#pat6_cLP = SetAssayData(object = pat6_cLP, assay = "CITE",slot = "data", new.data = normalized_matrix)
pat6_cLP <- subset(pat6_cLP, cells = cells)
pat6_cLP[["CITE"]] <- CreateAssayObject(data = normalized_matrix)

save(pat6_cLP,file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_6/analysis/Rdata/R1_pat6_wCITE_dsbnorm.Rdata")

Idents(pat6_cLP) <- pat6_cLP@meta.data$RNA_snn_res.0.4
DoHeatmap(pat6_cLP, features = rownames(pat6_cLP@assays$CITE@data), assay = "CITE", slot = "data", 
           raster = F, lines.width = 20, disp.max = 2)
