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
library(ggplot2)
library(viridis)
library(dsb)
library(Seurat)

#### Variables to use ####
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# project
project <- "citeseq_MNP"

#### Read in raw data to identify empty droplets/debris ####
## EXAMPLE OF DSB NORMALIZATION FOR ONE SAMPLE 
## ALL OF THESE ARE CONCATENATED AND THEN ADDED TO SEURAT OBJECTS
## (SEE LINE 96)
pat04_SILP_GEX <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/filtered_feature_bc_matrix/")
pat04_SILP_CITE <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/citeseq_data/umi_count/",gene.column=1)
rownames(pat04_SILP_CITE) <- c("Iso_IgG1","Iso_IgG2a","Iso_IgG2b","CD14","CD209","CD1c","CD55","CD5","CD207","CD206","CD11a","CD103","CD11c","unmapped")
colnames(pat04_SILP_CITE) <- paste0(colnames(pat04_SILP_CITE),"-1")
dim(pat04_SILP_GEX)
dim(pat04_SILP_CITE)

###### As Seurat objects ######

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

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,
                                        empty_drop_matrix = neg_adt_matrix,
                                        use.isotype.control = TRUE,
                                        isotype.control.name.vec = isotypes)

colnames(normalized_matrix) <- str_replace(colnames(normalized_matrix),"-1","")
#save the normalized matrix to be able to upload to other objects
saveRDS(normalized_matrix,file= "Pat4AllCells_CITEseq_dsbnorm.rds")

#rename cells in norm matrix according to object you want to fit to
load("R1_pat4_SILP.Rdata")
rownames(pat4_SILP@meta.data)
colnames(normalized_matrix) <- paste("SILP_",colnames(normalized_matrix), sep="")

#load in the data you want to plot with
length(colnames(pat4_SILP))
length(colnames(pat4_SILP)[colnames(pat4_SILP) %in% colnames(normalized_matrix)])
cells <- colnames(pat4_SILP)[colnames(pat4_SILP) %in% colnames(normalized_matrix)]
normalized_matrix <- normalized_matrix[,cells] #the data you want to pull 

pat4_SILP <- subset(pat4_SILP, cells = cells)

# now add the normalized dat back to the object (the singlets defined above as "object")
#pat4_SILP = SetAssayData(object = pat4_SILP, assay = "CITE",slot = "data", new.data = normalized_matrix)

pat4_SILP[["CITE"]] <- CreateAssayObject(data = normalized_matrix)

saveRDS(pat4_SILP,file="/R1_pat4_wCITE_dsbnorm.rds")

Idents(pat4_SILP) <- pat4_SILP@meta.data$RNA_snn_res.0.4
DoHeatmap(pat4_SILP, features = rownames(pat4_SILP@assays$CITE@data), assay = "CITE", slot = "data", raster = F, lines.width = 20, disp.max = 2)


#### Adding concatenated DSB norm. CITE matrix for all samples with CITE-seq to Seurat object ####
