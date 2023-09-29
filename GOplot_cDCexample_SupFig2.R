#########################################################################################################
########################## F2 & SF2 - GO ontology in mature cDC #########################################
#########################################################################################################
rm(list=ls())
setwd("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Figures/F2_DC_cLPSILP")

#### Libraries ####
library(gplots)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(viridis)
library(ggpubr)
library(ggradar)
library(tidyr)
library(tidyverse)

#### Variables ####
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
mycols<- myColorRamp(mycols, seq(1:50))
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

matcDC_col <- c("cDC1"="#F564E3","cDC2"="#9590FF","cDC3"="#D89000")

#### mature cDC ####
GO_cDC1 <- read.csv("22_11_cDC1_BiologicalProcess.csv", header = T)
GO_cDC1$pop <- "cDC1"
GO_cDC1$Term <- as.character(GO_cDC1$Term)
GO_cDC2 <- read.csv("22_11_cDC2_BiologicalProcess.csv", header = T)
GO_cDC2$pop <- "cDC2"
GO_cDC2$Term <- as.character(GO_cDC2$Term)
GO_cDC3 <- read.csv("22_11_cDC3_BiologicalProcess.csv", header = T)
GO_cDC3$pop <- "cDC3"
GO_cDC3$Term <- as.character(GO_cDC3$Term)


GO_cDC1_sig <- GO_cDC1[GO_cDC1$Adjusted.P.value<0.05,]$Term #66
GO_cDC2_sig <- GO_cDC2[GO_cDC2$Adjusted.P.value<0.1,]$Term #27
GO_cDC3_sig <- GO_cDC3[GO_cDC3$Adjusted.P.value<0.05,]$Term #27

terms <- c(GO_cDC1_sig[1:10],GO_cDC2_sig[1:10],GO_cDC3_sig[1:10])

## DF for rdarplot
GO_df <- GO_cDC1[GO_cDC1$Term %in% terms,c("pop","Adjusted.P.value","Term")]
missing_t <- terms[!terms %in% GO_cDC1$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cDC1",1,item))}

GO_df <- rbind(GO_df, GO_cDC2[GO_cDC2$Term %in% terms,c("pop","Adjusted.P.value","Term")])
missing_t <- terms[!terms %in% GO_cDC2$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cDC2",1,item))}

GO_df <- rbind(GO_df, GO_cDC3[GO_cDC3$Term %in% terms,c("pop","Adjusted.P.value","Term")])
missing_t <- terms[!terms %in% GO_cDC3$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cDC3",1,item))}

## Adjusting dataframe
head(GO_df)
GO_df$Adjusted.P.value <- as.numeric(GO_df$Adjusted.P.value)
rownames(GO_df) <-NULL

#each term needs to be a column
GO_df <- spread(GO_df, key = Term, value = Adjusted.P.value)

head(GO_df)
colmax <- ncol(GO_df) #2:ncol value in following 
GO_df[,2:colmax] <- sqrt(-log10(GO_df[,2:colmax]))
min(GO_df[,2:colmax])
max(GO_df[,2:colmax])

GO_df_ord <- GO_df[,c("pop",
                      #cDC1
                      "positive regulation of cytokine production (GO:0001819)", 
                      "response to unfolded protein (GO:0006986)",
                      "protein targeting to ER (GO:0045047)",                                                       
                      "SRP-dependent cotranslational protein targeting to membrane (GO:0006614)",                   
                      "cotranslational protein targeting to membrane (GO:0006613)",                                 
                      "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay (GO:0000184)",           
                      "cytoplasmic translation (GO:0002181)",                                                       
                      "nuclear-transcribed mRNA catabolic process (GO:0000956)",                                    
                      "peptide biosynthetic process (GO:0043043)",                                                  
                      "translation (GO:0006412)",
                      
                      #cDC2                                                  
                      "positive regulation of epithelial cell proliferation (GO:0050679)",                          
                      "antigen processing and presentation of lipid antigen via MHC class Ib (GO:0048003)",         
                      "antigen processing and presentation, endogenous lipid antigen via MHC class Ib (GO:0048006)",
                      "antigen processing and presentation, exogenous lipid antigen via MHC class Ib (GO:0048007)", 
                      "positive regulation of T-helper 2 cell differentiation (GO:0045630)",
                      "positive regulation of interleukin-2 production (GO:0032743)", 
                      "positive regulation of interleukin-4 production (GO:0032753)",                               
                      "regulation of interleukin-4 production (GO:0032673)",
                      "positive regulation of T-helper cell differentiation (GO:0045624)",
                      
                      #cDC3
                      "inflammatory response (GO:0006954)",                                                         
                      "cellular response to molecule of bacterial origin (GO:0071219)",                             
                      "cellular response to lipopolysaccharide (GO:0071222)",                                       
                      "response to lipopolysaccharide (GO:0032496)",                                                
                      "granulocyte chemotaxis (GO:0071621)",
                      "neutrophil degranulation (GO:0043312)",                                                      
                      "neutrophil activation involved in immune response (GO:0002283)",                             
                      "neutrophil mediated immunity (GO:0002446)",                                                  
                      "cellular response to cytokine stimulus (GO:0071345)",
                      "cytokine-mediated signaling pathway (GO:0019221)"
)]

pdf(paste(dato,"SF2_matcDC_GOBioprocesses_v1.pdf"), height = 4, width = 5)
ggradar(GO_df_ord,
        # Plot area - change according to highest and lowest values
        values.radar = c("1", "0.05", "lowest"),
        grid.min = 0, grid.mid = 1.140627, grid.max = 4.9, #0.1padj=1, 0.05padj.=1.140627
        # Polygons
        group.line.width = 1, 
        group.point.size = 0,
        group.colours = matcDC_col, 
        # Background and grid lines
        background.circle.colour = "white",
        gridline.mid.colour = "black", #if set to logFC threshold set  todifferent colour
        # Text and labels
        #base.size = 5,
        legend.position = "bottom",
        label.centre.y = FALSE,
        label.gridline.min = FALSE,
        label.gridline.mid = FALSE,
        label.gridline.max = FALSE,
        axis.label.size = 4,
        axis.label.offset =1.05)
dev.off()

#### mature cDC - selected ####
# sel_terms <- c("response to unfolded protein (GO:0006986)",
#   "protein targeting to ER (GO:0045047)",  
#   "translation (GO:0006412)", #[This and previous two are interesting because it suggests the DC1 are secreting lots of something?] 
#   "positive regulation of epithelial cell proliferation (GO:0050679)", 
#   "antigen processing and presentation of lipid antigen via MHC class Ib (GO:0048003)",
#   "positive regulation of interleukin-2 production (GO:0032743)", 
#   "positive regulation of interleukin-4 production (GO:0032753)",
#   "positive regulation of T-helper 2 cell differentiation (GO:0045630)",
#   "positive regulation of T-helper cell differentiation (GO:0045624)",
#   "inflammatory response (GO:0006954)",
#   "cellular response to molecule of bacterial origin (GO:0071219)",
#   "granulocyte chemotaxis (GO:0071621)",
#   "neutrophil activation involved in immune response (GO:0002283)", 
#   "cellular response to cytokine stimulus (GO:0071345)",
#   "positive regulation of cytokine production (GO:0001819)")
sel_terms <- c("interferon-gamma-mediated signaling pathway (GO:0060333)",
               "positive regulation of cytokine production (GO:0001819)",
               "cellular response to cytokine stimulus (GO:0071345)",
               "inflammatory response (GO:0006954)",
               "cellular response to molecule of bacterial origin (GO:0071219)",
               "granulocyte chemotaxis (GO:0071621)",
               "regulation of interleukin-1 beta production (GO:0032651)",
               "regulation of complement activation (GO:0030449)",
               "positive regulation of T-helper cell differentiation (GO:0045624)",
               
               "positive regulation of epithelial cell proliferation (GO:0050679)", 
               
               # new
               "Fc-epsilon receptor signaling pathway (GO:0038095)",
               "antigen processing and presentation of peptide antigen via MHC class I (GO:0002474)",
               
               "response to unfolded protein (GO:0006986)",
               "protein targeting to ER (GO:0045047)",  
               "translation (GO:0006412)" #[This and previous two are interesting because it suggests the DC1 are secreting lots of something?]
)


terms <- sel_terms

GO_df <- GO_cDC1[GO_cDC1$Term %in% terms,c("pop","Adjusted.P.value","Term")]
missing_t <- terms[!terms %in% GO_cDC1$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cDC1",1,item))}

GO_df <- rbind(GO_df, GO_cDC2[GO_cDC2$Term %in% terms,c("pop","Adjusted.P.value","Term")])
missing_t <- terms[!terms %in% GO_cDC2$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cDC2",1,item))}

GO_df <- rbind(GO_df, GO_cDC3[GO_cDC3$Term %in% terms,c("pop","Adjusted.P.value","Term")])
missing_t <- terms[!terms %in% GO_cDC3$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cDC3",1,item))}

## Adjusting dataframe
head(GO_df)
GO_df$Adjusted.P.value <- as.numeric(GO_df$Adjusted.P.value)
rownames(GO_df) <-NULL

#each term needs to be a column
GO_df <- spread(GO_df, key = Term, value = Adjusted.P.value)

head(GO_df)
colmax <- ncol(GO_df) #2:ncol value in following 
GO_df[,2:colmax] <- sqrt(-log10(GO_df[,2:colmax]))
min(GO_df[,2:colmax])
max(GO_df[,2:colmax])

GO_df_ord <- GO_df[,c("pop",sel_terms)]

pdf(paste(dato,"SF2_matcDC_GOBioprocessesSelTerms_v3.pdf"), height = 4, width = 5)
ggradar(GO_df_ord,
        # Plot area - change according to highest and lowest values
        values.radar = c("1", "0.05", "lowest"),
        grid.min = 0, grid.mid = 1.140627, grid.max = 4.9, #0.1padj=1, 0.05padj.=1.140627
        # Polygons
        group.line.width = 1, 
        group.point.size = 0,
        group.colours = matcDC_col, 
        # Background and grid lines
        background.circle.colour = "white",
        gridline.mid.colour = "black", #if set to logFC threshold set  todifferent colour
        # Text and labels
        #base.size = 5,
        legend.position = "bottom",
        label.centre.y = FALSE,
        label.gridline.min = FALSE,
        label.gridline.mid = FALSE,
        label.gridline.max = FALSE,
        axis.label.size = 4,
        axis.label.offset =1.05)
dev.off()

#### SILP vs cLP - selected terms ####
sel_terms <- c( #cLP all
  "regulation of cellular response to heat (GO:1900034)",
  "regulation of cellular response to stress (GO:0080135)",
  "response to unfolded protein (GO:0006986)",
  "interferon-gamma-mediated signaling pathway (GO:0060333)",
  "positive regulation of tumor necrosis factor-mediated signaling pathway (GO:1903265)",
  "positive regulation of cell-substrate adhesion (GO:0010811)",
  # cDC1 cLP
  "MyD88-dependent toll-like receptor signaling pathway (GO:0002755)",
  "cellular response to cytokine stimulus (GO:0071345)",
  #SILP
  "intestinal cholesterol absorption (GO:0030299)",
  "regulation of complement activation (GO:0030449)")

terms_ord <- c( #cLP all
  "regulation of cellular response to heat (GO:1900034)",
  "regulation of cellular response to stress (GO:0080135)",
  "response to unfolded protein (GO:0006986)",
  "interferon-gamma-mediated signaling pathway (GO:0060333)",
  "positive regulation of tumor necrosis factor-mediated signaling pathway (GO:1903265)",
  "positive regulation of cell-substrate adhesion (GO:0010811)",
  # cDC1 cLP
  "MyD88-dependent toll-like receptor signaling pathway (GO:0002755)",
  #SILP
  "intestinal cholesterol absorption (GO:0030299)",
  "regulation of complement activation (GO:0030449)",
  "cellular response to cytokine stimulus (GO:0071345)")

# Tom's suggestions
sel_terms <- c("protein targeting to ER (GO:0045047)", 
               "response to unfolded protein (GO:0006986)",   #  [This and previous are interesting because it suggests the SILP DC are secreting more of something?] 
               "humoral immune response mediated by circulating immunoglobulin (GO:0002455)",      
               "positive regulation of cell population proliferation (GO:0008284)",   
               "regulation of complement activation (GO:0030449)",   
               "intestinal cholesterol absorption (GO:0030299)",   
               "cellular response to cytokine stimulus (GO:0071345)",
               "interferon-gamma-mediated signaling pathway (GO:0060333)",  
               "neutrophil activation involved in immune response (GO:0002283)", 
               "negative regulation of apoptotic process (GO:0043066)",
               "regulation of cellular response to heat (GO:1900034)" )

terms_ord <- c("cellular response to cytokine stimulus (GO:0071345)",
               "humoral immune response mediated by circulating immunoglobulin (GO:0002455)",  
               "protein targeting to ER (GO:0045047)",
               "intestinal cholesterol absorption (GO:0030299)",
               "regulation of complement activation (GO:0030449)",
               "positive regulation of cell population proliferation (GO:0008284)",  
               
               "interferon-gamma-mediated signaling pathway (GO:0060333)",
               "response to unfolded protein (GO:0006986)",   #  [This and previous are interesting because it suggests the SILP DC are secreting more of something?] 
               
               "negative regulation of apoptotic process (GO:0043066)",
               "regulation of cellular response to heat (GO:1900034)",
               "neutrophil activation involved in immune response (GO:0002283)")

# updated terms based on Tom and Bill request
sel_terms <- c("protein targeting to ER (GO:0045047)", 
               "response to unfolded protein (GO:0006986)",   #  [This and previous are interesting because it suggests the SILP DC are secreting more of something?] 
               "positive regulation of cell population proliferation (GO:0008284)",   
               "regulation of complement activation (GO:0030449)",   
               "intestinal cholesterol absorption (GO:0030299)",   
               "interferon-gamma-mediated signaling pathway (GO:0060333)", 
               "negative regulation of apoptotic process (GO:0043066)",
               "negative regulation of cytokine production (GO:0001818)",
               "regulation of interleukin-12 production (GO:0032655)",
               "MyD88-dependent toll-like receptor signaling pathway (GO:0002755)",
               "cytokine-mediated signaling pathway (GO:0019221)",
               "regulation of cellular response to stress (GO:0080135)",
               "positive regulation of NF-kappaB transcription factor activity (GO:0051092)",
               "positive regulation of cytokine production (GO:0001819)",
               "interleukin-1-mediated signaling pathway (GO:0070498)")

terms_ord <- c("interferon-gamma-mediated signaling pathway (GO:0060333)",
               "MyD88-dependent toll-like receptor signaling pathway (GO:0002755)",
               "positive regulation of cell population proliferation (GO:0008284)",  #2
               "protein targeting to ER (GO:0045047)", #3
               "regulation of complement activation (GO:0030449)", #4
               "intestinal cholesterol absorption (GO:0030299)",  #5
               "negative regulation of cytokine production (GO:0001818)", #6
               
               
               "regulation of interleukin-12 production (GO:0032655)",
               "interleukin-1-mediated signaling pathway (GO:0070498)",
               "negative regulation of apoptotic process (GO:0043066)",
               "regulation of cellular response to stress (GO:0080135)",
               "response to unfolded protein (GO:0006986)",
               "cytokine-mediated signaling pathway (GO:0019221)",
               "positive regulation of NF-kappaB transcription factor activity (GO:0051092)",
               "positive regulation of cytokine production (GO:0001819)")



#### cDC1 ####
GO_cDC1SILP <- read.csv("22_11_SILPcDC1_BiologicalProcess.csv", header = T)
GO_cDC1SILP$tissue <- "SILP"
GO_cDC1cLP <- read.csv("22_11_cLPcDC1_BiologicalProcess.csv", header = T)
GO_cDC1cLP$tissue <- "cLP"

GO_cDC1SILP_sig <- GO_cDC1SILP[GO_cDC1SILP$Adjusted.P.value<0.05,]$Term
GO_cDC1cLP_sig <- GO_cDC1cLP[GO_cDC1cLP$Adjusted.P.value<0.05,]$Term

#terms <- c(GO_cDC1SILP_sig[1:10],GO_cDC1cLP_sig [1:10])
terms <- sel_terms

GO_df <- GO_cDC1SILP[GO_cDC1SILP$Term %in% terms,c("tissue","Adjusted.P.value","Term")]
missing_t <- terms[!terms %in% GO_cDC1SILP$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("SILP",1,item))}

GO_df <- rbind(GO_df, GO_cDC1cLP[GO_cDC1cLP$Term %in% terms,c("tissue","Adjusted.P.value","Term")])
missing_t <- terms[!terms %in% GO_cDC1cLP$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cLP",1,item))}

## Adjusting dataframe
head(GO_df)
GO_df$Adjusted.P.value <- as.numeric(GO_df$Adjusted.P.value)
rownames(GO_df) <-NULL

#each term needs to be a column
GO_df <- spread(GO_df, key = Term, value = Adjusted.P.value)

head(GO_df)
colmax <- ncol(GO_df) #2:ncol value in following 
GO_df[,2:colmax] <- sqrt(-log10(GO_df[,2:colmax]))
min(GO_df[,2:colmax])
max(GO_df[,2:colmax])

GO_df_ord_2 <- GO_df[,c("tissue",terms_ord)]

pdf(paste(dato,"SIvsLI_cDC1_Bioprocesses_selv3.pdf"), height = 4, width = 5)
ggradar(GO_df_ord_2,
        # Plot area - change according to highest and lowest values
        values.radar = c("1", "0.05", "lowest"),
        grid.min = 0, grid.mid =  1.140627, grid.max = 4.8,
        # Polygons
        group.line.width = 1, 
        group.point.size = 0,
        group.colours = rev(hue_pal()(2)),
        # Background and grid lines
        background.circle.colour = "white",
        gridline.mid.colour = "black", #if set to logFC threshold set  todifferent colour
        # Text and labels
        #base.size = 5,
        legend.position = "bottom",
        label.centre.y = FALSE,
        label.gridline.min = FALSE,
        label.gridline.mid = FALSE,
        label.gridline.max = FALSE,
        axis.label.size = 4,
        axis.label.offset =1.05
)#+ theme(plot.margin=unit(c(0,2,0,2),"cm"))
dev.off()

#### cDC2 ####
GO_cDC2SILP <- read.csv("22_11_SILPcDC2_BiologicalProcess.csv", header = T)
GO_cDC2SILP$tissue <- "SILP"
GO_cDC2cLP <- read.csv("22_11_cLPcDC2_BiologicalProcess.csv", header = T)
GO_cDC2cLP$tissue <- "cLP"

GO_cDC2SILP_sig <- GO_cDC2SILP[GO_cDC2SILP$Adjusted.P.value<0.05,]$Term
GO_cDC2cLP_sig <- GO_cDC2cLP[GO_cDC2cLP$Adjusted.P.value<0.05,]$Term

#terms <- c(GO_cDC2SILP_sig[1:10],GO_cDC2cLP_sig [1:10])
terms <- sel_terms

GO_df <- GO_cDC2SILP[GO_cDC2SILP$Term %in% terms,c("tissue","Adjusted.P.value","Term")]
missing_t <- terms[!terms %in% GO_cDC2SILP$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("SILP",1,item))}

GO_df <- rbind(GO_df, GO_cDC2cLP[GO_cDC2cLP$Term %in% terms,c("tissue","Adjusted.P.value","Term")])
missing_t <- terms[!terms %in% GO_cDC2cLP$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cLP",1,item))}

## Adjusting dataframe
head(GO_df)
GO_df$Adjusted.P.value <- as.numeric(GO_df$Adjusted.P.value)
rownames(GO_df) <-NULL

#each term needs to be a column
GO_df <- spread(GO_df, key = Term, value = Adjusted.P.value)

head(GO_df)
colmax <- ncol(GO_df) #2:ncol value in following 
GO_df[,2:colmax] <- sqrt(-log10(GO_df[,2:colmax]))
min(GO_df[,2:colmax])
max(GO_df[,2:colmax])

GO_df_ord_2 <- GO_df[,c("tissue",terms_ord)]

pdf(paste(dato,"SIvsLI_cDC2_Bioprocesses_selv3.pdf"), height = 4, width = 5)
ggradar(GO_df_ord_2,
        # Plot area - change according to highest and lowest values
        values.radar = c("1", "0.05", "lowest"),
        grid.min = 0, grid.mid =  1.140627, grid.max = 3.8,
        # Polygons
        group.line.width = 1, 
        group.point.size = 0,
        group.colours = rev(hue_pal()(2)),
        # Background and grid lines
        background.circle.colour = "white",
        gridline.mid.colour = "black", #if set to logFC threshold set  todifferent colour
        # Text and labels
        #base.size = 5,
        legend.position = "bottom",
        label.centre.y = FALSE,
        label.gridline.min = FALSE,
        label.gridline.mid = FALSE,
        label.gridline.max = FALSE,
        axis.label.size = 4,
        axis.label.offset =1.05
)#+ theme(plot.margin=unit(c(0,2,0,2),"cm"))
dev.off()

#### cDC3 ####
GO_cDC3SILP <- read.csv("22_11_SILPcDC3_BiologicalProcess.csv", header = T)
GO_cDC3SILP$tissue <- "SILP"
GO_cDC3cLP <- read.csv("22_11_cLPcDC3_BiologicalProcess.csv", header = T)
GO_cDC3cLP$tissue <- "cLP"

GO_cDC3SILP_sig <- GO_cDC3SILP[GO_cDC3SILP$Adjusted.P.value<0.05,]$Term
GO_cDC3cLP_sig <- GO_cDC3cLP[GO_cDC3cLP$Adjusted.P.value<0.05,]$Term

#terms <- c(GO_cDC3SILP_sig[1:10],GO_cDC3cLP_sig [1:10])
terms <- sel_terms

GO_df <- GO_cDC3SILP[GO_cDC3SILP$Term %in% terms,c("tissue","Adjusted.P.value","Term")]
missing_t <- terms[!terms %in% GO_cDC3SILP$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("SILP",1,item))}

GO_df <- rbind(GO_df, GO_cDC3cLP[GO_cDC3cLP$Term %in% terms,c("tissue","Adjusted.P.value","Term")])
missing_t <- terms[!terms %in% GO_cDC3cLP$Term]
for (item in missing_t){
  GO_df <- rbind(GO_df, c("cLP",1,item))}

## Adjusting dataframe
head(GO_df)
GO_df$Adjusted.P.value <- as.numeric(GO_df$Adjusted.P.value)
rownames(GO_df) <-NULL

#each term needs to be a column
GO_df <- spread(GO_df, key = Term, value = Adjusted.P.value)

head(GO_df)
colmax <- ncol(GO_df) #2:ncol value in following 
GO_df[,2:colmax] <- sqrt(-log10(GO_df[,2:colmax]))
min(GO_df[,2:colmax])
max(GO_df[,2:colmax])

GO_df_ord_2 <- GO_df[,c("tissue",terms_ord)]

pdf(paste(dato,"SIvsLI_cDC3_Bioprocesses_selv3.pdf"), height = 4, width = 5)
ggradar(GO_df_ord_2,
        # Plot area - change according to highest and lowest values
        values.radar = c("1", "0.05", "lowest"),
        grid.min = 0, grid.mid =  1.140627, grid.max = 3.6,
        # Polygons
        group.line.width = 1, 
        group.point.size = 0,
        group.colours = rev(hue_pal()(2)),
        # Background and grid lines
        background.circle.colour = "white",
        gridline.mid.colour = "black", #if set to logFC threshold set  todifferent colour
        # Text and labels
        #base.size = 5,
        legend.position = "bottom",
        label.centre.y = FALSE,
        label.gridline.min = FALSE,
        label.gridline.mid = FALSE,
        label.gridline.max = FALSE,
        axis.label.size = 4,
        axis.label.offset =1.05
)#+ theme(plot.margin=unit(c(0,2,0,2),"cm"))
dev.off()
