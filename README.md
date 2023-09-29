# FentonWulff_LP_MNP
This repository contains all the most essential scripts, which were used in analysing the single cell data for the preprint **Fenton, Wulff et al. BioRxiv 2021**.

QC of single cells was first assessed per sample based on number of read counts, feature/gene counts and mitochondrial gene percentage. Barcodes determined to come from debris or doublets by these parameters were removed before the samples were integrated with Seurat's anchor integration and scaled while removing effects of cell cycle, MT gene load, read and gene depth. Based on 15 PC's UMAP and clustering were calculated and contaminating cell types were removed. The remaining cells were rescaled, pca rerun, reclustered and a new UMAP was run (still on 15 PCs).
The supercluster of CD14+ and CD1C+ cells was again subsetted, rescaled, reclustered etc. at resolution 2.8.

**ProgenyDorotheaAnalysis.R** - example of how analysis with Dorothea and Progeny were run.

**Fig4_putprecursors.R** - Analysis for Figure 4 and Supplementary Figure 4.
