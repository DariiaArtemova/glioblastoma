
## код не передбачено для повного запуску, натомість використовувався для інтерактивної роботи та аналізу в блокноті RMarkdown. Перепрошую за незручності.
#building libraries
install.packages("BiocManager")
install.packages("SeuratObject")
install.packages("SeuratObject")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
install.packages("enrichR")
install.packages("enrichplot")
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
remotes::install_github("satijalab/seurat-object", ref = "develop")
BiocManager::install("EnhancedVolcano")

library(tidyverse)
library(readxl)
library(Seurat)
library(SeuratObject)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichR)
library(enrichplot)

#establishing dataset from tumor

setwd("D:/timely/Artests/Project")
obj <- readRDS("01_initial_processing.Rds")

#intro
testobj <- obj
DimPlot(testobj, reduction = "umap", group.by = "Location")+ 
  ggtitle("Клітини зразка за локалізацією") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, hjust = 1.5), margin = margin(b = 10))

# demonstrating tumor cell clusters using basic markers

VlnPlot(testobj, features = c("SOX9", "EGFR"))
FeaturePlot(testobj, features = c("EGFR","SOX9", "GFAP", "SOX2"))

# just in case 

testobj <- SCTransform(testobj)
testobj <- RunPCA(object = testobj)
testobj <- FindNeighbors(object = testobj, dims = 1:30)
testobj <- FindClusters(object = testobj)
testobj <- RunUMAP(object = testobj, dims = 1:30)

# subsetting potential tumor cells

Idents(testobj) <- "seurat_clusters"
total_tumor <- subset(x = testobj, idents = c("2", "3", "6", "11", "13"), subset = PTPRC == 0)

FeaturePlot(total_tumor, features = c("EGFR", "PTPRC"))
DimPlot(total_tumor, reduction = "umap", group.by = "Selection")

total_tumor <- SCTransform(total_tumor)
total_tumor <- ScaleData(object = total_tumor)
total_tumor <- RunPCA(object = total_tumor)
total_tumor <- FindNeighbors(object = total_tumor, dims = 1:30)
total_tumor <- FindClusters(object = total_tumor)
total_tumor <- RunUMAP(object = total_tumor, dims = 1:30)

# identifying clusters
tumor_cells <- PrepSCTFindMarkers(testobj)
tumor_cells <- FindAllMarkers(tumor_cells, features = c("EGFR", "SOX9"))
tumor_cells <- subset(tumor_cells, subset = avg_log2FC > 0)

# clusters 2, 3, 6, 11, 13 identified

DimPlot(total_tumor, reduction = "umap", group.by = "Selection")
DimPlot(total_tumor, reduction = "umap", group.by = "Location")
DimPlot(total_tumor, reduction = "umap", group.by = "Patient")

FeaturePlot(total_tumor, features = c("GFAP", "SOX2", "SOX9", "EGFR"))
VlnPlot(testobj, features = c("GFAP", "SOX2", "EGFR"))

FeaturePlot(testobj, features = c("CD133", "CD44", "CD15", "CD70", "S100A4", "ALDH1A3", "NANOG", "OCT4", "NESTIN"))
FeaturePlot(total_tumor, features = c("CD133", "CD44", "CD15", "CD70", "S100A4", "ALDH1A3", "NANOG", "OCT4", "NESTIN"))

#compare total tumor cell expression to all cells given
Idents(testobj) <- "combined_clusters"
testobj$combined_clusters <- as.character(testobj$seurat_clusters)
testobj$combined_clusters[testobj$seurat_clusters %in% c("2", "3", "6", "11", "13")] <- "GBM_cluster"
testobj$combined_clusters[!testobj$seurat_clusters %in% c("2", "3", "6", "11", "13")] <- "Others_cluster"
Idents(testobj) <- "combined_clusters"

tot_to_all <- PrepSCTFindMarkers(testobj) # for upregulated all tumor to total
Idents(testobj) <- "combined_clusters"
tot_to_all <- FindMarkers(tot_to_all, 
                          ident.1 = "GBM_cluster",
                          ident.2 = NULL,
                          assay = "SCT",
                          logfc.threshold = 0, 
                          min.pct = 0.1, 
                          test.use = "wilcox", 
                          only.pos = FALSE )


tot_to_all <- PrepSCTFindMarkers(testobj) # for downregulated all tumor to total
tot_to_all <- FindMarkers(tot_to_all, 
                          ident.1 = "GBM_cluster",
                          ident.2 = NULL,
                          assay = "SCT",
                          logfc.threshold = 0, 
                          min.pct = 0.1, 
                          test.use = "wilcox", 
                          only.pos = FALSE )

Idents(testobj) <- "combined_clusters"
EnhancedVolcano(tot_to_all,
                lab = rownames(tot_to_all),
               title = "Диференційна експресія генів неопластичного кластера", 
               selectLab = c("EGFR","SOX9", "SOX2", "GFAP", #basics
                             "CHL1", "POUSF2", "OLIG2", "SALL2", "NFIB", #verifiably upregulated
                             "MBP", "OPALIN", #oligodendrocites
                             "GPR17", #OPCs
                             "L1CAM", #neurons
                             "ALDH1L1", "WIF1","NTSR2", #astrocites
                             "CD133", "CD44", "CD15", "CD70", "S100A4", "ALDH1A3", "NANOG", "OCT4", "NESTIN", #other glioblastoma-associated markers
                             "CNN3", "PTPRC",
                             "PDGFB", "CLU", "TNC"),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                x = 'avg_log2FC',
                y = 'p_val')

Idents(testobj) <- "combined_clusters"
EnhancedVolcano(tot_to_all_down,
                lab = rownames(tot_to_all_down),
                x = 'avg_log2FC',
                y = 'p_val')


FeaturePlot(total_tumor, features = c("GFAP", "SOX2", "SOX9", "EGFR"))


#compare core tumor cell expression to all cells given

testobj$tumor_location_clusters <- ifelse(testobj$Location == "Tumor" & 
                                        testobj$combined_clusters == "GBM_cluster", 
                                      "core_tumor", 
                                      "other")
Idents(testobj) <- "tumor_location_clusters"
DimPlot(testobj, reduction = "umap", group.by = "tumor_location_clusters")
DimPlot(testobj, reduction = "umap", group.by = "combined_clusters")

tot_to_core <- PrepSCTFindMarkers(testobj)
tot_to_core <- FindMarkers(tot_to_core, 
                          ident.1 = "core_tumor",
                          ident.2 = NULL,
                          assay = "SCT",
                          logfc.threshold = 0, 
                          min.pct = 0.1, 
                          test.use = "wilcox", 
                          only.pos = FALSE )

Idents(testobj) <- "tumor_location_clusters"

EnhancedVolcano(tot_to_core,
                lab = rownames(tot_to_core),
                title = "Диференційна експресія неопластичних клітин кору",
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                #max.overlaps = 10,
                x = 'avg_log2FC',
                y = 'p_val')

EnhancedVolcano(tot_to_core,
                lab = rownames(tot_to_core),
                title = "Диференційна експресія неопластичних клітин кору",
                selectLab = c("EGFR","SOX9", "SOX2", "GFAP", #basics
                              "CHL1", "POUSF2", "OLIG2", "SALL2", "NFIB", #verifiably upregulated
                              "MBP", "OPALIN", #oligodendrocites
                              "GPR17", #OPCs
                              "L1CAM", #neurons
                              "ALDH1L1", "WIF1","NTSR2", #astrocites
                              "CD133", "CD44", "CD15", "CD70", "S100A4", "ALDH1A3", "NANOG", "OCT4", "NESTIN", #other glioblastoma-associated markers
                              "CNN3", "PTPRC"),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                max.overlaps = 20,
                x = 'avg_log2FC',
                y = 'p_val')


#compare periphery tumor cell expression to all cells given

testobj$pery.tumor_location_clusters <- ifelse(testobj$Location == "Periphery" & 
                                            testobj$combined_clusters == "GBM_cluster", 
                                          "pery_tumor", 
                                          "other")

Idents(testobj) <- "pery.tumor_location_clusters"

DimPlot(testobj, reduction = "umap", group.by = "pery.tumor_location_clusters")
DimPlot(testobj, reduction = "umap", group.by = "combined_clusters")

tot_to_pery <- PrepSCTFindMarkers(testobj)
tot_to_pery <- FindMarkers(tot_to_pery, 
                          ident.1 = "pery_tumor",
                          ident.2 = NULL,
                          assay = "SCT",
                          logfc.threshold = 0, 
                          min.pct = 0.1, 
                          test.use = "wilcox", 
                          only.pos = FALSE )

Idents(testobj) <- "tumor_location_clusters"
Idents(testobj) <- "combined_clusters"


EnhancedVolcano(tot_to_pery,
                lab = rownames(tot_to_pery),
                title = "Диференційна експресія неопластичних клітин периферії",
                selectLab = c("EGFR","SOX9", "SOX2", "GFAP", #basics
                              "CHL1", "POUSF2", "OLIG2", "SALL2", "NFIB", #verifiably upregulated
                              "MBP", "OPALIN", #oligodendrocites
                              "GPR17", #OPCs
                              "L1CAM", #neurons
                              "ALDH1L1", "WIF1","NTSR2", #astrocites
                              "CD133", "CD44", "CD15", "CD70", "S100A4", "ALDH1A3", "NANOG", "OCT4", "NESTIN", #other glioblastoma-associated markers
                              "CNN3", "PTPRC"),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                x = 'avg_log2FC',
                y = 'p_val')

EnhancedVolcano(tot_to_pery,
                lab = rownames(tot_to_pery),
                title = "Диференційна експресія неопластичних клітин периферії",
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                #max.overlaps = 10,
                x = 'avg_log2FC',
                y = 'p_val')


## deining the most overexpressed 

top150total <- subset(tot_to_all, subset = avg_log2FC > 2.5 & p_val < 0.05) #subsetting total tumor
the150total <- top150total[order(top150total$p_val_adj), ][1:353, ]
the150total <- head(the150total, 150)
top150total_genes <- rownames(the150total)
#top150total_genes <- paste(top150total_genes, collapse = " ")


the150downtotal <- subset(tot_to_all, subset = avg_log2FC < -2.5 & p_val < 0.05)
the150downtotal <- the150downtotal[order(the150downtotal$p_val_adj), ][1:637, ]
the150downtotal <- head(the150downtotal, 150)
top150downtotal_genes <- rownames(the150downtotal)
#top150downtotal_genes <- paste(top150downtotal_genes, collapse = " ")


top150core <- subset (tot_to_core, subset = avg_log2FC > 2.5 & p_val < 0.05) #subsetting core tumor 
the150core <- top150core[order(top150core$p_val_adj), ][1:294, ]
top150core <- head(the150core, 150)
top150core_genes <- rownames(top150core)
#top150core_genes <- paste(top150core_genes, collapse = " ")

tot_to_core_down <- subset(tot_to_core, subset = avg_log2FC < -2.5 & p_val < 0.05)
the150downcore <- tot_to_core_down[order(tot_to_core_down$p_val_adj), ][1:635, ]
the150downcore <- head(the150downcore, 150)
top150downcore_genes <- rownames(the150downcore)
#top150downcore_genes <- paste(top150downcore_genes, collapse = " ")


top150pery <- subset (tot_to_pery, subset = avg_log2FC > 1.5 & p_val < 0.05) #subsetting periphery tumor
the150pery <- top150pery[order(top150pery$p_val_adj), ][1:605, ]
top150pery <- head(the150pery, 150)
top150pery_genes <- rownames(top150pery)
#top150pery_genes <- paste(top150pery_genes, collapse = " ")

top150downpery <- subset (tot_to_pery, subset = avg_log2FC > -1.5 & p_val < 0.05) #subsetting periphery tumor
the150downpery <- top150downpery[order(top150downpery$p_val_adj), ][1:93, ]
top150downpery <- head(the150downpery, 150)
top150downpery_genes <- rownames(top150downpery)
#top150downpery_genes <- paste(top150downpery_genes, collapse = " ")

##functional analysis 

#entrez ids for them
topcore_entrez_ids <- mapIds(org.Hs.eg.db, keys = top150core_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
topcoredown_entrez_ids <- mapIds(org.Hs.eg.db, keys = top150downcore_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

toppery_entrez_ids <- mapIds(org.Hs.eg.db, keys = top150pery_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
topperydown_entrez_ids <- mapIds(org.Hs.eg.db, keys = top150downpery_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

toptotal_entrez_ids <- mapIds(org.Hs.eg.db, keys = top150total_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
toptotaldown_entrez_ids <- mapIds(org.Hs.eg.db, keys = top150downtotal_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# GO enrichment 
corego_result <- enrichGO(gene = topcore_entrez_ids, #core to all upreg
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
coredown_go_result <- enrichGO(gene = topcoredown_entrez_ids,  #core to all downreg
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)


perygo_result <- enrichGO(gene = toppery_entrez_ids, #periphery to all upreg
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
perydown_go_result <- enrichGO(gene = topperydown_entrez_ids,  #periphery to all downreg
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)


totalgo_result <- enrichGO(gene = toptotal_entrez_ids, #total tumor to all upreg
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
totaldown_go_result <- enrichGO(gene = toptotaldown_entrez_ids,  #total tumorto all downreg
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2)

# Visualize

barplot(corego_result, showCategory = 10, title = "GO аналіз клітин кору пухлини (Top Upregulated)")
barplot(coredown_go_result, showCategory = 10, title = "GO аналіз клітин кору пухлини (Top Downregulated)")

barplot(perygo_result, showCategory = 10, title = "GO аналіз клітин периферії пухлини (Top Upregulated)")
barplot(perydown_go_result, showCategory = 10, title = "GO аналіз клітин периферії пухлини (Top Downregulated)")

barplot(totalgo_result, showCategory = 10, title = "GO аналіз всіх клітин пухлини (Top Upregulated)")
barplot(totaldown_go_result, showCategory = 10, title = "GO аналіз всіх клітин пухлини (Top Downregulated)")

#adding come clarity by extra visualisation

testobj$characteristics <- ifelse(testobj$Location == "Periphery" & testobj$combined_clusters == "GBM_cluster", 
                                  "periphery",
                                  ifelse(testobj$Location == "Tumor" & testobj$combined_clusters == "GBM_cluster", 
                                         "core", 
                                         "other"))

Idents(testobj) <- "characteristics"

DimPlot(testobj, reduction = "umap", group.by = "characteristics")+ 
  ggtitle("Локалізація кластерів досліджуваних клітин") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 25, hjust = -0.5), margin = margin(b = 10))

ident_counts <- table(testobj$characteristics)
