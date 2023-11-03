try(dev.off())
rm(list = ls())

library(DoubletFinder)
library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
datalocation <- args[1]
outdir <- args[2]
# profile <- args[3]

# load the data
# seur <- Read10X(paste0(datalocation, '/outs/filtered_feature_bc_matrix/'))

seur <- Read10X(datalocation)

# if (profile == 'Multiome') {
#   seur <- seur$`Gene Expression`
# }

seur <- CreateSeuratObject(seur)

seur %<>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:10) %>%
  FindNeighbors(reduction = 'pca', dims = 1:20) %>%
  FindClusters(resolution = 0.4)

# DoubletFinder wants reasonably cleaned data.
# Set some default parameters to clean the very obvious doublets
seur <- subset(seur, nFeature_RNA > 100 & nFeature_RNA < 3500 & nCount_RNA < 7500)

sweep.res <- paramSweep_v3(seur)
sweep.summary <- summarizeSweep(sweep.res, GT=F)
res <- find.pK(sweep.summary)
res %>% arrange(desc(BCmetric)) %>% head(n=1) %>% pull(pK) %>% as.character() %>% as.numeric() -> pK.choose

pdf(paste0(outdir, '/dis_pK_BCmetric.pdf'), width = 5, height = 2)
ggplot(res, mapping = aes(x = pK, y = BCmetric)) +
  geom_point() +
  geom_line(mapping = aes(group =1)) +
  theme_classic() +
  theme(aspect.ratio = .5) +
  geom_vline(xintercept = pK.choose, lty=2, col='red')
dev.off()

seur <- doubletFinder_v3(seur,
                         PCs = 1:10,
                         pN = 0.25,
                         pK = pK.choose,
                         nExp = modelHomotypic(seur$seurat_clusters),
                         reuse.pANN = FALSE, sct = FALSE)

colnames(seur@meta.data)[grepl('DF.classification', colnames(seur@meta.data))] -> colname
seur = seur[, seur@meta.data[, colname] == "Singlet"]

cat(colnames(seur), file = paste0(outdir, '/singlet_barcodes.txt'), sep = '\n')
