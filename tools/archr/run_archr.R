install.packages('Seurat',repos = "http://cran.us.r-project.org")


library(ArchR)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
fragments <- args[1]
barcodes <- args[2]
pool <- args[3]
outdir <- args[4]

addArchRGenome("hg38")
addArchRThreads(threads=1)

barcodes = read.csv(barcodes, header=F)$V1

print("Fragments file: ")
print(fragments)

ArrowFiles <- createArrowFiles(inputFiles = fragments,
                               sampleNames = pool,
                               outputNames = pool,
                               QCDir = "QualityControl",
                               threads = 1,
                               addTileMat=T,
                               addGeneScoreMat=F,
                               subThreading = F,
                               validBarcodes = barcodes,
                               force = T, cleanTmp = F)

doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 10,
                               knnMethod = "UMAP",
                               LSIMethod = 1)


bcs <- names(doubScores[[1]]$doubletEnrich[which(doubScores[[1]]$doubletEnrich < 1)])
lapply(bcs, function(x) {str_split(x, pattern = '#')[[1]][2]}) %>% unlist() -> bcs

cat(bcs, file = paste0(outdir, 'archr_singlets.txt'), sep = '\n')
