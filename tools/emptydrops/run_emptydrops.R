library(DropletUtils)
library(Matrix)
library(dplyr)
library(ggplot2)

THRESHOLD = 0.001 # FDR of 0.1% for droplet calling

args = commandArgs(trailingOnly=TRUE)
datalocation <- args[1]
outdir <- args[2]
profile <- args[3]
datatype <- args[4]

# 
# doED <- function(sce, outdir) {
# 
#   bcrank <- barcodeRanks(counts(sce))
#   png(paste0(outdir, '/ranktest.png'), width = 6, height = 6, units = 'in', res=300)
#   plot(bcrank$rank, bcrank$total, log="xy", xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
# 
#   abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
#   abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
# 
#   legend("bottomleft", legend=c("Inflection", "Knee"),col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
#   dev.off()
# 
#   ed.out <- emptyDrops(counts(sce))
#   print(head(ed.out))
#   print('....')
#   print(tail(ed.out))
#   cells <- ed.out %>% as.data.frame() %>% filter(!is.na(FDR)) %>% filter(FDR < THRESHOLD) %>% filter(Total > 500) %>% rownames(.)
# 
#   print(head(cells))
#   print('____')
#   print(tail(cells))
# 
# 
#   cat(cells, file = paste0(outdir, '/emptydrops_droplets.txt'), sep = '\n')
# }
# 

doED <- function(in.sce, barcodes, threshold, outdir) {
  
  # bcrank <- barcodeRanks(counts(in.sce))
  print(head(counts(in.sce)))
  print(dim(counts(in.sce)))
  ed.out <- emptyDrops(counts(in.sce))
  
  ed.out$bc <- barcodes
  ed.out$iscell <- ifelse(ed.out$FDR < threshold, T, F)
  
  write.csv(ed.out, file = paste0(outdir, '/edout.csv'))
  
  cells <- ed.out %>% as.data.frame() %>% filter(iscell) %>% pull(bc)
  cat(cells, file = paste0(outdir, '/emptydrops_droplets.txt'), sep = '\n')
}

if (profile == 'Multiome'){
  features = data.table::fread(paste0(datalocation, '/outs/raw_feature_bc_matrix/features.tsv.gz'))
  barcodes = data.table::fread(paste0(datalocation, '/outs/raw_feature_bc_matrix/barcodes.tsv.gz'), header=F)$V1
  sce <- read10xCounts(paste0(datalocation, '/outs/raw_feature_bc_matrix/'))
  colnames(sce) <- barcodes
  
  if (datatype == 'rna') {
    sce.gex <- sce[features %>% filter(V3 == 'Gene Expression') %>% pull(V1), ]
    doED(sce.gex, barcodes, THRESHOLD, outdir)
  } else if (datatype == 'atac') {
    sce.atac <- sce[features %>% filter(V3 == 'Peaks') %>% pull(V1), ]
    doED(sce.atac, barcodes, THRESHOLD, outdir)
  }
} else if (profile == 'RNA') {
  print('RNA only')
  sce <- read10xCounts(paste0(datalocation, '/outs/raw_feature_bc_matrix/'))
  colnames(sce) <- data.table::fread(paste0(datalocation, '/outs/raw_feature_bc_matrix/barcodes.tsv.gz'), header=F)$V1
  # doED(sce, outdir)
  doED(sce, colnames(sce), THRESHOLD, outdir)


} else if (profile == 'ATAC') {
  print('ATAC only')
  Matrix::readMM(paste0(datalocation, '/outs/raw_peak_bc_matrix/matrix.mtx')) -> mtx
  barcodes <- read.table(paste0(datalocation, '/outs/raw_peak_bc_matrix/barcodes.tsv'), header=F)$V1
  features <- data.table::fread(cmd = paste0('cat ', datalocation, '/outs/raw_peak_bc_matrix/peaks.bed | tr "\t" "_" '), header=F)$V1

  colnames(mtx) <- barcodes
  rownames(mtx) <- features

  sce <- SingleCellExperiment(assays = list(counts = mtx))
  # doED(sce, outdir)
  doED(sce, barcodes, THRESHOLD, outdir)
  
}
