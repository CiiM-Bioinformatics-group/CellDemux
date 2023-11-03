library(dplyr)
library(magrittr)
library(Seurat)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
modality <- args[3]

print('Parameters: ')
print(paste0('input: ', input))
print(paste0('output: ', output))
print(paste0('Modality: ', modality))

ReadMtx(
  mtx = paste0(input, '/outs/raw_feature_bc_matrix/matrix.mtx.gz'),
  features = paste0(input, '/outs/raw_feature_bc_matrix/features.tsv.gz'),
  cells =  paste0(input, '/outs/raw_feature_bc_matrix/barcodes.tsv.gz')) -> data

features <- data.table::fread(paste0(input, '/outs/raw_feature_bc_matrix/features.tsv.gz')) -> features
barcodes <- data.table::fread(paste0(input, '/outs/raw_feature_bc_matrix/barcodes.tsv.gz'), header=F)$V1 -> barcodes

features <- features %>% filter(V3 == modality)
data <- data[features$V2, ]

cat(barcodes, sep = '\n', file = paste0(output, 'barcodes.tsv'))
writeMM(data, file = paste0(output, 'matrix.mtx'))
write.table(features, file = paste0(output, 'genes.tsv'), sep = '\t', quote=F, row.names = F, col.names = F)
