library(SoupX)
library(Seurat)
args = commandArgs(trailingOnly=TRUE)

# SoupX takes only GEX by default. No need to filter
sc <- load10X(paste0(args[1], '/outs/'))

png(paste0(args[2], 'est_ambientrna.png'), height = 5, width = 8, units = 'in', res = 300)
sc = autoEstCont(sc, forceAccept = T, soupQuantile = 0.1, tfidfMin = 0.0)
dev.off()

cat(sc$fit$rhoEst, file = paste0(args[2], 'rho_estimate.txt'))