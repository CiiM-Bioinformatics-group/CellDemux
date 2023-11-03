try(dev.off())
rm(list = ls())

library(dplyr)
library(magrittr)
library(ggplot2)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

outlocation = args[1]

fread(paste0(outlocation, '/tmp/tmpvcfstats'), fill = T) %>% as.data.frame() -> sc
fread(paste0(outlocation, '/tmp/tmprefvcfstats'), fill = T) %>% as.data.frame() -> ref

refsamples <- read.csv(paste0(outlocation, '/tmp/tmprefvcfsamples'), header=F) %>% pull(V1)
sc_samples <- read.csv(paste0(outlocation, '/tmp/tmpvcfsamples'), header=F) %>% pull(V1)

# Check if we have "chr" in font of the chromosomes. If not, add
if(!any(grepl('chr', sc$CHROM))) { sc$CHROM <- paste0('chr', sc$CHROM) }
if(!any(grepl('chr', ref$CHROM))) { ref$CHROM <- paste0('chr', ref$CHROM)}

sc$newID <- paste0(sc$CHROM, '_', sc$POS, '_', sc$REF, '_', sc$ALT)
ref$newID <- paste0(ref$CHROM, '_', ref$POS, '_', ref$REF, '_', ref$ALT)

# We match on chr, bp, ref and alt
overlaps <- intersect(sc$newID, ref$newID)
print(paste0('Nr of overlapping SNPs based on CHR, POS, REF and ALT: ', length(overlaps)))

# We retain only the SNPs that match based on CHR, POS, REF and ALT.
# Take care of SNPS that
sc %<>% filter(newID %in% overlaps) %>% distinct(newID, .keep_all = T) %>% arrange(match(newID, overlaps))
ref %<>% filter(newID %in% overlaps) %>% distinct(newID, .keep_all = T) %>% arrange(match(newID, overlaps))
stopifnot(all(sc$newID == ref$newID))

# Replace "/" to "|" if any to make equal
sc[, as.character(sc_samples)] <- apply(X = sc[, as.character(sc_samples)], MARGIN = 2, gsub, pattern = '/', replacement = '|')
ref[, refsamples] <- apply(ref[, refsamples], 2, gsub, pattern = '/', replacement = '|')


doSum <- function(x) {
  # Function to collapse the 1|1, 0|1 notations to the sum of alternate alleles
  strsplit(x = x, split = '|')[[1]] -> alleles
  ref = alleles[[1]]
  alt = alleles[[3]]
  if(ref == '.' | alt == '.'){
    return(NA)
  } else {
    return(sum(as.numeric(ref), as.numeric(alt)))
  }
}


# Matching
counts <- matrix(0, nrow = length(refsamples), ncol = length(sc_samples), dimnames = list(refsamples, sc_samples))

for (sc_cluster in sc_samples) {
  for (refsample in refsamples) {
# 
#     # Compare all SNPs at the same time
#     # By comparing this way, all genotypes that are equal between reference and
#     # souporcell are TRUE (=1), those that dont match are false (=0).
#     # Taking the sum gives us the sum of TRUEs = matching genotypes
#     # SNPs with missing information are included in the sc_genotype, and will produce FALSE because they dont match.
#     sc_genotype <- sc %>% pull(as.character(sc_cluster))
#     ref_genotype <- ref  %>% pull(refsample)
# 
#     counts[refsample, as.character(sc_cluster)] <- sum(sc_genotype == ref_genotype)
# 
#     
#     
    sc_genotype <- sc %>% pull(as.character(sc_cluster))
    ref_genotype <- ref  %>% pull(refsample)
    
    scgen <- lapply(sc_genotype, doSum) %>% unlist()
    refgen <- lapply(ref_genotype, doSum) %>% unlist()
    
    counts[refsample, as.character(sc_cluster)] <- sum(na.omit(scgen == refgen))
    
  }
}

reshape2::melt(counts) -> plotdf
plotdf %>% group_by(Var2) %>% arrange(desc(value)) %>% slice(1) -> label

# Grubss test
library(outliers)

doGrubss <- function(v) {
  if (var(v) == 0) {return(list('p' = NA, 'sample' = NA, 'sig' = NA))} # If there is nog variation in the nr of SNPs being tested whatsoever..
  
  highest.sample <- sort(v, decreasing = T) %>% head(1) %>% names()

  # We need to check if we need to apply the opposite argument.
  mean(v) -> m
  
  lapply(v, function(x) {abs(x - m)} ) %>% do.call(rbind, .) %>%
    as.data.frame() %>% set_colnames(c('dev')) %>%
    mutate(test = ifelse(dev == max(dev), 'min', 'other')) %>%
    filter(test == 'min') %>% rownames(.) -> tested.sample
  
  if (length(tested.sample > 1)) {
    tested.sample <- tested.sample[1]
    print('Warning: multiple samples with same nr of SNPs. Testing random first one')
  }
  
  if(tested.sample != highest.sample) {
    opposite = T
  } else {
    opposite = F
  }

  x <- grubbs.test(v, type=10, opposite = opposite)

  return(list(
    'p' = x$p.value,
    'sample' = names(x$p.value),
    'sig' = ifelse(x$p.value < 0.05, T, F)
  ))
}

apply(counts, 2, doGrubss) %>% do.call(rbind, .) %>% as.data.frame() %>% mutate(cluster = rownames(.)) %>% apply(., 2, unlist) %>% as.data.frame() -> grubbs.res
write.csv(x = counts, file = paste0(outlocation, '/counts_matching.csv'), quote=F, row.names = T, col.names = T)
write.csv(x=grubbs.res, file = paste0(outlocation, '/statistics_matching.csv'), quote = F, row.names=F)

pdf(paste0(outlocation, '/match_refgen.pdf'), width = 6, height = 3)
ggplot(plotdf) +
  geom_boxplot(aes(x=as.character(Var2), y=value), outlier.shape=NA) +
  geom_point(aes(x=as.character(Var2), y=value), position = position_jitter(width = 0.1, height = 0)) +
  labs(x = 'Souporcell', y = 'Nr. SNPs') +
  theme_bw() +
  ggrepel::geom_text_repel(data=label,aes(x=as.character(Var2),y=value, label=Var1))
dev.off()




#### To integrate:

# counts2 <- matrix(0, nrow = length(refsamples), ncol = length(sc_samples), dimnames = list(refsamples, sc_samples))


# for (sc_cluster in sc_samples) {
#   for (refsample in refsamples) {
#     
#     # Compare all SNPs at the same time
#     # By comparing this way, all genotypes that are equal between reference and
#     # souporcell are TRUE (=1), those that dont match are false (=0).
#     # Taking the sum gives us the sum of TRUEs = matching genotypes
#     # SNPs with missing information are included in the sc_genotype, and will produce FALSE because they dont match.
#     sc_genotype <- sc %>% pull(as.character(sc_cluster))
#     ref_genotype <- ref  %>% pull(refsample)
#     scgen <- lapply(sc_genotype, doSum) %>% unlist()
#     refgen <- lapply(ref_genotype, doSum) %>% unlist()
#     counts2[refsample, as.character(sc_cluster)] <- sum(na.omit(scgen == refgen))
#     
#     # counts2[refsample, as.character(sc_cluster)] <- sum(sc_genotype == ref_genotype)
#     
#   }
# }
