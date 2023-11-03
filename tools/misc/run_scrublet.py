#!/usr/bin/env python

import argparse
import scipy.io
import numpy as np
import scrublet as scr
import pandas as pd

parser = argparse.ArgumentParser(description = 'Run scrublet')
parser.add_argument('-m', required=True, help = 'Matrix.mtx.gz')
parser.add_argument('-g', required=True, help = 'features.tsv')
parser.add_argument('-b', required=True, help = 'barcodes.tsv.gz')
parser.add_argument('-o', required=True, help = 'outdir')
args = parser.parse_args()

# Code from: https://github.com/swolock/scrublet/blob/master/examples/scrublet_basics.ipynb
counts_matrix = scipy.io.mmread(args.m).T.tocsc()
genes = np.array(scr.load_genes(args.g, delimiter='\t', column=1))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)

# Export scrublet results
df = pd.DataFrame({
    'doublet_score': scrub.doublet_scores_obs_,
    'predicted_doublet': scrub.predicted_doublets_
})
df.to_csv("{}/scrublet_output.csv".format(args.o), index=False)

# Invert doublets to singlets
singlets = [not elem for elem in predicted_doublets]

# Read the barcodes, subset and export to a list of singlets
bc = pd.read_csv(args.b, header=None)
bc = bc[singlets]
bc.to_csv("{}/singlet_barcodes.txt".format(args.o), header=False, index=False)
