import yaml
import os
import glob
from itertools import product
import json
import sys
import re
import pandas as pd

# Jobs that dont need to be submitted to the cluster
localrules: prep_files, define_confident_nonempty_ATAC

############################################################################################
# User parameters:
reffasta = ''       # Path to the CellRanger reference fasta file that the reads were mapped to
reffastafai = ''    # Path to the fasta index (.fai) file of the above fasta file
outlocation = ""    # Path to the location of the output
samplesheet = ""    # Path to the excel input sheet. See README for details.

############################################################################################

# Samples to demultiplex with the right information
df = pd.read_excel(samplesheet)

SAMPLES = df['pool']
DATALOCATION = dict(zip(df.loc[:, 'pool'], df.loc[:, 'datalocation']))
DONORS = dict(zip(df.loc[:, 'pool'], df.loc[:, 'donors']))
VCF = dict(zip(df.loc[:, 'pool'], df.loc[:, 'vcf']))
PROFILE = dict(zip(df.loc[:, 'pool'], df.loc[:, 'profile']))

# Make an output directory for every pool
for sample in df['pool']:
    outdir = f'{outlocation}/{sample}'
    if not os.path.exists(outdir): os.makedirs(outdir)

def get_output(df):
    out = []

    # Symlinked files
    out.extend(expand(outlocation + '{sample}/inputfiles/done', sample = SAMPLES))

    for index, row in df.iterrows():

        profile = row['profile']
        sample = row['pool']


        if profile == 'RNA':

            # RNA tools output
            out.extend([outlocation + sample + '/rna/qc/soupx/est_ambientrna.png'])
            out.extend([outlocation + sample + '/rna/vireo/donor_ids.tsv'])
            out.extend([outlocation + sample + '/rna/demuxlet/popscle_demuxlet.best'])
            out.extend([outlocation + sample + '/rna/scrublet/singlet_barcodes.txt'])
            out.extend([outlocation + sample + '/rna/doubletfinder/singlet_barcodes.txt'])

            # # Souporcell
            out.extend([outlocation + sample + '/rna/souporcell/clusters.tsv'])

            # CellBender
            out.extend([outlocation + sample + '/rna/qc/cellbender/cellbenderout_cell_barcodes.csv'])

            # EmtyDrops
            out.extend([outlocation + sample + '/rna/qc/emptydrops/emptydrops_droplets.txt'])

            #Overlap cellbender + emptydrops
            out.extend([outlocation + sample + '/rna/qc/manually_flt_data/barcodes.tsv'])

            # Confident singlets
            out.extend([outlocation + sample + '/rna/confident_rna_singlets.txt'])

            # Singlets
            out.extend([outlocation + sample + '/rna/singlet_data_only/flt_bam.bam'])

            # Final demultiplexing results
            out.extend([outlocation + sample + '/rna/demultiplex/clusters.tsv'])

            # Matching to reference genotype
            out.extend([outlocation + sample + '/rna/demultiplex/match_refgen.pdf'])

        elif profile == 'ATAC':

            # AMULET output
            out.extend([outlocation + sample + '/atac/amulet/singlet_barcodes.txt'])

            # # Souporcell
            out.extend([outlocation + sample + '/atac/souporcell/clusters.tsv'])

            # ArchR doublets
            out.extend([outlocation + sample + '/atac/archr/archr_singlets.txt'])

            # Confident singlets
            out.extend([outlocation + sample + '/atac/confident_atac_singlets.txt'])

            # Singlets
            out.extend([outlocation + sample + '/atac/singlet_data_only/flt_bam.bam'])

            # Matching to reference genotype
            out.extend([outlocation + sample + '/atac/souporcell/match_refgen.pdf'])

        elif profile == 'Multiome':

            out.extend([outlocation + sample + '/rna/qc/soupx/est_ambientrna.png'])
            out.extend([outlocation + sample + '/rna/vireo/donor_ids.tsv'])
            out.extend([outlocation + sample + '/rna/demuxlet/popscle_demuxlet.best'])
            out.extend([outlocation + sample + '/rna/scrublet/singlet_barcodes.txt'])
            out.extend([outlocation + sample + '/rna/doubletfinder/singlet_barcodes.txt'])

            # AMULET output
            out.extend([outlocation + sample + '/atac/amulet/singlet_barcodes.txt'])

            # CellBender
            out.extend(expand(outlocation + sample + '/{datatype}/qc/cellbender/cellbenderout_cell_barcodes.csv', datatype=['rna']))

            # Souporcell
            out.extend(expand(outlocation + sample + '/{datatype}/souporcell/clusters.tsv', datatype=['rna', 'atac']))

            # EmtyDrops
            out.extend(expand(outlocation + sample + '/{datatype}/qc/emptydrops/emptydrops_droplets.txt', datatype=['rna']))

            # Overlap cellbender + emptydrops
            # out.extend(expand(outlocation + sample + '/{datatype}/qc/manually_flt_data/barcodes.tsv', datatype=['rna', 'atac']))
            out.extend(expand(outlocation + sample + '/{datatype}/qc/manually_flt_data/barcodes.tsv', datatype=['rna']))

            # ArchR doublets
            out.extend([outlocation + sample + '/atac/archr/archr_singlets.txt'])

            # Confident singlets
            out.extend([outlocation + sample + '/confident_multiome_singlets.txt'])

            # Singlets
            out.extend(expand(outlocation + sample + '/rna/singlet_data_only/flt_bam.bam', datatype = ['rna', 'atac']))

            # Final demultiplexing results
            out.extend(expand(outlocation + sample + '/{datatype}/demultiplex/clusters.tsv', datatype=['rna', 'atac']))

            # Matching to reference genotype
            out.extend(expand(outlocation + sample + '/{datatype}/demultiplex/match_refgen.pdf', datatype=['rna', 'atac']))

        else:
            print('Wrong profile')
            sys.exit(1)
    return out

# Actual pipeline
rule all:
    input:
        get_output(df)

rule prep_files:
    output:
        outlocation + '{sample}/inputfiles/done'
    params:
        datalocation = lambda wcs: DATALOCATION[wcs.sample],
        vcf = lambda wcs: VCF[wcs.sample],
        profile = lambda wcs: PROFILE[wcs.sample],
    shell:
        """
        # If already exists, delete
        if [ -d {outlocation}/{wildcards.sample}/inputfiles ]; then
            rm -rf {outlocation}/{wildcards.sample}/inputfiles
        fi

        mkdir -p {outlocation}/{wildcards.sample}/inputfiles

        # Vcf needed for all profiles
        ln -s {params.vcf} {outlocation}{wildcards.sample}/inputfiles/refvcf.vcf

        if [ {params.profile} = 'Multiome' ]; then
            ln -s {params.datalocation}outs/atac_possorted_bam.bam {outlocation}{wildcards.sample}/inputfiles/atacbam.bam
            ln -s {params.datalocation}outs/atac_possorted_bam.bam.bai {outlocation}{wildcards.sample}/inputfiles/atacbam.bam.bai

            ln -s {params.datalocation}outs/gex_possorted_bam.bam {outlocation}{wildcards.sample}/inputfiles/rnabam.bam
            ln -s {params.datalocation}outs/gex_possorted_bam.bam.bai {outlocation}{wildcards.sample}/inputfiles/rnabam.bam.bai

            ln -s {params.datalocation}outs/per_barcode_metrics.csv {outlocation}{wildcards.sample}/inputfiles/metrics.csv

            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/matrix.mtx.gz {outlocation}{wildcards.sample}/inputfiles/matrix.mtx.gz
            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/features.tsv.gz {outlocation}{wildcards.sample}/inputfiles/features.tsv.gz
            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/barcodes.tsv.gz {outlocation}{wildcards.sample}/inputfiles/barcodes.tsv.gz
            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/barcodes.tsv {outlocation}{wildcards.sample}/inputfiles/barcodes.tsv

        elif [ {params.profile} = 'RNA' ]; then
            ln -s {params.datalocation}outs/possorted_genome_bam.bam {outlocation}{wildcards.sample}/inputfiles/rnabam.bam
            ln -s {params.datalocation}outs/possorted_genome_bam.bam.bai {outlocation}{wildcards.sample}/inputfiles/rnabam.bam.bai

            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/matrix.mtx.gz {outlocation}{wildcards.sample}/inputfiles/matrix.mtx.gz
            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/features.tsv.gz {outlocation}{wildcards.sample}/inputfiles/features.tsv.gz
            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/barcodes.tsv.gz {outlocation}{wildcards.sample}/inputfiles/barcodes.tsv.gz
            ln -s {params.datalocation}outs/filtered_feature_bc_matrix/barcodes.tsv {outlocation}{wildcards.sample}/inputfiles/barcodes.tsv

        elif [ {params.profile} = 'ATAC' ]; then
            ln -s {params.datalocation}outs/possorted_bam.bam {outlocation}{wildcards.sample}/inputfiles/atacbam.bam
            ln -s {params.datalocation}outs/possorted_bam.bam.bai {outlocation}{wildcards.sample}/inputfiles/atacbam.bam.bai

            ln -s {params.datalocation}outs/singlecell.csv {outlocation}{wildcards.sample}/inputfiles/metrics.csv

            ln -s {params.datalocation}outs/filtered_peak_bc_matrix/matrix.mtx.gz {outlocation}{wildcards.sample}/inputfiles/matrix.mtx.gz
            ln -s {params.datalocation}outs/filtered_peak_bc_matrix/peaks.bed.gz {outlocation}{wildcards.sample}/inputfiles/features.tsv.gz
            ln -s {params.datalocation}outs/filtered_peak_bc_matrix/barcodes.tsv.gz {outlocation}{wildcards.sample}/inputfiles/barcodes.tsv.gz
            ln -s {params.datalocation}outs/filtered_peak_bc_matrix/barcodes.tsv {outlocation}{wildcards.sample}/inputfiles/barcodes.tsv
        else
            echo "Wrong profile?"
            exit 1
        fi

        # Output file for snakemake
        touch {outlocation}{wildcards.sample}/inputfiles/done
        """

rule emptydrops:
    input: outlocation + '{sample}/inputfiles/done',
    output:
        outlocation + '{sample}/{datatype}/qc/emptydrops/emptydrops_droplets.txt'
    params:
        profile = lambda wcs: PROFILE[wcs.sample],
        datalocation = lambda wcs: DATALOCATION[wcs.sample]
    threads: 1
    resources: mem='15G', time='2:00:00'
    shell:  "mkdir -p {outlocation}{wildcards.sample}/{wildcards.datatype}/qc/emptydrops; Rscript tools/emptydrops/run_emptydrops.R {params.datalocation} {outlocation}{wildcards.sample}/{wildcards.datatype}/qc/emptydrops/ {params.profile} {wildcards.datatype}"

rule define_confident_nonempty_ATAC:
    input:
        outlocation + '{sample}/inputfiles/done'
    params:
        outdir = outlocation + '{sample}/atac/qc/manually_flt_data',
        bam = outlocation + '{sample}/inputfiles/atacbam.bam',
        barcodes = outlocation + '{sample}/inputfiles/barcodes.tsv',
        features = outlocation + '{sample}/inputfiles/features.tsv.gz'
    output:
        bam = outlocation + '{sample}/atac/qc/manually_flt_data/atacbam.bam',
        barcodes = outlocation + '{sample}/atac/qc/manually_flt_data/barcodes.tsv',
        features = outlocation + '{sample}/atac/qc/manually_flt_data/features.tsv.gz'
    shell:
        """
        # We rely on 10X's ATACseq cell calling algorithms.
        [ -d {params.outdir} ] && rm -rf {params.outdir}
        mkdir -p {params.outdir}

        ln -s {params.bam} {output.bam}
        ln -s {params.bam}.bai {output.bam}.bai

        ln -s {params.barcodes} {output.barcodes}
        ln -s {params.features} {output.features}
        """

rule define_confident_nonempty_RNA:
    input:
        ed_nonempty = outlocation + '{sample}/rna/qc/emptydrops/emptydrops_droplets.txt',
        cb_nonempty = outlocation + '{sample}/rna/qc/cellbender/cellbenderout_cell_barcodes.csv'
    params:
        outdir = outlocation + '{sample}/rna/qc/manually_flt_data/',
        datalocation = lambda wcs: DATALOCATION[wcs.sample],
        profile = lambda wcs: PROFILE[wcs.sample],
        complete_bam = outlocation + '{sample}/inputfiles/rnabam.bam'
    output:
        flt_bam = outlocation + '{sample}/rna/qc/manually_flt_data/bam.bam',
        barcodes = outlocation + '{sample}/rna/qc/manually_flt_data/barcodes.tsv',
        features = outlocation + '{sample}/rna/qc/manually_flt_data/genes.tsv',
        matrix = outlocation + '{sample}/rna/qc/manually_flt_data/matrix.mtx',
        outdir = directory(outlocation + '{sample}/rna/qc/manually_flt_data/')
    # threads: 15
    resources: mem='25G', time='10:00:00', cpus_per_task=15
    shell:
        """
        mkdir -p {outlocation}{wildcards.sample}/rna/qc/manually_flt_data/

        # Barcodes that we think are non-empty
        comm -12 <(sort {input.ed_nonempty}) <(sort {input.cb_nonempty}) > {outlocation}{wildcards.sample}/rna/qc/overlap_ed_cb_droplets.txt

        # Read in the raw feature matrix, subset the right barcodes for other tools to use
        # If exists, first delete to stop DropletUtils from complaining
        if test -d {params.outdir}; then rm -r {params.outdir}; fi

        Rscript tools/flt_feature_bc_matrix.R {params.datalocation} {outlocation}{wildcards.sample}/rna/qc/overlap_ed_cb_droplets.txt rna {params.profile} {params.outdir}

        # Test if filtered bam exists and remove if so
        if test -f {output.flt_bam}; then rm {output.flt_bam}; fi

        # Subset the bam file and index
        subset-bam --bam {params.complete_bam} \
            --cell-barcodes {output.barcodes} \
            --out-bam {output.flt_bam} \
            --cores 15
        samtools index {output.flt_bam}

        # Barcodes.tsv is the file we use, this is superfluous now.
        rm {outlocation}{wildcards.sample}/rna/qc/overlap_ed_cb_droplets.txt
        """

rule define_confident_singlets_ATAC:
    input:
        amulet_singlets = outlocation + '{sample}/atac/amulet/singlet_barcodes.txt',
        archr_singlets = outlocation + '{sample}/atac/archr/archr_singlets.txt'
    output:
        singlets = outlocation + '{sample}/atac/confident_atac_singlets.txt'
    shell:
        """
        # We take as ATAC confident singlets the singlets predicted by Amulet
        ln -s {input.amulet_singlets} {output.singlets}
        """

rule define_confident_singlets_RNA:
    input:
        souporcell_singlets = outlocation + '{sample}/rna/souporcell/clusters.tsv',
        scrublet_singlets = outlocation + '{sample}/rna/scrublet/singlet_barcodes.txt',
        doubletfinder_singlets = outlocation + '{sample}/rna/doubletfinder/singlet_barcodes.txt'
    output:
        singlets = outlocation + '{sample}/rna/confident_rna_singlets.txt'
    shell:
        """
        # We take as RNA confident singlets anything that is predicted either by Scrublet or Souporcell or DoubletFinder

        Rscript tools/filter_rna_confident_singlets.R {wildcards.sample} {outlocation}
        """

rule define_confident_singlets_multiome:
    input:
        amulet_singlets = outlocation + '{sample}/atac/amulet/singlet_barcodes.txt',
        souporcell_singlets = outlocation + '{sample}/rna/souporcell/clusters.tsv',
        scrublet_singlets = outlocation + '{sample}/rna/scrublet/singlet_barcodes.txt',
        doubletfinder_singlets = outlocation + '{sample}/rna/doubletfinder/singlet_barcodes.txt'
    output:
        singlets = outlocation + '{sample}/confident_multiome_singlets.txt'
    shell:
        """
        Rscript tools/filter_multiome_confident_singlets.R {wildcards.sample} {outlocation}
        """

rule soupx:
    output:
        est_ambientrna = outlocation + '{sample}/rna/qc/soupx/est_ambientrna.png'
    params:
        datalocation = lambda wcs: DATALOCATION[wcs.sample]
    threads: 1
    resources: mem='15G', time='2:00:00'
    shell:  "mkdir -p '{outlocation}{wildcards.sample}/rna/qc/soupx'; Rscript tools/soupx/run_soupx.R {params.datalocation} '{outlocation}{wildcards.sample}/rna/qc/soupx/'"

rule cellbender:
    input: outlocation + '{sample}/inputfiles/done',
    output:
        outlocation + '{sample}/{datatype}/qc/cellbender/cellbenderout_cell_barcodes.csv'
    params:
        datalocation = lambda wcs: DATALOCATION[wcs.sample],
        profile = lambda wcs: PROFILE[wcs.sample],
        modality = lambda wcs: "Gene Expression" if wcs.datatype == 'rna' else 'Peaks'
    resources:
        partition='gpu', extra_options='--gres=gpu:t4:1', mem='35G', time='15:00:00'
    shell:
        """

        mkdir -p '{outlocation}{wildcards.sample}/{wildcards.datatype}/qc/cellbender'

        if [ {params.profile} = 'Multiome' ]; then
            Rscript tools/cellbender/load_subsetmodality.R {params.datalocation} '{outlocation}{wildcards.sample}/{wildcards.datatype}/qc/cellbender/' '{params.modality}'
            input={outlocation}{wildcards.sample}/{wildcards.datatype}/qc/cellbender/
        elif [ {wildcards.datatype} = 'rna' ]; then
            input={params.datalocation}/outs/raw_feature_bc_matrix/
        elif [ {wildcards.datatype} = 'atac' ]; then
            ln -sf {params.datalocation}/outs/raw_peak_bc_matrix/peaks.bed {params.datalocation}/outs/raw_peak_bc_matrix/genes.tsv
            input={params.datalocation}/outs/raw_peak_bc_matrix/
        fi

        singularity exec --cleanenv --nv /vol/projects/BIIM/resources/tools/singularity/cellbender.sif cellbender remove-background \
          --input $input \
          --output {outlocation}{wildcards.sample}/{wildcards.datatype}/qc/cellbender/cellbenderout.h5 \
          --expected-cells 8000 \
          --cuda \
          --cells-posterior-reg-calc 10 \
          --posterior-batch-size 5 \
          --epochs 100
         """

rule archr:
    input:
        barcodes = outlocation + '{sample}/atac/qc/manually_flt_data/barcodes.tsv',
    output:
        singlets_archr = outlocation + '{sample}/atac/archr/archr_singlets.txt'
    conda: "archr"
    params:
        outdir = outlocation + '{sample}/atac/archr/',
        fragmentsfile = lambda wcs: DATALOCATION[wcs.sample] + '/outs/atac_fragments.tsv.gz' if PROFILE[wcs.sample] == 'Multiome' else DATALOCATION[wcs.sample] + '/outs/fragments.tsv.gz'
    threads: 1
    resources: mem = '15G', time='2:00:00'
    shell:  "mkdir -p {params.outdir}; Rscript tools/archr/run_archr.R {params.fragmentsfile} {input.barcodes} {wildcards.sample} {params.outdir}"

rule doubletfinder:
    input:
        indir = outlocation + '{sample}/rna/qc/manually_flt_data/'
    output:
        singlets_doubletfinder = outlocation + '{sample}/rna/doubletfinder/singlet_barcodes.txt'
    params:
        # datalocation = outlocation + '{sample}/inputfiles/',
        # profile = lambda wcs: PROFILE[wcs.sample],
        outdir = outlocation + '{sample}/rna/doubletfinder/'
    threads: 1
    resources: mem = '25G', time='1:00:00'
    shell:  "mkdir -p {params.outdir}; Rscript tools/doubletfinder/run_doubletfinder.R {input.indir} {params.outdir}"

rule amulet:
    input:
        bam = outlocation + '{sample}/atac/qc/manually_flt_data/atacbam.bam',
        # metrics = outlocation +'{sample}/inputfiles/metrics.csv',
        barcodes = outlocation +'{sample}/atac/qc/manually_flt_data/barcodes.tsv'
    output:
        singletbarcodes = outlocation + '{sample}/atac/amulet/singlet_barcodes.txt'
    params:
        outdir = outlocation + '{sample}/atac/amulet/',
        bcidx = 0,
        iscellidx = lambda wcs: 3 if PROFILE[wcs.sample] == 'Multiome' else 9,
        metrics = outlocation +'{sample}/inputfiles/metrics.csv'
    conda: "amulet"
    threads: 1
    resources: mem = '15G', time='3:00:00'
    shell:
        """
        mkdir -p {params.outdir}

        tools/AMULET/AMULET.sh {input.bam} \
            {params.metrics} \
            tools/AMULET/human_autosomes.txt \
            tools/AMULET/hg38-blacklist.v2.bed \
            {params.outdir} \
            tools/AMULET/ \
            --forcesorted \
            --bcidx {params.bcidx} \
            --iscellidx {params.iscellidx}

        # AMULET writes out only the doublets. Write file containing singlets to make life easier
        grep -f "{outlocation}{wildcards.sample}/atac/amulet/MultipletBarcodes_01.txt" <(cat {input.barcodes}) | grep -vf- <(cat {input.barcodes}) > "{outlocation}{wildcards.sample}/atac/amulet/singlet_barcodes.txt"
        """

rule map_minimap2:
    "First souporcell rule for RNA libraries. Prepares the reads to be mapped and maps using minimap2. Bam files are retagged, sorted and indexed"
    input:
        bam = outlocation + '{sample}/rna/qc/manually_flt_data/bam.bam',
        barcodes = outlocation + '{sample}/rna/qc/manually_flt_data/barcodes.tsv',
    output:
        # outdir = directory(outlocation + '{sample}/rna/souporcell'),
        maptagged = outlocation + '{sample}/rna/souporcell/maptagged_sorted.bam'
    params:
        outdir = directory(outlocation + '{sample}/rna/souporcell/')
    conda: "souporcell"
    resources: mem = '50G', time='80:00:00', qos='long', cpus_per_task = 15
    shell:
        """
        mkdir -p {params.outdir}

        tools/souporcell/renamer.py --bam {input.bam} --barcodes {input.barcodes} --out {params.outdir}/fq.fq

        # mapper conditional on what data we have.
        minimap2 -ax splice -t 15 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no {reffasta} {params.outdir}/fq.fq > {params.outdir}/map.sam

        tools/souporcell/retag.py --sam {params.outdir}/map.sam --out {params.outdir}/maptagged.bam
        samtools sort {params.outdir}/maptagged.bam {params.outdir}/maptagged_sorted
        samtools index {params.outdir}/maptagged_sorted.bam

        # Clean up
        rm {params.outdir}/fq.fq && rm {params.outdir}/map.sam && rm {params.outdir}/maptagged.bam
        """

rule map_bwamem:
    "Rule to remap reads using BWA-MEM. Only for ATAC pool"
    input:
        # singletbarcodes = rules.define_confident_singlets.output.outfile,
        singletbarcodes = outlocation + '{sample}/atac/confident_atac_singlets.txt',
        bam = outlocation + '{sample}/inputfiles/atacbam.bam'
    params:
        outdir = directory(outlocation + '{sample}/atac/souporcell')
    output:
        maptagged = outlocation + '{sample}/atac/souporcell/maptagged_sorted.bam'
    conda: "demultiplex_atac"
    resources: mem = '75G', time='40:00:00', cpus_per_task = 15
    shell:
        """
        mkdir -p {params.outdir}

        # Get the reads associated to singlet barcodes from AMULET.
        # Try to delete BAM file if it exists first so that subset-bam does not complain.
        if test -f {params.outdir}/filtered_bam_singlets.bam; then
         rm {params.outdir}/filtered_bam_singlets.bam
        fi

        subset-bam --bam {input.bam} \
          --cell-barcodes {input.singletbarcodes} \
          --out-bam {params.outdir}/filtered_bam_singlets.bam \
          --cores 15

        # Rename, remap and retag the reads
        tools/souporcell/renamer.py --no_umi True \
            --bam {params.outdir}/filtered_bam_singlets.bam \
            --barcodes {input.singletbarcodes} \
            --out {params.outdir}/fq.fq

        bwa mem {reffasta} {params.outdir}/fq.fq -t 15 > {params.outdir}/map.sam
        tools/souporcell/retag.py --no_umi True --sam {params.outdir}/map.sam --out {params.outdir}/maptagged.bam

        # Sort, index using samtools
        samtools sort {params.outdir}/maptagged.bam {params.outdir}/maptagged_sorted
        samtools index {params.outdir}/maptagged_sorted.bam

        # Clean up
        rm {params.outdir}/filtered_bam_singlets.bam && rm  {params.outdir}/fq.fq && rm {params.outdir}/map.sam && rm {params.outdir}/maptagged.bam
        """

rule callsnps:
    "Rule to call SNPs in parallel with Freebayes. Done on ATAC and RNA libraries alike"
    input:
        bam = outlocation + '{sample}/{datatype}/souporcell/maptagged_sorted.bam',
    output:
        freebayesvcf = outlocation + '{sample}/{datatype}/souporcell/freebayes.vcf'
    conda: "souporcell"
    resources: mem = '45G', time='60:00:00', qos='long', cpus_per_task = 25
    shell: "tools/freebayes/scripts/freebayes-parallel <(tools/freebayes/scripts/fasta_generate_regions.py {reffastafai} 100000) 25 -f {reffasta} -iXu -C 2 -q 30 -n 3 -E 1 -m 30 --min-coverage 6 {input.bam} > {output.freebayesvcf}"

rule run_souporcellmodel:
    "Run the Souporcell machine learning models to cluster cells based on called SNPs and remapped reads. Rule used for RNA and ATAC both"
    input:
        snps = rules.callsnps.output.freebayesvcf,
        bam = outlocation + '{sample}/{datatype}/souporcell/maptagged_sorted.bam',
        vcf = outlocation + '{sample}/inputfiles/refvcf.vcf'
    output:
        clusters = outlocation + '{sample}/{datatype}/souporcell/clusters.tsv',
        cluster_genotypes = outlocation + '{sample}/{datatype}/souporcell/cluster_genotypes.vcf'
    params:
        ndonors = lambda wcs: DONORS[wcs.sample],
        barcodes = lambda wcs: outlocation + wcs.sample + '/rna/qc/manually_flt_data/barcodes.tsv' if wcs.datatype == 'rna' else outlocation + wcs.sample +'/atac/amulet/singlet_barcodes.txt',
        umiflag = lambda wcs: "" if wcs.datatype == 'atac' else '--umi',
        outdir = directory(outlocation + '{sample}/{datatype}/souporcell')
    conda: "souporcell"
    resources: mem = '50G', time='50:00:00', qos='long', cpus_per_task = 25
    shell:
        """
        # To stop vartrix from complaining.
        if test -f {params.outdir}/ref.mtx; then rm {params.outdir}/ref.mtx; fi
        if test -f {params.outdir}/alt.mtx; then rm {params.outdir}/alt.mtx; fi

        vartrix --mapq 30 \
            {params.umiflag} \
            -b {input.bam} \
            -c {params.barcodes} \
            --scoring-method coverage \
            --threads 25 \
            --ref-matrix {params.outdir}/ref.mtx \
            --out-matrix {params.outdir}/alt.mtx \
            -v {input.snps} \
            --fasta {reffasta}

        tools/souporcell/souporcell.py --min_alt 6 --min_ref 6 --threads 25 \
            -a {params.outdir}/alt.mtx \
            -r {params.outdir}/ref.mtx \
            -b {params.barcodes} \
            -k {params.ndonors}  \
            --out {params.outdir}/clusters_tmp.tsv

        tools/souporcell/troublet/target/release/troublet -a {params.outdir}/alt.mtx \
            -r {params.outdir}/ref.mtx \
            --clusters {params.outdir}/clusters_tmp.tsv > {params.outdir}/clusters.tsv

        tools/souporcell/consensus.py -c {params.outdir}/clusters.tsv \
            -a {params.outdir}/alt.mtx \
            -r {params.outdir}/ref.mtx \
            --soup_out {params.outdir}/soup.txt \
            -v {params.outdir}/freebayes.vcf \
            --vcf_out {params.outdir}/cluster_genotypes.vcf \
            --output_dir {params.outdir}

        # Visualize souporcell results
        Rscript tools/souporcell/souporcell_pcaclusters.R {params.outdir}/clusters.tsv {params.outdir}/pca_souporcell_clusters.png
        """

rule match_genotype:
    "Match the genotypes of the souporcell clusters versus the reference genotype"
    input:
        cluster_genotypes = outlocation + '{sample}/{datatype}/demultiplex/cluster_genotypes.vcf',
        vcf = outlocation + '{sample}/inputfiles/refvcf.vcf'
    params:
        outdir = directory(outlocation + '{sample}/{datatype}/demultiplex/')
    output:
        matchref = outlocation + '{sample}/{datatype}/demultiplex/match_refgen.pdf'
    threads: 1
    resources: mem = '25G', time='1:00:00'
    shell:
        """
        # Requirements:
          # bcftools
          # Both vcf files assumed to contain the GT field

        mkdir -p {params.outdir}/tmp

        echo -e "ID\tCHROM\tPOS\tREF\tALT\t$(bcftools query -l {input.vcf} | tr "\n" "\t")" > {params.outdir}/tmp/tmprefvcfstats
        bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input.vcf} >> {params.outdir}/tmp/tmprefvcfstats
        bcftools query -l {input.vcf} > {params.outdir}/tmp/tmprefvcfsamples

        echo -e "ID\tCHROM\tPOS\tREF\tALT\t$(bcftools query -l {input.cluster_genotypes} | tr "\n" "\t")" > {params.outdir}/tmp/tmpvcfstats
        bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input.cluster_genotypes} >> {params.outdir}/tmp/tmpvcfstats
        bcftools query -l {input.cluster_genotypes} > {params.outdir}/tmp/tmpvcfsamples

        Rscript tools/souporcell/matchvcf.R {params.outdir}

        # Clean output
        rm -rf {params.outdir}/tmp
        """

rule match_genotypeATAC:
    "Match the genotypes of the souporcell clusters versus the reference genotype. ONLY USED FOR ATAC, TO REPLACE"
    input:
        cluster_genotypes = outlocation + '{sample}/{datatype}/souporcell/cluster_genotypes.vcf',
        vcf = outlocation + '{sample}/inputfiles/refvcf.vcf'
    params:
        outdir = directory(outlocation + '{sample}/{datatype}/souporcell/')
    output:
        matchref = outlocation + '{sample}/{datatype}/souporcell/match_refgen.pdf'
    threads: 1
    resources: mem = '25G', time='1:00:00'
    shell:
        """
        # Requirements:
          # bcftools
          # Both vcf files assumed to contain the GT field

        mkdir -p {params.outdir}/tmp

        echo -e "ID\tCHROM\tPOS\tREF\tALT\t$(bcftools query -l {input.vcf} | tr "\n" "\t")" > {params.outdir}/tmp/tmprefvcfstats
        bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input.vcf} >> {params.outdir}/tmp/tmprefvcfstats
        bcftools query -l {input.vcf} > {params.outdir}/tmp/tmprefvcfsamples

        echo -e "ID\tCHROM\tPOS\tREF\tALT\t$(bcftools query -l {input.cluster_genotypes} | tr "\n" "\t")" > {params.outdir}/tmp/tmpvcfstats
        bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input.cluster_genotypes} >> {params.outdir}/tmp/tmpvcfstats
        bcftools query -l {input.cluster_genotypes} > {params.outdir}/tmp/tmpvcfsamples

        Rscript tools/souporcell/matchvcf.R {params.outdir}

        # Clean output
        rm -rf {params.outdir}/tmp
        """

rule vireo:
    input:
        bam = outlocation + '{sample}/rna/qc/manually_flt_data/bam.bam',
        barcodes = outlocation + '{sample}/rna/qc/manually_flt_data/barcodes.tsv',
        vcf = outlocation + '{sample}/inputfiles/refvcf.vcf'
    output:
        donorids = outlocation + '{sample}/rna/vireo/donor_ids.tsv'
    threads: 15
    conda: "vireo"
    resources: mem='150G', time='40:00:00'
    shell:
        """
        mkdir -p "{outlocation}{wildcards.sample}/rna/vireo/"

        cellsnp-lite -s {input.bam} -b {input.barcodes} -O {outlocation}{wildcards.sample}/rna/vireo/ -R {input.vcf} -p 15 --minMAF 0.01 --minCOUNT 20 --gzip --UMItag UB
        vireo -c {outlocation}{wildcards.sample}/rna/vireo/ -d {input.vcf} -o {outlocation}{wildcards.sample}/rna/vireo/ -t GT
        """

rule demuxlet:
    input:
        bam = outlocation + '{sample}/rna/qc/manually_flt_data/bam.bam',
        barcodes = outlocation + '{sample}/rna/qc/manually_flt_data/barcodes.tsv',
        vcf = outlocation + '{sample}/inputfiles/refvcf.vcf'
    output:
        demuxbest = outlocation + '{sample}/rna/demuxlet/popscle_demuxlet.best'
    threads: 5
    conda: "demuxlet"
    resources: mem="150G", time='48:00:00'
    shell:
        """
        mkdir -p "{outlocation}{wildcards.sample}/rna/demuxlet/"

        source tools/popscle_tools.sh

        # Get the correct BCFtools plugins to get popscle_tools working.
        export BCFTOOLS_PLUGINS=/home/mzoodsma/bin/bcftools-1.16/plugins/

        # Rename the vcf to correct chr notation
        tools/bcftools annotate --rename-chrs tools/chr_names.txt {input.vcf} -Ov -o "{outlocation}{wildcards.sample}/rna/demuxlet/vcf_renamed.vcf"

        # Compress with bgzip to stop bcftools complaining
        bgzip -c "{outlocation}{wildcards.sample}/rna/demuxlet/vcf_renamed.vcf" > "{outlocation}{wildcards.sample}/rna/demuxlet/bgzip_renamed_vcf.vcf.gz"
        tools/bcftools index "{outlocation}{wildcards.sample}/rna/demuxlet/bgzip_renamed_vcf.vcf.gz"

        # Pre-process the vcf file
        tools/bcftools view "{outlocation}{wildcards.sample}/rna/demuxlet/bgzip_renamed_vcf.vcf.gz" -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 |
          only_keep_snps |
          filter_out_mutations_homozygous_reference_in_all_samples |
          filter_out_mutations_homozygous_in_all_samples |
          calculate_AF_AC_AN_values_based_on_genotype_info > "{outlocation}{wildcards.sample}/rna/demuxlet/final.vcf"

        # Pre-process bam file
        samtools view -b {input.bam} chr{{1..22}} > "{outlocation}{wildcards.sample}/rna/demuxlet/flt.bam"

        bash /vol/projects/CIIM/resources/tools/demuxlet/sort_demuxlet.sh "{outlocation}{wildcards.sample}/rna/demuxlet/flt.bam" "{outlocation}{wildcards.sample}/rna/demuxlet/final.vcf" > "{outlocation}{wildcards.sample}/rna/demuxlet/final_sorted.vcf"

        # New, filter bam file for speed?
        sh tools/popscle_filterbam.sh "{outlocation}{wildcards.sample}/rna/demuxlet/flt.bam" \
            {input.barcodes} \
            "{outlocation}{wildcards.sample}/rna/demuxlet/final_sorted.vcf" \
            "{outlocation}{wildcards.sample}/rna/demuxlet/fltbam_forpopscle.bam"

        popscle dsc-pileup --sam "{outlocation}{wildcards.sample}/rna/demuxlet/fltbam_forpopscle.bam" \
            --vcf "{outlocation}{wildcards.sample}/rna/demuxlet/final_sorted.vcf" \
            --out "{outlocation}{wildcards.sample}/rna/demuxlet/popscle_pileup"

        popscle demuxlet --plp "{outlocation}{wildcards.sample}/rna/demuxlet/popscle_pileup" \
            --vcf "{outlocation}{wildcards.sample}/rna/demuxlet/final_sorted.vcf" \
            --field GT \
            --out "{outlocation}{wildcards.sample}/rna/demuxlet/popscle_demuxlet"

        # Clean up the mess we made
        rm {outlocation}{wildcards.sample}/rna/demuxlet/bgzip_renamed_vcf.vcf.gz
        rm {outlocation}{wildcards.sample}/rna/demuxlet/vcf_renamed.vcf
        rm {outlocation}{wildcards.sample}/rna/demuxlet/bgzip_renamed_vcf.vcf.gz.csi
        rm {outlocation}{wildcards.sample}/rna/demuxlet/final_sorted.vcf
        rm {outlocation}{wildcards.sample}/rna/demuxlet/fltbam_forpopscle.bam*
        """

rule scrublet:
    input:
        barcodes = outlocation + '{sample}/rna/qc/manually_flt_data/barcodes.tsv',
        features = outlocation + '{sample}/rna/qc/manually_flt_data/genes.tsv',
        matrix = outlocation + '{sample}/rna/qc/manually_flt_data/matrix.mtx'
    output:
        singlets = outlocation + '{sample}/rna/scrublet/singlet_barcodes.txt'
    threads: 1
    conda: "scrublet"
    resources: mem='50G', time='20:00:00'
    shell:
        """
        mkdir -p {outlocation}{wildcards.sample}/rna/scrublet/

        python tools/misc/run_scrublet.py -m {input.matrix} \
          -g {input.features} \
          -b {input.barcodes} \
          -o {outlocation}{wildcards.sample}/rna/scrublet/
        """

rule flt_singletsonly_RNA:
    input:
        singlets = lambda wcs: outlocation + '{sample}/confident_multiome_singlets.txt' if PROFILE[wcs.sample] == 'Multiome' else outlocation + '{sample}/rna/confident_rna_singlets.txt',
        bam = outlocation + '{sample}/rna/qc/manually_flt_data/bam.bam',
        genes = outlocation + '{sample}/rna/qc/manually_flt_data/genes.tsv'
    output:
        outdir = directory(outlocation + '{sample}/rna/singlet_data_only/'),
        bam = outlocation + '{sample}/rna/singlet_data_only/flt_bam.bam',
        barcodes = outlocation + '{sample}/rna/singlet_data_only/barcodes.tsv',
        genes = outlocation + '{sample}/rna/singlet_data_only/genes.tsv'
    resources: mem='50G', time='5:00:00'
    shell:
        """
        mkdir -p {output.outdir}

        # Test if filtered bam exists and remove if so
        if test -f {output.bam}; then
            rm {output.bam}
        fi

        # Subset the bam file and index
        subset-bam --bam {input.bam} \
            --cell-barcodes {input.singlets} \
            --out-bam {output.bam} \
            --cores 15
        samtools index {output.bam}
        ln -s {input.singlets} {output.barcodes}
        ln -s {input.genes} {output.genes}

        """

rule flt_singletsonly_ATAC:
    input:
        # singlets = outlocation + '{sample}/atac/confident_rna_singlets.txt',
        singlets = lambda wcs: outlocation + '{sample}/confident_multiome_singlets.txt' if PROFILE[wcs.sample] == 'Multiome' else outlocation + '{sample}/atac/confident_atac_singlets.txt',
        bam = outlocation + '{sample}/atac/qc/manually_flt_data/atacbam.bam',
        genes = outlocation + '{sample}/atac/qc/manually_flt_data/features.tsv.gz'
    output:
        outdir = directory(outlocation + '{sample}/atac/singlet_data_only/'),
        bam = outlocation + '{sample}/atac/singlet_data_only/flt_bam.bam',
        barcodes = outlocation + '{sample}/atac/singlet_data_only/barcodes.tsv',
        genes = outlocation + '{sample}/atac/singlet_data_only/features.tsv.gz'
    resources: mem='50G', time='5:00:00'
    shell:
        """
        mkdir -p {output.outdir}


        # Test if filtered bam exists and remove if so
        if test -f {output.bam}; then rm {output.bam}; fi

        # Subset the bam file and index
        subset-bam --bam {input.bam} \
            --cell-barcodes {input.singlets} \
            --out-bam {output.bam} \
            --cores 15
        samtools index {output.bam}
        ln -s {input.singlets} {output.barcodes}
        ln -s {input.genes} {output.genes}

        """

use rule map_minimap2 as dmx_map_minimap2 with:
    input:
        bam = outlocation + '{sample}/rna/singlet_data_only/flt_bam.bam',
        barcodes = outlocation + '{sample}/rna/singlet_data_only/barcodes.tsv'
    output:
        maptagged = outlocation + '{sample}/rna/demultiplex/maptagged_sorted.bam'
    params:
        outdir = directory(outlocation + '{sample}/rna/demultiplex/')

use rule map_bwamem as dmx_map_bwamem with:
    input:
        bam = outlocation + '{sample}/atac/singlet_data_only/flt_bam.bam',
        singletbarcodes = outlocation + '{sample}/confident_multiome_singlets.txt'
    output:
        maptagged = outlocation + '{sample}/atac/demultiplex/maptagged_sorted.bam'
    params:
        outdir = directory(outlocation + '{sample}/atac/demultiplex/')

use rule callsnps as dmx_callsnps with:
    input:
        bam = outlocation + '{sample}/{datatype}/singlet_data_only/flt_bam.bam'
    output:
        freebayesvcf = outlocation + '{sample}/{datatype}/demultiplex/freebayes.vcf'

use rule run_souporcellmodel as dmx_run_souporcellmodel with:
    input:
        snps = rules.dmx_callsnps.output.freebayesvcf,
        bam = outlocation + '{sample}/{datatype}/demultiplex/maptagged_sorted.bam'
    output:
        clusters = outlocation + '{sample}/{datatype}/demultiplex/clusters.tsv',
        cluster_genotypes = outlocation + '{sample}/{datatype}/demultiplex/cluster_genotypes.vcf'
    params:
        barcodes = outlocation + '{sample}/{datatype}/singlet_data_only/barcodes.tsv',
        umiflag = lambda wcs: "" if wcs.datatype == 'atac' else '--umi',
        outdir = directory(outlocation + '{sample}/{datatype}/demultiplex'),
        ndonors = lambda wcs: DONORS[wcs.sample]
