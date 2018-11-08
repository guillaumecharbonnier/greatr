#!/usr/bin/env Rscript
#------------------------------------------------------------------------------
# Guillaume Charbonnier
# Created: 2018-10-01 13:54:22
#------------------------------------------------------------------------------

"A script to produces heatmaps from multiple samples automated queries to GREAT. 

Usage: 
    GREATR.R [-i <indir>] [-f <files>] [-o <outdir>] [-a <assembly>] [-g <ontology>] [-m <filterMetrics>] [-s <filterGreaterLowerThans>] [-t <filterThresholds>] [-r <subsampleReplicates>]

Options:
    -i <indir>      Input directory. [default: .]
    -f <files>      Input files as comma-separated list. Combined with <indir>. If not provided, all bed files inside <indir> will be taken.
    -o <outdir>     Output directory. [default: .]
    -a <assembly>   Assembly (hg19,mm9,mm10) [default: hg19]
    -g <ontology>   Ontology. Possible values and aliases (GOBP: GO Biological Process, MSDBP: MSigDB Pathway), [default: MSigDB Pathway]
    -r <subsampleReplicates>    Number of replicates for subsampling, no subsampling if unspecified.
    -m <filterMetrics>  Comma separated list of metrics to use for filtering. Defaults try to match GREAT default settings. See together <filterGreaterLowerThans> and <filterThresholds>. [default: Binom_Fold_Enrichment,Binom_Adjp_BH,Hyper_Adjp_BH,Post_Filter_Binom_Rank]
    -s <filterGreaterLowerThans>    Vector of 'lower' or 'greater' to apply fiterThresholds values to filterMetrics. [default: greater,lower,lower,lower]
    -t <filterThresholds>   Vector of thresholds to apply to filterMetrics. Values matching treshold are kept, i.e. '<=' and '>=' are used for comparison. [default: 2,0.05,0.05,10]
    " -> doc


# -----------------------------------------------------------------------------
# Arg parsing
# -----------------------------------------------------------------------------
library(docopt)
opts <- docopt(doc)

opts$m <- unlist(strsplit(x=opts$m, split=','))
opts$s <- unlist(strsplit(x=opts$s, split=','))
opts$t <- as.numeric(unlist(strsplit(x=opts$t, split=',')))
print(opts)
## -----------------------------------------------------------------------------
## Function calling
## -----------------------------------------------------------------------------

library(greatr)
#source('opt/greatr/R/great_heatmap.R')
plot_facet_heatmaps(indir=opts$i,
#make_preset_heatmaps(indir=opts$i,
                     files=opts$f,
                     outdir=opts$o,
                     assembly=opts$a,
                     filterMetrics=opts$m,
                     filterGreaterLowerThans=opts$s,
                     filterThresholds=opts$t,
                     subsampleReplicates=opts$r,
                     ontology=opts$g)


