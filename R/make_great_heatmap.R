library(tools) # file_path_sans_ext
library(rtracklayer) # Loading BED
library(rGREAT) # Queries
library(ggplot2)
library(ggpubr)
library(reshape) # melt
library(scales)
##library(ggdendro)
##library(mvtnorm) 
#library(gplots) # Need heatmap.2
#library(RColorBrewer) # Need a better palette for the heatmap
##library(klaR) # kmode
#library(cluster) # PAM gower

#' Produce heatmaps of ontology enrichment for a list of bed files, or a list of bed objects
#'
#' @param outdir The output directory for plots and data.
#' @param outfile The filename of the heatmap. If NULL, a smart name is given based on chosen settings.
#' @param indir The input directory for bed files, will be combined with ``files'' parameter to get full path.
#' @param files A comma-separated list of filenames
#' @param beds A GRangesList object of the bed files to query to GREAT. If provided, ``indir''  and ``files'' are ignored.
#' @param enrichmentTables A list of queried tables. If provided, ``indir'', ``files'' and ``beds'' are ignored.
#' @param sampleLabels A comma-separated list of labels to use for filenames
#' @param smartLabels If TRUE, labels will be taken from filename without extension and dirname.
#' @param showedMetric Define the metric displayed in the heatmap. Can be modified by transformation parameter.
#' @param subsampleSize Specify the number of features to draw from beds. Keep it to NULL and the function will use the number of features of the sample with lowest number of features.
#' @param subsample_replicates The number of drawing replicates for each sample.
#' @param transformation A vector of transformations that can be done on the showed metric. Note that if "mlog10" and "Zscore" are given, then log transformation is always applied before Zscore.
#' @param filterMetric The metric used to limit the number of rows displayed in the heatmap.
#' @param filterThreshold GO terms where at least one sample has its "filterMetric" value below "filterThreshold" will be taken for display in the heatmap.
#' @param assembly The genome assembly of the provided input files.
#' @param ontology The displayed ontology. 'MSigDB Pathway' by default. Choose only one from available choices. Can accept ontology with '_' in name.
#' @param goLabels Use 'ID', 'name' or the combination of both 'ID-name' from ontologies as labels on heatmaps.
#' @return enrichmentTables that can be used to shortcut bed loading and GREAT querying for subsequent use of ``make_great_heatmap'' for different metric and transformation.
#'
#' @export
make_great_heatmap <- function(outdir='.',
                               outfile=NULL,
                               indir='.',
                               files=NULL,
                               beds=NULL,
                               enrichmentTables=NULL,
                               sampleLabels,
                               smartLabels=TRUE,
                               subsampleReplicates=NULL,
                               subsampleSize=NULL,
                               showedMetric=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Binom_Adjp_BH','Binom_Raw_PValue','Hyper_Bonf_PValue','Hyper_Adjp_BH','Hyper_Fold_Enrichment','Hyper_Raw_PValue','Binom_Rank','Hyper_Rank', 'mlog10_Binom_Bonf_PValue', 'mlog10_Hyper_Bonf_PValue','sum_Norm_mlog10_Binom_Bonf_PValue','sum_Norm_mlog10_Hyper_Bonf_PValue','Post_Filter_Binom_Rank','Post_Filter_Hyper_Rank','sum_Signif_Norm_mlog10_Binom_Bonf_PValue','sum_Signif_Norm_mlog10_Hyper_Bonf_PValue','sum_Displayed_Norm_mlog10_Binom_Bonf_PValue','sum_Displayed_Norm_mlog10_Hyper_Bonf_PValue'),
                               transformation=c('mlog10','Zscore','x100','unsignifAsNa'),
                               filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                               filterGreaterLowerThans=c('greater','lower','lower','lower'),
                               filterThresholds=c('1.5','0.05','0.05','10'),
                               assembly=c('hg19','mm9','mm10'),
                               ontology=c('MSigDB Pathway','GO Molecular Function','GO Biological Process','PANTHER Pathway','Disease Ontology',
                                          'MSigDB_Pathway','GO_Molecular_Function','GO_Biological_Process','PANTHER_Pathway','Disease_Ontology'),
                               goLabels=c('name','ID-name','ID')){
    showedMetric <- match.arg(showedMetric) # Default to the first one.
    #filterMetric <- match.arg(filterMetric)
    assembly <- match.arg(assembly)
    ontology <- match.arg(ontology)
    goLabels <- match.arg(goLabels)

    ontology <- gsub(pattern='_',
                     replacement=' ',
                     x=ontology)
    if (is.null(enrichmentTables)){
        if (is.null(beds)){
            beds <- load_beds(outdir,
                              indir,
                              files,
                              assembly='hg19')
            if (!is.null(subsampleReplicates)){
                beds <- subsample_beds(beds,
                                       subsampleSize=subsampleSize,
                                       subsampleReplicates=subsampleReplicates)
            }
        }
        enrichmentTables <- query_great(outdir,
                                         assembly=assembly,
                                         beds=beds)
        enrichmentTables <- add_slim_ontologies(enrichmentTables,
                                                slimList)
        enrichmentTables <- add_metrics_to_enrichment_tables(enrichmentTables,
                                                             filterMetrics=filterMetrics,
                                                             filterGreaterLowerThans=filterGreaterLowerThans,
                                                             filterThresholds=filterThresholds)
    }
    data_for_heatmap <- prepare_data_for_heatmap(enrichmentTables=enrichmentTables,
                                                 showedMetric=showedMetric,
                                                 transformation=transformation,
                                                 filterMetrics=filterMetrics,
                                                 filterGreaterLowerThans=filterGreaterLowerThans,
                                                 filterThresholds=filterThresholds,
                                                 ontology=ontology,
                                                 goLabels=goLabels)
    p <- plot_heatmap(data_for_heatmap,
                 outdir=outdir,
                 showedMetric=showedMetric,
                 transformation=transformation)
    great_heatmap <- list(enrichmentTables=enrichmentTables,
                          data_for_heatmap=data_for_heatmap,
                          p=p)
    return(great_heatmap)
}


