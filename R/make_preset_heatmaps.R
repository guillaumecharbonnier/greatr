#' Make multiple heatmaps with standard metrics and their transformation.
#'
#' @export
make_preset_heatmaps <- function(indir='.',
                                 outdir='.',
                                 files=NULL,
                                 filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                 filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                 filterThresholds=c('1.5','0.05','0.05','10'),
                                 subsampleReplicates=NULL,
                                 assembly=c('hg19','mm9','mm10'),
                                 ontology='GO Biological Process'){
    if (ontology=='GOBP'){
        ontology <- 'GO Biological Process'
    }
    if (ontology=='MSDBP'){
        ontology <- 'MSigDB Pathway'
    }
    great_heatmap <- make_great_heatmap(indir=indir,
                                           outdir=outdir,
                                           ontology=ontology,
                                           files=files,
                                           filterMetrics=filterMetrics,
                                           filterGreaterLowerThans=filterGreaterLowerThans,
                                           filterThresholds=filterThresholds,
                                           subsampleReplicates=subsampleReplicates,
                                           showedMetric='Binom_Fold_Enrichment',
                                           transformation='')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Fold_Enrichment',
                       transformation='unsignifAsNa')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Post_Filter_Binom_Rank',
                       transformation='')
    ## Post_Filter_Binom_Rank is computed after NA transformation.
    ## Thus it can not work in the current implementation.
    ## It is not very useful so I am not going to make effort to make it works.
    #    enrichmentTables <- make_great_heatmap(enrichmentTables=enrichmentTables, 
    #                       outdir=outdir, 
    #                       ontology=ontology,
    #                       filterMetrics=filterMetrics,
    #                       filterGreaterLowerThans=filterGreaterLowerThans,
    #                       filterThresholds=filterThresholds,
    #                       showedMetric='Post_Filter_Binom_Rank',
    #                       transformation='unsignifAsNa')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Bonf_PValue',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Hyper_Bonf_PValue',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Adjp_BH',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Hyper_Adjp_BH',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Bonf_PValue',
                       transformation=c('mlog10','unsignifAsNa'))
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Bonf_PValue',
                       transformation=c('Zscore','mlog10'))
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Fold_Enrichment',
                       transformation='Zscore')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue',
                       transformation=c('x100','unsignifAsNa'))
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Signif_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Signif_Norm_mlog10_Binom_Bonf_PValue',
                       transformation=c('x100','unsignifAsNa'))
}

