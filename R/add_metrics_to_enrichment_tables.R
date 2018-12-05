#' Note: Post filter binom rank is somewhat deprecated as a filtering option as adjusting semantic distance should be better in any cases.
add_metrics_to_enrichment_tables <- function(enrichment_tables,
                                             filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                             filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                             filterThresholds=c('2','0.05','0.05','30')){
    for (sample_label in names(enrichment_tables)){
        for (ontology in names(enrichment_tables[[sample_label]])){
            enrichment_tables[[sample_label]][[ontology]] <- compute_additional_metrics(enrichment_tables[[sample_label]][[ontology]],
                                                                                        filterMetrics=filterMetrics,
                                                                                        filterGreaterLowerThans=filterGreaterLowerThans,
                                                                                        filterThresholds=filterThresholds)
        }
    }

    attr(enrichment_tables,'filterMetrics') <- filterMetrics
    attr(enrichment_tables,'filterGreaterLowerThans') <- filterGreaterLowerThans
    attr(enrichment_tables,'filterThresholds') <- filterThresholds

    return(enrichment_tables)
}


