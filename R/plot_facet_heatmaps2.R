plot_facet_heatmaps2 <- function(d,
                                 outdir='.',
                                 files=NULL,
                                 slimList=NULL){
    
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
    melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted$metric <- 'Binom Fold\nEnrichment'
    melted$scaled <- scales::rescale(melted$value)

    #After the first heatmap we can get an idea of the dimension of the final plot
    nr <- nrow(great_heatmap$data_for_heatmap)
    nc <- ncol(great_heatmap$data_for_heatmap)
    height_heatmap <- 0.2 * nr
    height_top <- 0.5 # I may improve this looking at how much '\n' I have in 'metric' and multiply this number by 0.2
    height_bottom <- 0.2 + 0.15 * max(nchar(colnames(great_heatmap$data_for_heatmap)))
    height <- height_top + height_bottom + height_heatmap
    width_heatmap <- nc*0.25
    width_labels <- 0.15* max(nchar(as.character(great_heatmap$data_for_heatmap$label)))
    # Manually set n_heatmaps here but can be improved later
    n_heatmaps <- 7
    width <- n_heatmaps * width_heatmap + width_labels

    #great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
    #                   outdir=outdir,
    #                   ontology=ontology,
    #                   filterMetrics=filterMetrics,
    #                   filterGreaterLowerThans=filterGreaterLowerThans,
    #                   filterThresholds=filterThresholds,
    #                   showedMetric='Binom_Fold_Enrichment',
    #                   transformation='unsignifAsNa')
    #plots[['Binom_Fold_Enrichment']] <- plots[['Binom_Fold_Enrichment']] + theme(axis.text.y = element_blank())

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Fold_Enrichment',
                       transformation='Zscore')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'Zscore by GO term\nBinom Fold Enrich'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Rank',
                       transformation='')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'Binom Rank'
    append_to_melted$scaled <- scales::rescale(-append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Post_Filter_Binom_Rank',
                       transformation='')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'Post-Filter\nBinom Rank'
    append_to_melted$scaled <- scales::rescale(-append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Adjp_BH',
                       transformation='mlog10')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'mlog10 Binom\nAdjp_BH'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Hyper_Adjp_BH',
                       transformation='mlog10')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'mlog10 Hyper\nAdjp_BH'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- '% of sum of mlog10\nBonf PValue in sample'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Signif_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- '% of sum of signif mlog10\nBonf PValue in sample'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                                        outdir=outdir, 
                                        ontology=ontology,
                                        filterMetrics=filterMetrics,
                                        filterGreaterLowerThans=filterGreaterLowerThans,
                                        filterThresholds=filterThresholds,
                                        showedMetric='sum_Displayed_Norm_mlog10_Binom_Bonf_PValue',
                                        transformation='x100')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- '% of sum of displayed mlog10\nBonf PValue in sample'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    n_heatmaps <- 4
    print(unique(melted$metric))
    melted_filters <- melted[melted$metric %in% c('Post-Filter\nBinom Rank',
                                                  'mlog10 Binom\nAdjp_BH',
                                                  'mlog10 Hyper\nAdjp_BH',
                                                  'Binom Fold\nEnrichment'),]
    outpath <- file.path(outdir,'facet_filters.pdf')
    plot_melted_data(melted_filters,
                     outpath)

    melted_filters <- melted[melted$metric %in% c('Zscore by GO term\nBinom Fold Enrich',
                                                  'Binom Fold\nEnrichment'),]
    outpath <- file.path(outdir,'facet_fold_enrich.pdf')
    plot_melted_data(melted_filters,
                     outpath)

    melted_filters <- melted[melted$metric %in% c('% of sum of displayed mlog10\nBonf PValue in sample',
                                                  '% of sum of signif mlog10\nBonf PValue in sample',
                                                  '% of sum of mlog10\nBonf PValue in sample',
                                                  'Binom Fold\nEnrichment'),]
    outpath <- file.path(outdir,'facet_sum_norm_mlog10.pdf')
    plot_melted_data(melted_filters,
                     outpath)

    outpath <- file.path(outdir,'facet_all.pdf')
    plot_melted_data(melted,
                     outpath)
}


