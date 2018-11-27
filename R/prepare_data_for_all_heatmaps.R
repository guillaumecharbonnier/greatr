prepare_data_for_all_heatmaps <- function(enrichmentTables=enrichmentTables,
                             showedMetric=showedMetric,
                             transformation=transformation,
                             filterMetrics=filterMetrics,
                             filterGreaterLowerThans=filterGreaterLowerThans,
                             filterThresholds=filterThresholds,
                             ontology=ontology,
                             goLabels=goLabels){
    for (ontology in names(enrichment_tables[[1]])){
        data_for_all_heatmaps <- prepare_data_for_heatmap(enrichmentTables=enrichmentTables,
                                                          showedMetric='Binom_Fold_Enrichment',
                                                          transformation='',
                                                          filterMetrics=filterMetrics,
                                                          filterGreaterLowerThans=filterGreaterLowerThans,
                                                          filterThresholds=filterThresholds)
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
prepare_data_for_heatmap <- function(enrichmentTables,
                                     showedMetric='Binom_Bonf_PValue',
                                     transformation=c('mlog10','Zscore','x100','unsignifAsNa'),
                                     filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                     filterThresholds=c(1.5,0.05,0.05,5),
                                     filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                     ontology='MSigDB Pathway',
                                     orderGOTerms=FALSE,
                                     goLabels=c('name','ID','ID-name')){
    filterThresholds <- as.numeric(filterThresholds)
    goLabels <- match.arg(goLabels) # Default to the first one.
    data_for_heatmap <- enrichmentTables[[1]][[ontology]]
    browser()
    print(paste("debug data_for_heatmap", str(data_for_heatmap)))
    col_to_remove_after_merge <- colnames(data_for_heatmap)

    #col_to_remove_after_merge <- col_to_remove_after_merge[! col_to_remove_after_merge %in% c('ID','name')]
    for (sample_label in names(enrichmentTables)){
        to_merge <- enrichmentTables[[sample_label]][[ontology]]
        if ('unsignifAsNa' %in% transformation){
            to_merge$unsignifAsNaShowedMetric <- to_merge[,showedMetric]
            to_merge$unsignifAsNaShowedMetric[!to_merge$pass_signif_tests] <- NA
        }
        data_for_heatmap <- merge(data_for_heatmap,
                                  to_merge,
                                  by=c('ID','name'),
                                  suffix=c('',
                                           paste0('.',
                                                  sample_label)))
    }
    if (goLabels == 'ID-name'){
        labels_for_heatmap <- paste(data_for_heatmap[,'ID'],
                                    data_for_heatmap[,'name'],
                                    sep=' - ')
    } else {
        labels_for_heatmap <- data_for_heatmap[,goLabels]
    }
    data_for_heatmap[,col_to_remove_after_merge] <- NULL
    # Deduplicating labels with suffix number.
    # It occurs for some databases when using only "name".
    labels_for_heatmap <- ifelse(duplicated(labels_for_heatmap) | duplicated(labels_for_heatmap, fromLast=TRUE), 
                                 paste(labels_for_heatmap,
                                       ave(labels_for_heatmap, 
                                           labels_for_heatmap, 
                                           FUN=seq_along),
                                       sep='_'), 
                                 labels_for_heatmap)
    rownames(data_for_heatmap) <- labels_for_heatmap
    print(paste('Number of GO terms before filtering: ',nrow(data_for_heatmap)))
    data_for_heatmap$any_sample_pass_signif_tests <- apply(data_for_heatmap[,grep(pattern="pass_signif_tests",
                                                                                  x=colnames(data_for_heatmap))],
                                                           1,
                                                           any)
    data_for_heatmap$any_sample_pass_post_filter_tests <- apply(data_for_heatmap[,grep(pattern="pass_post_filter_tests",
                                                                                       x=colnames(data_for_heatmap))],
                                                                1,
                                                                any)
    data_for_heatmap <- data_for_heatmap[data_for_heatmap$any_sample_pass_signif_tests & data_for_heatmap$any_sample_pass_post_filter_tests,]
    #for (filterNumber in 1:length(filterMetrics)){
    #    pattern <- paste0('^',
    #                      filterMetrics[filterNumber])
    #    threshold <- filterThresholds[filterNumber]
    #    print(paste('filterMetric pattern:',
    #                  pattern))
    #    print(paste('threshold:',
    #                threshold))

    #    if (filterGreaterLowerThans[filterNumber] == 'lower'){
    #        minFilterMetric <- apply(X=data_for_heatmap[,grep(x=colnames(data_for_heatmap),
    #                                                          pattern=pattern)],
    #                                 MARGIN=1,
    #                                 FUN=min)
    #        data_for_heatmap <- data_for_heatmap[minFilterMetric <= threshold,]
    #        #if (pattern == '^Post_Filter_Binom_Rank'){
    #        #    browser()
    #        #}
    #    } else if (filterGreaterLowerThans[filterNumber] == 'greater'){
    #        maxFilterMetric <- apply(X=data_for_heatmap[,grep(x=colnames(data_for_heatmap),
    #                                                          pattern=pattern)],
    #                                 MARGIN=1,
    #                                 FUN=max)
    #        data_for_heatmap <- data_for_heatmap[maxFilterMetric >= threshold,]
    #    }
    #    print(paste('Number of GO terms after filtering for',
    #                filterMetrics[filterNumber],
    #                filterGreaterLowerThans[filterNumber],
    #                'than',
    #                filterThresholds[filterNumber],
    #                ':',
    #                nrow(data_for_heatmap)))

    #    #print(apply(data_for_heatmap[,grep(pattern='Post_Filter_Binom_Rank', x=colnames(data_for_heatmap))],1,min))
    #    # Computing post filter ranks to be able to select only top N terms if Post_Filter_(Binom|Hyper)_Rank is used in next iteration.
    #    post_filter_ranks <- apply(data_for_heatmap[,grep(pattern='^(Binom|Hyper)_Raw_PValue',
    #                                                           x=colnames(data_for_heatmap))], 
    #                                    MARGIN=2, 
    #                                    FUN=rank)
    #    colnames(post_filter_ranks) <- gsub(pattern='^(Binom|Hyper)_Raw_PValue',
    #                                             replacement="Post_Filter_\\1_Rank",
    #                                             x=colnames(post_filter_ranks))
    #    data_for_heatmap[,grep(pattern='Post_Filter_(Binom|Hyper)_Rank',
    #                           x=colnames(data_for_heatmap))] <- NULL
    #    data_for_heatmap <- cbind(data_for_heatmap,
    #                              post_filter_ranks)
    #}
    showedMetricPattern <- paste0('^',
                                  showedMetric)
    if ('unsignifAsNa' %in% transformation){
        showedMetricPattern <- paste0('^',
                                      'unsignifAsNaShowedMetric')
    }
    print(paste("debug data_for_heatmap",str(data_for_heatmap)))
    signif_binom <- data_for_heatmap[,grep(x=colnames(data_for_heatmap),
                                           pattern="^signif_binom")]
    colnames(signif_binom) <- gsub(pattern=paste0("^signif_binom",'.'),
                                  replacement='',
                                  x=colnames(signif_binom))
    signif_hyper <- data_for_heatmap[,grep(x=colnames(data_for_heatmap),
                                           pattern="^signif_hyper")]
    colnames(signif_hyper) <- gsub(pattern=paste0("^signif_hyper",'.'),
                                  replacement='',
                                  x=colnames(signif_hyper))
    data_for_heatmap <- data_for_heatmap[,grep(x=colnames(data_for_heatmap),
                                               pattern=showedMetricPattern)]
    colnames(data_for_heatmap) <- gsub(pattern=paste0(showedMetricPattern,'.'),
                                       replacement='',
                                       x=colnames(data_for_heatmap))

    data_for_heatmap <- merge_replicates(data_for_heatmap)

    if ('mlog10' %in% transformation){
        data_for_heatmap <- -log10(data_for_heatmap)
        data_for_heatmap[data_for_heatmap == Inf] <- 333
    }
    if ('pseudolog2' %in% transformation){
        data_for_heatmap <- log2(data_for_heatmap + 1)  
        #data_for_heatmap[data_for_heatmap == Inf] <- 333
    }
    if ('Zscore' %in% transformation){
        data_for_heatmap <- t(scale(t(data_for_heatmap)))
        attr(data_for_heatmap,"scaled:center") <- NULL
        attr(data_for_heatmap,"scaled:scale") <- NULL 
        data_for_heatmap[is.nan(data_for_heatmap)] <- 0
        data_for_heatmap <- as.data.frame(data_for_heatmap)
    }
    if ('x100' %in% transformation){
        data_for_heatmap <- 100 * data_for_heatmap
    }
    labels_for_heatmap <- rownames(data_for_heatmap)
    if(orderGOTerms){
        hclust_GOterms <- hclust(dist(data_for_heatmap,
                                      method = "euclidean"),
                                 method = "ward.D" )
        labels_for_heatmap <- factor(x=labels_for_heatmap,
                                     levels=labels_for_heatmap[hclust_GOterms$order])
    }
    metric <- cbind(label=labels_for_heatmap,
                                  data_for_heatmap)
    signif_binom <- cbind(label=labels_for_heatmap,
                          signif_binom)
    signif_hyper <- cbind(label=labels_for_heatmap,
                          signif_hyper)
    data_for_heatmap <- list(metric=metric,
                             signif_binom=signif_binom,
                             signif_hyper=signif_hyper)
    return(data_for_heatmap)
}
