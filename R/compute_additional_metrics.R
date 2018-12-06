#' Add various additional metrics to the provided enrichment table.
compute_additional_metrics <- function(enrichment_table,
                                       filterMetrics=c('Binom_Fold_Enrichment','Binom_Adjp_BH','Hyper_Adjp_BH','Post_Filter_Binom_Fold_Enrichment_Rank'),
                                       filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                       filterThresholds=c('2','0.05','0.05','10')){
    enrichment_table$Binom_Rank <- rank(enrichment_table$Binom_Raw_PValue)
    enrichment_table$Hyper_Rank <- rank(enrichment_table$Hyper_Raw_PValue)
    enrichment_table$Binom_Bonf_PValue <- p.adjust(enrichment_table$Binom_Raw_PValue, method='bonferroni')
    enrichment_table$Hyper_Bonf_PValue <- p.adjust(enrichment_table$Hyper_Raw_PValue, method='bonferroni')
    # I think before v1.12, rGREAT API had by default the adjusted p-value as these colnames.
    # It is computed here in case recent version is used.
    enrichment_table$Binom_Adjp_BH <- p.adjust(enrichment_table$Binom_Raw_PValue, method='BH')
    enrichment_table$Hyper_Adjp_BH <- p.adjust(enrichment_table$Hyper_Raw_PValue, method='BH')
    enrichment_table$mlog10_Binom_Bonf_PValue <- -log10(enrichment_table$Binom_Bonf_PValue)
    enrichment_table$mlog10_Hyper_Bonf_PValue <- -log10(enrichment_table$Hyper_Bonf_PValue)
    enrichment_table$mlog10_Binom_Bonf_PValue[enrichment_table$mlog10_Binom_Bonf_PValue == Inf] <- 333 
    enrichment_table$mlog10_Hyper_Bonf_PValue[enrichment_table$mlog10_Hyper_Bonf_PValue == Inf] <- 333 
    enrichment_table$mlog10_Binom_BH_PValue <- -log10(enrichment_table$Binom_Adjp_BH)
    enrichment_table$mlog10_Hyper_BH_PValue <- -log10(enrichment_table$Hyper_Adjp_BH)
    enrichment_table$mlog10_Binom_BH_PValue[enrichment_table$mlog10_Binom_BH_PValue == Inf] <- 333 
    enrichment_table$mlog10_Hyper_BH_PValue[enrichment_table$mlog10_Hyper_BH_PValue == Inf] <- 333 

    filterPostFilter <- grepl(pattern='^Post_Filter', x=filterMetrics)
    filterPreFilter <- !grepl(pattern='^Post_Filter', x=filterMetrics)
    filterPreFilterMetrics <- filterMetrics[filterPreFilter]
    filterPreFilterThresholds <- filterThresholds[filterPreFilter]
    filterPreFilterGreaterLowerThans <- filterGreaterLowerThans[filterPreFilter]
    filterPostFilterMetrics <- filterMetrics[filterPostFilter]
    filterPostFilterThresholds <- filterThresholds[filterPostFilter]
    filterPostFilterGreaterLowerThans <- filterGreaterLowerThans[filterPostFilter]

    # These two tests should be used to set font to bold and or italic in plots:
    enrichment_table$signif_binom <- pass_tests(enrichmentTable=enrichment_table,
                                                filterMetrics='Binom_Adjp_BH',
                                                filterThresholds='0.05',
                                                filterGreaterLowerThans='lower')
    enrichment_table$signif_hyper <- pass_tests(enrichmentTable=enrichment_table,
                                                filterMetrics='Hyper_Adjp_BH',
                                                filterThresholds='0.05',
                                                filterGreaterLowerThans='lower')

    enrichment_table$pass_signif_tests <- pass_tests(enrichmentTable=enrichment_table,
                                                            filterMetrics=filterPreFilterMetrics,
                                                            filterThresholds=filterPreFilterThresholds,
                                                            filterGreaterLowerThans=filterPreFilterGreaterLowerThans)

    #enrichment_table$pass_similarity_test <- pass_similarity_test(enrichment_table, ...)

    enrichment_table[enrichment_table$pass_signif_tests, 'Post_Filter_Binom_PValue_Rank'] <- rank(enrichment_table[enrichment_table$pass_signif_tests,'Binom_Raw_PValue'])
    enrichment_table[enrichment_table$pass_signif_tests, 'Post_Filter_Hyper_PValue_Rank'] <- rank(enrichment_table[enrichment_table$pass_signif_tests,'Hyper_Raw_PValue'])
    enrichment_table[enrichment_table$pass_signif_tests, 'Post_Filter_Binom_Fold_Enrichment_Rank'] <- rank( - enrichment_table[enrichment_table$pass_signif_tests,'Binom_Fold_Enrichment'])
    enrichment_table[enrichment_table$pass_signif_tests, 'Post_Filter_Hyper_Fold_Enrichment_Rank'] <- rank( - enrichment_table[enrichment_table$pass_signif_tests,'Hyper_Fold_Enrichment'])

#    enrichment_table$Post_Filter_Binom_Rank <- ifelse(enrichment_table$pass_signif_tests,
#                                                      rank(enrichment_table[enrichment_table$pass_signif_tests,'Binom_Raw_PValue']),
#                                                      NA)
#    enrichment_table$Post_Filter_Hyper_Rank <- ifelse(enrichment_table$pass_signif_tests,
#                                                      rank(enrichment_table[enrichment_table$pass_signif_tests,'Hyper_Raw_PValue']),
#                                                      NA)
    enrichment_table$sum_Norm_mlog10_Binom_Bonf_PValue <- enrichment_table$mlog10_Binom_Bonf_PValue / sum(enrichment_table$mlog10_Binom_Bonf_PValue, na.rm=T)
    enrichment_table$sum_Norm_mlog10_Hyper_Bonf_PValue <- enrichment_table$mlog10_Hyper_Bonf_PValue / sum(enrichment_table$mlog10_Hyper_Bonf_PValue, na.rm=T)

    enrichment_table$sum_Signif_Norm_mlog10_Binom_Bonf_PValue <- ifelse(enrichment_table$pass_signif_tests,
                                                                        enrichment_table$mlog10_Binom_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_signif_tests,'mlog10_Binom_Bonf_PValue']),
                                                                        NA)
    enrichment_table$sum_Signif_Norm_mlog10_Hyper_Bonf_PValue <- ifelse(enrichment_table$pass_signif_tests,
                                                                        enrichment_table$mlog10_Hyper_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_signif_tests,'mlog10_Hyper_Bonf_PValue']),
                                                                        NA)
    # Set post_filter as true in case there is none of them.
    #enrichment_table$pass_post_filter_tests <- TRUE
    if (length(filterPostFilterMetrics) > 0){
        enrichment_table$pass_post_filter_tests <- pass_tests(enrichmentTable=enrichment_table,
                                                              filterMetrics=filterPostFilterMetrics,
                                                              filterThresholds=filterPostFilterThresholds,
                                                              filterGreaterLowerThans=filterPostFilterGreaterLowerThans)
    } else {
        enrichment_table$pass_post_filter_tests <- TRUE
    }
    enrichment_table$pass_tests <- enrichment_table$pass_signif_tests & enrichment_table$pass_post_filter_tests

    enrichment_table$sum_Displayed_Norm_mlog10_Binom_Bonf_PValue <- ifelse(enrichment_table$pass_post_filter_tests,
                                                                           enrichment_table$mlog10_Binom_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_post_filter_tests,'mlog10_Binom_Bonf_PValue']),
                                                                           NA)
    enrichment_table$sum_Displayed_Norm_mlog10_Hyper_Bonf_PValue <- ifelse(enrichment_table$pass_post_filter_tests,
                                                                           enrichment_table$mlog10_Hyper_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_post_filter_tests,'mlog10_Hyper_Bonf_PValue']),
                                                                           NA)


    return(enrichment_table)
}

