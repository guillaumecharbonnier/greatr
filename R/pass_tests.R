#' Look for the filterMetrics and return a boolean vector flagging GO failing to pass significance test.
#' @return Could return any metric column from enrichmentTable. I should rewrite this function to work directly on the vector. Later...
pass_tests <- function(enrichmentTable,
                       filterMetrics=c('Binom_Fold_Enrichment','Binom_Adjp_BH','Hyper_Adjp_BH'),
                       filterThresholds=c(2,0.05,0.05),
                       filterGreaterLowerThans=c('greater','lower','lower')){
    filterThresholds <- as.numeric(filterThresholds)
    filterMetricsDf <- enrichmentTable[, filterMetrics, drop=FALSE]
    filterMetricsBool <- filterMetricsDf
    for (filterNumber in 1:length(filterMetrics)){
        threshold <- filterThresholds[filterNumber]
        if (filterGreaterLowerThans[filterNumber] == 'lower'){
            filterMetricsBool[,filterNumber] <- filterMetricsDf[,filterNumber] <= threshold
        } else if (filterGreaterLowerThans[filterNumber] == 'greater'){
            filterMetricsBool[,filterNumber] <- filterMetricsDf[,filterNumber] >= threshold
        }
    }
    d <- apply(filterMetricsBool,1,all)
    d[is.na(d)] <- FALSE
    return(d)
}

