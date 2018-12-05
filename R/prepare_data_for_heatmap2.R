#' Prepare data from all enrichment_tables and apply filters relying on the evaluation of all samples.
#'
#' @param  clusterTermsBy string indicating a metric to use to cluster GO terms.
prepare_data_for_heatmap2 <- function(enrichmentTables,
                                      clusterTermsBy=NULL,
                                      goLabels=c('name','ID','ID-name')){
    d <- do.call(Map, c(f=rbind, enrichmentTables))
    d <- do.call(rbind, d)
    splitted_rownames <- matrix(unlist(strsplit(x=rownames(d), split='.', fixed=TRUE)), ncol=3, byrow=TRUE)
    colnames(splitted_rownames) <- c('Ontology', 'Sample','Index')
    d <- cbind(splitted_rownames, d)
    rownames(d) <- NULL
    d$Sample <- as.character(d$Sample)

    # Add n_query_size to sample name
    for (sample in unique(d$Sample)){
        d$Sample[d$Sample == sample] <- paste0(sample, ' (', attributes(enrichmentTables[[sample]])$n_queried_regions,')')
    }
    
    d$uniqueId <- paste(d$Ontology, d$ID, d$name, sep='_')
    
    # Add pass_test info
    pass_tests_for_all_samples <- reshape2::dcast(d, uniqueId ~ Sample, value.var = "pass_tests")
    samples_cols <- 2:ncol(pass_tests_for_all_samples)
    pass_tests_for_all_samples$pass_tests_for_any_samples <- apply(pass_tests_for_all_samples[,samples_cols],
                                        MARGIN=1,
                                        FUN=any)
    pass_tests_for_all_samples[,samples_cols] <- NULL

    d <- merge(d,pass_tests_for_all_samples)
    d <- d[d$pass_tests_for_any_samples,]
    # Working but super-slow.
    #pass_tests_for_any_samples <- function(ID='GO:0005515',
    #                                       d){
    #    test <- any(d[d$ID == ID,'pass_signif_tests']) & any(d[d$ID == ID,'pass_post_filter_tests'])
    #    print(ID)
    #    return(test)
    #}
    
    # Add Z-score transformation for Binom_Fold_Enrichment
    tmp <- reshape2::dcast(d, uniqueId ~ Sample, value.var = "Binom_Fold_Enrichment")
    samples_cols <- 2:ncol(tmp)
    tmp[,samples_cols] <- t(scale(t(tmp[,samples_cols])))
    tmp <- reshape2::melt(tmp, id.vars='uniqueId', value.name='Zscore_Binom_Fold_Enrichment', variable.name='Sample')
    tmp$'Zscore_Binom_Fold_Enrichment'[is.nan(tmp$'Zscore_Binom_Fold_Enrichment')] <- 0
    d <- merge(d,tmp)


    if (goLabels == 'ID-name'){
        d$label <- paste(d$ID, d$name)
    } else {
        d$label <- d[,goLabels]
    }

    # Maybe try fastcluster here because hclust find the data too big.
    # Also clusters ontologies separately first to reduce memory footprint.
    if(!is.null(clusterTermsBy)){
        tmp <- reshape2::dcast(d, uniqueId ~ Sample, value.var = clusterTermsBy)
        
        hclust_GOterms <- hclust(dist(tmp[,2:ncol(tmp)],
                                      method = "euclidean"),
                                 method = "ward.D" )
        d$label <- factor(x=d$uniqueId,
                          levels=d$uniqueId[hclust_GOterms$order],
                          labels=d$label)
    }

    id_vars <- c("Sample","ID","name","Ontology","Index","signif_binom","signif_hyper","pass_signif_tests", "pass_post_filter_tests","pass_tests","uniqueId","label")
    
    measure_vars <- c("Zscore_Binom_Fold_Enrichment",
                      #"Binom_Genome_Fraction",
                      #"Binom_Expected",
                      "Binom_Observed_Region_Hits",
                      "Binom_Fold_Enrichment",
                      #"Binom_Region_Set_Coverage",
                      #"Binom_Raw_PValue",
                      #"Binom_Adjp_BH",
                      #"Hyper_Total_Genes",
                      #"Hyper_Expected",
                      "Hyper_Observed_Gene_Hits",
                      "Hyper_Fold_Enrichment",
                      #"Hyper_Gene_Set_Coverage",
                      #"Hyper_Term_Gene_Coverage",
                      #"Hyper_Raw_PValue",
                      #"Hyper_Adjp_BH",
                      #"Binom_Rank",
                      #"Hyper_Rank",
                      #"Binom_Bonf_PValue",
                      #"Hyper_Bonf_PValue",
                      "mlog10_Binom_Bonf_PValue",
                      "mlog10_Hyper_Bonf_PValue",
                      "mlog10_Binom_BH_PValue",
                      "mlog10_Hyper_BH_PValue",
                      #"Post_Filter_Binom_Rank",
                      #"Post_Filter_Hyper_Rank",
                      "sum_Norm_mlog10_Binom_Bonf_PValue",
                      "sum_Norm_mlog10_Hyper_Bonf_PValue",
                      "sum_Signif_Norm_mlog10_Binom_Bonf_PValue",
                      "sum_Signif_Norm_mlog10_Hyper_Bonf_PValue",
                      "sum_Displayed_Norm_mlog10_Binom_Bonf_PValue",
                      "sum_Displayed_Norm_mlog10_Hyper_Bonf_PValue")


    scaled <- apply(X=d[,measure_vars],
                 MARGIN=2,
                 FUN=scales::rescale)
    scaled <- data.frame(d[,id_vars], scaled, stringsAsFactors=FALSE)
    scaled <- reshape2::melt(scaled,
                             id.vars = id_vars,
                             measure.vars = measure_vars,
                             variable.name='metric',
                             value.name='scaled')
    #scaled$metric <- as.character(scaled$metric)
    #scaled$uniqueIdMetric <- paste0(scaled$uniqueId, scaled$metric)

    values <- reshape2::melt(d, 
                             id.vars = id_vars, 
                             measure.vars=measure_vars,
                             variable.name='metric')
    #melted$uniqueIdMetric <- paste0(melted$uniqueId, melted$metric)
    melted <- merge(values, scaled)

    # Replacing underscore for better-looking labels in heatmap
    melted$metric <- gsub(pattern='_',
                          replacement=' ',
                          x=melted$metric)

    #tmp <- merge(melted, scaled, by='uniqueIdMetric')
    #dim(tmp)

    #library(dplyr)
    #tmp <- left_join(melted, scaled, by = c('uniqueId', 'metric'))

    #melted <- merge(melted, scaled)
    return(melted)
}
