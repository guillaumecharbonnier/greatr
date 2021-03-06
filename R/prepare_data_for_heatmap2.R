#' Prepare data from all enrichment_tables and apply filters relying on the evaluation of all samples.
#'
#' @param enrichmentTables list of rGREAT outputs modified by add_additional_metrics function.
#' @param outdir path for the tabulatated file containing result for all samples.
#' @param clusterTermsBy string indicating a metric to use to cluster GO terms. If 'ontology_order' is given, terms will be displayed in their order in enrichmentTables.
#' @param goLabels Specify if name (ex: immune response), ID (ex: GO:0006955) or both together are used for labels in heatmaps.
#' @return melted dataframe to use with ggplot2
#' @export
prepare_data_for_heatmap2 <- function(enrichmentTables,
                                      outdir='.',
                                      clusterTermsBy='ontology_order',
                                      goLabels=c('name','ID','ID-name')){
    # Here I should add Sample and Ontology as columns of each table
    for (sample in names(enrichmentTables)){
        for (ontology in names(enrichmentTables[[sample]])){
            et <- enrichmentTables[[sample]][[ontology]]
            et$Sample <- sample
            et$Ontology <- ontology
            et$uniqueId <- paste(et$Ontology, et$ID, et$name, sep='_')
            # setting default order
            #et$uniqueId <- factor(x=et$uniqueId, levels=et$uniqueId[1:nrow(et)], labels=et$uniqueId)
            # reverse order seems to be better as ids will be displayed in the order of the ontology:
            et$uniqueId <- factor(x=et$uniqueId, levels=rev(et$uniqueId[1:nrow(et)]))

            if (goLabels == 'ID-name'){
                et$label <- paste(et$ID, et$name)
            } else {
                et$label <- et[,goLabels]
            }

            # Using factor levels is not enough to ensure row orders in heatmap because some label are shared by multiple ontologies where their order is different.
            # Thus, label_order store the order of the label for a given ontology and is used during plotting.
            # See: https://stackoverflow.com/questions/42112100/ggplot-facet-wrap-with-specific-order-of-variables-in-each-facet
            #et$label_order <- 1:nrow(et)
            #et$label <- factor(x=et$uniqueId,
            #                   levels=et$uniqueId[1:nrow(et)],
            #                   labels=et$label)
            enrichmentTables[[sample]][[ontology]] <- et
            ### NOTE TODO:
            ## EVEN WITH THAT IDs ARE displayed alphabetically and not with factor order,
            # go back here and debug.
            # 

        }
    }
    # check if this call ruin label factor levels:
    d <- do.call(Map, c(f=rbind, enrichmentTables))
    d <- do.call(rbind, d)
    # Command below will fail if '.' is containained in ontology name (which is the case for similarity-filtered ones).
    # That
    #splitted_rownames <- matrix(unlist(strsplit(x=rownames(d), split='.', fixed=TRUE)), ncol=3, byrow=TRUE)
    #colnames(splitted_rownames) <- c('Ontology', 'Sample','Index')
    #d <- cbind(splitted_rownames, d)
    rownames(d) <- NULL
    #d$Sample <- as.character(d$Sample)
    # Add n_query_size to sample name
    for (sample in unique(d$Sample)){
        d$Sample[d$Sample == sample] <- paste0(sample, ' (', attributes(enrichmentTables[[sample]])$n_queried_regions,')')
    }
    # We want to preserve the order of samples provided as input instead of alphabetical order
    d$Sample <- factor(x=d$Sample, levels=unique(d$Sample))
    
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
    
    # Maybe try fastcluster here because hclust find the data too big.
    # Also clusters ontologies separately first to reduce memory footprint.
    if (clusterTermsBy != 'ontology_order'){
        tmp <- reshape2::dcast(d, uniqueId ~ Sample, value.var = clusterTermsBy)
        
        hclust_GOterms <- hclust(dist(tmp[,2:ncol(tmp)],
                                      method = "euclidean"),
                                 method = "ward.D" )
        d$uniqueId <- factor(x=d$uniqueId,
                          levels=tmp$uniqueId[hclust_GOterms$order],
                          labels=tmp$uniqueId[hclust_GOterms$order])
    }

    # Write tables
    #outdir_tables <- file.path(outdir, 'tables')
    #dir.create(outdir_tables, recursive=TRUE)
    write.table(x=d,
                file=file.path(outdir,'all.tsv'), 
                sep='\t', 
                row.names=FALSE)
    for (metric in c('mlog10_Hyper_BH_PValue', 'Binom_Fold_Enrichment')){
        for (ontology in unique(d$Ontology)){
            dr <- reshape2::dcast(d[d$Ontology == ontology,], uniqueId ~ Sample, value.var = metric)
            # Rows have to be reverted in order to match heatmap order:
            dr <- dr[nrow(dr):1,]
            write.table(x=dr,
                        file=file.path(outdir,
                                       make.names(paste0('ont_',
                                              ontology,
                                              '__met_',metric,
                                              '.tsv'))),
                        sep='\t',
                        row.names=FALSE)
        }
    }
    
    id_vars <- c("Sample","ID","name","Ontology",
                 #"Index",
                 "signif_binom","signif_hyper","pass_signif_tests", "pass_post_filter_tests","pass_tests","uniqueId","label")
    
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
