#' Subset data for heatmap prepared by prepare_data_for_heatmap.
#'
#' @param metrics
#' @param ontologies
subset_data_for_heatmap <- function(d,
                                    metrics='Binom_Fold_Enrichment',
                                    ontologies='MSigDB Pathway'){
    d <- d[d$metric %in% metrics & d$Ontology %in% ontologies,]
    return(d)
}
