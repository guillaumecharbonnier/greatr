#' Look for all ID in enrichment_tables and create a "slim" version of an ontology which have at least one ID matching from "slimList")
#' slimList <- 'out/sed/extract_go_id_from_obo/wget/http/www.geneontology.org/ontology/subsets/goslim_generic.txt'
#' enrichment_tables: 
add_slim_ontologies <- function(enrichment_tables,
                                slimList){
    slimList <- scan(slimList, character(), quote = "")
    for (sample_label in names(enrichment_tables)){
        for (ontology in names(enrichment_tables[[sample_label]])){
            if (!grepl(pattern='^Slim ', x=ontology)){
                inSlimList <- enrichment_tables[[sample_label]][[ontology]]$ID %in% slimList
                slim_ontology <- enrichment_tables[[sample_label]][[ontology]][inSlimList,]
                if (nrow(slim_ontology) > 0){
                    enrichment_tables[[sample_label]][[paste('Slim', ontology)]] <- slim_ontology
                }
            }
        }
    }
    return(enrichment_tables)
}



