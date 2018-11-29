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

add_best_in_go_group_ontology <- function(enrichment_tables){
    for (sample_label in names(enrichment_tables)){
        for (ontology in c('GO Cellular Component', 'GO Biological Process', 'GO Molecular Function'){
             id_pval <- enrichment_tables[[sample_label]][[ontology]][,c('ID', 'Binom_Raw_PValue')]

        as.list(GOCCCHILDREN['GO:0006412']
    }
}

get_best_in_group_ids <- function(id_pval,
                                  ontology=c('GO Cellular Component', 'GO Biological Process', 'GO Molecular Function')){
    ids <- lapply(X=id_pval$ID, FUN=check_if_best_in_group, id_pval = id_pval, ontology = ontology)
    print(ids)
}

check_if_best_in_group <- function(id_to_check = 'GO:0044424',
                                   id_pval,
                                   ontology=c('GO Cellular Component', 'GO Biological Process', 'GO Molecular Function')){

                                   #relationship_to_eval=c('GOCCOFFSPRING','GOCCCHILDREN','GOCCPARENTS','GOCCANCESTOR','GOCCGROUP',
                                   #                       'GOMFOFFSPRING','GOMFCHILDREN','GOMFPARENTS','GOMFANCESTOR','GOMFGROUP',
                                   #                       'GOBPOFFSPRING','GOBPCHILDREN','GOBPPARENTS','GOBPANCESTOR','GOBPGROUP')
                                   #){
    # First use GOSYNONYM to find term GOID then use GOID to get the ANCESTOR and OFFSPRING
    # trycatch errors in this situation :
    #> as.list(GOSYNONYM['GO:1990904'])
    #Error in .checkKeys(value, Lkeys(x), x@ifnotfound) : 
    #  value for "GO:1990904" not found
    # Then use GOTERM instead.
    # Get Ontology fro GOT

    if (ontology == 'GO Cellular Component'){
    id_to_check_ancestor_and_offspring <-  c(id_to_check,
                                             as.list(GOCCANCESTOR[id_to_check])[[id_to_check]],
                                             as.list(GOCCOFFSPRING[id_to_check])[[id_to_check]])
    } else if (ontology == 'GO Biological Process'){
        id_to_check_ancestor_and_offspring <-  c(id_to_check,
                                                 as.list(GOBPANCESTOR[id_to_check])[[id_to_check]],
                                                 as.list(GOBPOFFSPRING[id_to_check])[[id_to_check]])
    } else if (ontology == 'GO Molecular Function'){
        id_to_check_ancestor_and_offspring <-  c(id_to_check,
                                                 as.list(GOMFANCESTOR[id_to_check])[[id_to_check]],
                                                 as.list(GOMFOFFSPRING[id_to_check])[[id_to_check]])
    }

    if (min(id_pval[id_pval$ID %in% id_to_check_ancestor_and_offspring, 'Binom_Raw_PValue']) == id_pval[id_pval$ID == id_to_check, 'Binom_Raw_PValue']){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

