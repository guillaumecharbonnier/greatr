#' Look for all ID in enrichment_tables and create a "slim" version of an ontology which have at least one ID matching from "slimList")
#' slimList <- 'out/sed/extract_go_id_from_obo/wget/http/www.geneontology.org/ontology/subsets/goslim_generic.txt'
#' enrichment_tables: 
#' @export
add_slim_ontologies <- function(enrichment_tables,
                                slimList=NULL){
    if (is.null(slimList)){
        print('No slim list provided, using goslim_generic instead')
        slimList <- data(goslim_generic, package='greatr')
    } else  if (tools::file_ext(slimList) == 'txt'){
        slimList <- scan(slimList, character(), quote = "")
    } else if (tools::file_ext(slimList) == 'yaml'){
        print('yaml as slimList is not implemented yet')
        slimList <- yaml::read_yaml(file=slimList)
    }

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

#' Add custom ontologies provided by a yaml file.
#' @param yaml path to a yaml containing custom ontologies, or the list of the loaded yaml itself
add_custom_ontologies <- function(enrichment_tables,
                                  yaml){
    if (typeof(yaml) == "character"){
        custom_ontologies <- yaml::read_yaml(file=yaml)
    } else if (typeof(yaml) == "list"){
        custom_ontologies <- yaml
    }

    if (any(!custom_ontologies$source %in% names(enrichment_tables[[1]]))){
        error('Specified ontology does not exist in enrichement tables. Check source in yaml')
    }
    for (custom_ontology in custom_ontologies$ontologies){
        print(custom_ontology$name)
        for (sample_label in names(enrichment_tables)){
            res <- enrichment_tables[[sample_label]][[1]]
            res$selId <- NA 
            res <- res[FALSE,]

            for (ontology in custom_ontology$source){
                enrichment_table <- enrichment_tables[[sample_label]][[ontology]]

                enrichment_table$selId <- switch(custom_ontology$id_type,
                                 "ID" = enrichment_table[,"ID"],
                                 "name" = enrichment_table[,"name"],
                                 "ID_name" = paste(enrichment_table[,"ID"], 
                                                   enrichment_table[,"name"], 
                                                   sep='_'),
                                 "ontology_ID_name" = paste(ontology,
                                                            enrichment_table[,"ID"], 
                                                            enrichment_table[,"name"],
                                                            sep='_'))

                inCustomList <- enrichment_table$selId %in% custom_ontology$ids

                res <- rbind(res,
                             enrichment_table[inCustomList,])

            }
            # Reordering 
            ordered <- data.frame(selId=custom_ontology$ids)
            res <- merge(ordered, res, by='selId', sort=FALSE)
            res$selId <- NULL

            # Probably obsolete because I always order by list order regardless:
            #if (custom_ontology$order){
            #    res$order <- # Add here ids as factors using custom_ontology$ids as levels. 
            #                d$label <- factor(x=d$uniqueId,
            #                                                            levels=d$uniqueId[hclust_GOterms$order],
            #                                                                                      labels=d$label)
            #}
            if (is.null(custom_ontology$recompute_adjusted_pval)){
                custom_ontology$recompute_adjusted_pval <- FALSE
            }
            if (custom_ontology$recompute_adjusted_pval){
                res <- compute_additional_metrics(res,
                                                  filterMetrics=attributes(enrichment_tables)$filterMetrics,
                                                  filterGreaterLowerThans=attributes(enrichment_tables)$filterGreaterLowerThans,
                                                  filterThresholds=attributes(enrichment_tables)$filterThresholds)
            }
            enrichment_tables[[sample_label]][[custom_ontology$name]] <- res
        }
    }
    return(enrichment_tables)
}
    



#' enrichment_tables should already contains 'pass_signif_tests' column to greatly reduce computational cost of this step as similiarity distance will only be computed on significant terms.
#' Then metrics have to be computed again so postfilter doesnot ruin similarity selection.
#' @export
add_similarity_filtered_ontologies <- function(enrichment_tables){
    for (ontology in c('GO Cellular Component', 'GO Biological Process', 'GO Molecular Function')){
        if (ontology == 'GO Cellular Component'){
            ont <- 'CC'
        } else if (ontology == 'GO Biological Process'){
            ont <- 'BP'
        } else if (ontology == 'GO Molecular Function'){
            ont <- 'MF'
        }

        if (attributes(enrichment_tables)$assembly == 'hg19'){
            org <- 'org.Hs.eg.db'
        } else if (attributes(enrichment_tables)$assembly %in% c('mm9', 'mm10')){
                   org <- 'org.Mm.eg.db'
        } else if (attributes(enrichment_tables)$assembly == 'danRer7'){
            org <- 'org.Dr.eg.db'
        }

        semData <- GOSemSim::godata(org, ont=ont)

        message('Adding best-in-similarity-group ontology for ', ontology)

        for (simCutoff in 1:9/10){
            for (selVar in c('Binom_Fold_Enrichment', 'Binom_Raw_PValue')){
                id_best_in_group <- data.frame(ID = enrichment_tables[[1]][[ontology]][,'ID'], best_in_group = NA, stringsAsFactors=FALSE)
                for (sample_label in names(enrichment_tables)){
                    message('for ', sample_label)
                    enrichment_table <- enrichment_tables[[sample_label]][[ontology]]
                    id_pval <- enrichment_table[enrichment_table$pass_signif_tests,c('ID', 'Binom_Raw_PValue','Binom_Fold_Enrichment')]
                    # Considering edge case where there is no significant term in sample:
                    if (nrow(id_pval) == 0){
                        id_best_in_group[,sample_label] <- FALSE
                    } else {
                        # NOTE: I may also try to ask distance GO1=onlyOneID, GO2=all_id_pval$ID. May be faster, or not.
                        #sim_dist <- GOSemSim::mgoSim(id_pval$ID, id_pval$ID, semData=GO_data, measure="Wang", combine=NULL)

                        id_pval$best_in_group <- unlist(lapply(X=id_pval$ID, FUN=check_if_best_in_simgroup, id_pval = id_pval, semData = semData, simCutoff=simCutoff, selVar=selVar))
                        #id_pval$best_in_group <- unlist(mclapply(X=id_pval$ID, FUN=check_if_best_in_group, id_pval = id_pval, mc.cores = detectCores(all.tests = FALSE, logical = TRUE)))
                        id_pval[,c('Binom_Raw_PValue','Binom_Fold_Enrichment')] <- NULL

                        id_best_in_group <- merge(id_best_in_group, id_pval, by='ID', suffixes=c('',sample_label), all=TRUE)
                    }
                }
                id_best_in_group$best_in_group <- NULL
                id_best_in_group[is.na(id_best_in_group)] <- FALSE
                #browser()
                id_best_in_group$best_in_group_in_any_samples <- apply(X=id_best_in_group[,-1], MARGIN=1, FUN=any)
                ids_to_keep <- id_best_in_group[id_best_in_group$best_in_group_in_any_samples,'ID']

                new_ontology_name <- paste0('Similarity filtered (Best ',selVar, '; Wang ',simCutoff,') ', ontology)
                for (sample_label in names(enrichment_tables)){
                    et <-  enrichment_tables[[sample_label]][[ontology]]
                    enrichment_tables[[sample_label]][[new_ontology_name]] <- et[et$ID %in% ids_to_keep,]
                }

                # here add compute additionnal metrics again so post-filter is updated with similarity filtering.
                enrichment_tables[[sample_label]][[new_ontology_name]] <- compute_additional_metrics(enrichment_tables[[sample_label]][[new_ontology_name]],
                                                                                                     filterMetrics=attributes(enrichment_tables)$filterMetrics,
                                                                                                     filterGreaterLowerThans=attributes(enrichment_tables)$filterGreaterLowerThans,
                                                                                                     filterThresholds=attributes(enrichment_tables)$filterThresholds)
            }
        }
    }
    return(enrichment_tables)
}

#' Check if a GOID is the one from its hierarchical branch with the lowest P-value.
#'
#' @param id_to_check The GOID to check. If this is a synonym, comparison will be done with the main ID for ancestor and offspring.
#' @param id_pval A two-column data.frame with all IDs and P-values. Currently colnames have to be 'ID' and 'Binom_Raw_PValue'.
#' @param selVar The variable to use to select the best in group. 
#' @return Return TRUE or FALSE. If it is an obsolete or unknown ID, return FALSE.
check_if_best_in_simgroup <- function(id_to_check = 'GO:0044424',
                                      id_pval,
                                      simCutoff=0.5,
                                      measure='Wang',
                                      selVar=c('Binom_Fold_Enrichment','Binom_Raw_PValue'),
                                      semData){
    # if id_to_check is a synonym, use main id to look for ancestor and offspring.
    #message('Looking if ', id_to_check, ' is best in group')
    #print(paste('Looking if ', id_to_check, ' is best in similarity group'))
    sim_dist <- GOSemSim::mgoSim(id_to_check, id_pval$ID, semData=semData, measure=measure, combine=NULL)
    ids_to_check <- colnames(sim_dist)[sim_dist > simCutoff]
    if (selVar == 'Binom_Fold_Enrichment'){
        checkTRUE <- identical(max(id_pval[id_pval$ID %in% ids_to_check, selVar]),
                               id_pval[id_pval$ID == id_to_check, selVar])
    } else {
        #print(paste('id_to_check:', id_to_check))
        #if (id_to_check == 'GO:0042611'){browser()}
        #print(paste(min(id_pval[id_pval$ID %in% ids_to_check, selVar])))
        checkTRUE <- identical(min(id_pval[id_pval$ID %in% ids_to_check, selVar]),
                               id_pval[id_pval$ID == id_to_check, selVar])
    }
    #print(checkTRUE)
    return(checkTRUE)
}


