convert_collapse_formula_to_yaml <- function(collapseSamples){
}

#' -c <collapseSamples>          An optional commma-separated list of group of samples to merge for additional analysis. Samples from the same group should be semicolon-separeted, e.g. Group1Sample1;Group1Sample2,Group2Sample1.
collapse_samples <- function(enrichment_tables,
                             yaml_path){
    yaml_content <- yaml::read_yaml(file=yaml_path)
    if (! 'groups' %in% names(yaml_content)){
        print('No groups found in yaml, no sample collapsing done.')
        ets <- enrichment_tables
    } else {
        ets <- list()
        attr(ets,'assembly') <- attributes(enrichment_tables)$assembly
        for (group in names(yaml_content$groups)){
            samples_in_group <- yaml_content[['groups']][[group]]
            if (length(samples_in_group) == 1){
                ets[[group]] <- enrichment_tables[[samples_in_group]]
            } else {
            print(group)
            enrichment_tables_to_merge <- enrichment_tables[samples_in_group]

            d <- reshape2::melt(enrichment_tables_to_merge)
            names(d)[names(d) == 'L2'] <- 'ontology'
            names(d)[names(d) == 'L1'] <- 'sample'
            d$uniqueId <- paste(d$ontology, d$ID, d$name, sep='_')
            mapper <- unique(d[,c('uniqueId','ontology', 'ID', 'name')])

            a <- reshape2::acast(d, formula = uniqueId ~ sample ~ variable)
            b <- apply(a, c(1,3), median)

            g <- merge(mapper, 
                       data.frame(uniqueId=row.names(b), b, stringsAsFactors=FALSE),
                       by='uniqueId')

            #e <- split(g, d$sample)
            #f <- lapply(e, split, unique(d$ontology))
            e <- split(g, g$ontology)
            e <- lapply(e, function(x) x[!(names(x) %in% c("uniqueId", "ontology"))])

            ets[[group]] <- e
            attr(ets[[group]],'n_queried_regions') <- paste(lapply(enrichment_tables_to_merge, function(x) attributes(x)$n_queried_regions), collapse=';')
            }
        }
    }
    return(ets)
}
