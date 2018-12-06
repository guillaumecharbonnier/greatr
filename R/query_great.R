#' Query GREAT for the list of loaded bed files.
#'
#' @param beds A list of bed files loaded by load_bed function.
#' @param outdir The output directory for plots and data.
#' @param files A comma-separated list of filenames
#' @param labels A comma-separated list of labels to use for filenames
#' @param assembly The genome assembly of the provided input files.
#' @param ontologies The queried ontologies. All available ontologies are queried by default.
#' @return enrichment_tables from GREAT.
#'
#' @export
query_great <- function(beds,
                        outdir='.',
                        saveTables=TRUE,
                        loadTables=TRUE,
                        assembly=c('hg19','mm9','mm10','danRer7'),
                        ontologies=NULL){
    assembly <- match.arg(assembly)
    enrichment_tables <- list()
    enrichment_tables_path <- file.path(outdir, 'enrichment_tables.Rdata')
    if (loadTables & file.exists(enrichment_tables_path)){
        print('Loading existing tables to save time. Set loadTables to FALSE if you want to force new queries to GREAT.')
        load(file=enrichment_tables_path)
    }

    for (sample in names(beds)){
        # Avoid redoing analysis if already done and loaded from save file.
        if (!sample %in% names(enrichment_tables)){
            print(paste0('Running GREAT analysis for ',sample,'.'))
            job = rGREAT::submitGreatJob(beds[[sample]], species = assembly, request_interval=5)
            if (is.null(ontologies)){
                ontologies <- availableOntologies(job)
            }
            enrichment_tables[[sample]] = rGREAT::getEnrichmentTables(job, ontology = ontologies)
            attr(enrichment_tables[[sample]],'n_queried_regions') <-  length(beds[[sample]])
        }
    }
    ## TODO: Check that adding this attr here does not mess downstream analysis.
    attr(enrichment_tables,'assembly') <- assembly

    if (saveTables){
        print('Saving tables to save time in future uses. Set saveTables to FALSE to disable.')
        save(enrichment_tables, file=enrichment_tables_path)
    }
    return(enrichment_tables)
}


