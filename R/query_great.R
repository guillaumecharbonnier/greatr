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
            job <- rGREAT::submitGreatJob(beds[[sample]], species = assembly, request_interval=5)
            if (is.null(ontologies)){
                ontologies <- rGREAT::availableOntologies(job)
            }
            enrichment_table <- rGREAT::getEnrichmentTables(job, ontology = ontologies)

            enrichment_tables[[sample]] <- enrichment_table
            attr(enrichment_tables[[sample]],'n_queried_regions') <-  length(beds[[sample]])
            single_sample_path <- paste0(outdir,
                                         '/single_sample/',
                                         sample)
            #browser()
            dir.create(single_sample_path, recursive=T, showWarnings=F)
            width <- 10
            height <- 4 # Test that with this height, labels do not overlap in 'single_sample/LB/GO.Molecular.Function/GO:0003824/regionGeneAssociationGraphs.pdf'
            pdf(file.path(single_sample_path, 'regionGeneAssociationGraphs.pdf'),
                width=width, 
                height=height)
            par(mfrow = c(1, 3))
            res <- plotRegionGeneAssociationGraphs(job)
            dev.off()
            write.table(data.frame(res), file.path(single_sample_path,'regionGeneAssociation.tsv'), col.names = TRUE, sep = "\t")
            for (ontology in names(enrichment_table)){
                single_sample_ontology_path <- file.path(single_sample_path, make.names(ontology))
                dir.create(single_sample_ontology_path, recursive=T, showWarnings=F)
                write.table(enrichment_table[[ontology]], file.path(single_sample_ontology_path,'enrichment_table.tsv'), col.names = TRUE, sep = "\t")
                
                # Require too much queries to GREAT. Slow, not really useful and will ban you.
                #for (id in enrichment_table[[ontology]]$ID){
                #    single_sample_ontology_id_path <- file.path(single_sample_ontology_path, id)
                #    dir.create(single_sample_ontology_id_path, recursive=T, showWarnings=F)
                #    pdf(file.path(single_sample_ontology_id_path, 'regionGeneAssociationGraphs.pdf'),
                #        width=width,
                #        height=height)
                #    par(mfrow = c(1, 3))
                #    res <- plotRegionGeneAssociationGraphs(job, ontology=ontology, termID=id)
                #    dev.off()
                #    write.table(data.frame(res), file.path(single_sample_ontology_id_path,'regionGeneAssociation.tsv'), col.names = TRUE, sep = "\t")
                #}
            }
        }
    }
    attr(enrichment_tables,'assembly') <- assembly

    if (saveTables){
        print('Saving tables to save time in future uses. Set saveTables to FALSE to disable.')
        save(enrichment_tables, file=enrichment_tables_path)
    }
    return(enrichment_tables)
}


