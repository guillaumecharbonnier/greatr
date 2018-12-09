#' A function to run all steps provided by greatr package. Is called by greatr command line.
#'
#' @param indir                    Input directory. [default: .]
#' @param files                    Input sample bed files as comma-separated list. Combined with <indir>. If not provided, all bed files inside <indir> will be taken.
#' @param background               An optional additional bed file containing regions not found in <files> to add to the background model.
#' @param outdir                   Output directory. If enrichment_tables.Rdata is present from a previous or interrupted run, new input files will be appended and the processing of existing ones will be skipped. All plots will be redrawn. [default: .]
#' @param assembly                 Assembly (hg19,mm9,mm10, danRer7) [default: hg19]
#' @param collapseSamples          An optional commma-separated list of group of samples to merge for additional analysis. Samples from the same group should be semicolon-separeted, e.g. Group1Sample1;Group1Sample2,Group2Sample1.
#' @param subsampleReplicates      Number of replicates for subsampling, no subsampling if unspecified.
#' @param filterMetrics          Comma-separated list of metrics to use for filtering. Defaults try to match GREAT default settings. See together <filterGreaterLowerThans> and <filterThresholds>. [default: Binom_Fold_Enrichment,Binom_Adjp_BH,Hyper_Adjp_BH]
#' @param filterGreaterLowerThans  Vector of 'lower' or 'greater' to apply fiterThresholds values to filterMetrics. [default: greater,lower,lower]
#' @param filterThresholds         Vector of thresholds to apply to filterMetrics. Values matching threshold are kept, i.e. '<=' and '>=' are used for comparison. [default: 2,0.05,0.05]
#' @param slimList                 An optional yaml file containing IDs to limit analysis for each ontology in order to keep output readable and non-redundant. Defaults use GO Slim Generic (http://www.geneontology.org).
#' @export

run_pipeline <- function(indir='.',
                         files=NULL,
                         background=NULL,
                         outdir='.',
                         assembly=c('hg19','mm9','mm10','danRer7'),
                         collapseSamples=NULL,
                         subsampleReplicates=NULL,
                         filterMetrics=c('Binom_Fold_Enrichment','Binom_Adjp_BH','Hyper_Adjp_BH'),
                         filterGreaterLowerThans=c('greater','lower','lower'),
                         filterThresholds=c('2','0.05','0.05'),
                         slimList=NULL){
    beds <- load_beds(indir=indir,
                      outdir=outdir,
                      files=files,
                      smartLabels=TRUE,
                      assembly=assembly)

    enrichment_tables <- query_great(beds=beds,
                                     outdir=outdir,
                                     saveTables=TRUE,
                                     loadTables=FALSE,
                                     assembly=assembly,
                                     ontologies=NULL)
    save(enrichment_tables, file=file.path(outdir,'enrichment_tables.Rdata'))

    enrichment_tables_with_slim <- add_slim_ontologies(enrichment_tables,
                                                       slimList='goslim_generic.txt')

    enrichment_tables_with_additional_metrics <- add_metrics_to_enrichment_tables(enrichment_tables_with_slim,
                                                                                  filterMetrics=c('Binom_Fold_Enrichment','Binom_Adjp_BH','Hyper_Adjp_BH'),
                                                                                  filterGreaterLowerThans=c('greater','lower','lower'),
                                                                                  filterThresholds=c(2,0.05,0.05))
    
    enrichment_tables_with_similarity_filtered_ontologies <- add_similarity_filtered_ontologies(enrichment_tables_with_additional_metrics)

    data_for_heatmap2 <- prepare_data_for_heatmap2(enrichmentTables = enrichment_tables_with_similarity_filtered_ontologies,
                                                   outdir=outdir,
                                                   clusterTermsBy=NULL,
                                                   goLabels='name')

    plot_all_heatmaps(d=data_for_heatmap2,
                      outdir=outdir,
                      device='pdf')
    # Should solve this issue for big files
    #cairo error 'invalid value (typically too big) for the size of the input (surface, pattern, etc.)'
    #test_query_great.R:114: error: all heatmaps can be produced in one call
    #test_that("all heatmaps can be produced in one call", {
    #          load('data_for_heatmap2.Rdata')
    #          plot_all_heatmaps(d = data_for_heatmap2,
    #                            device = 'png')
    #          expect_equal(1, 1)
    #})
}


