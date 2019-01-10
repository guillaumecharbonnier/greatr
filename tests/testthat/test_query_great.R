## file.remove(list.files(pattern='*.(png|pdf|tsv)', recursive=TRUE)); build();test()
#
context('Test queries to GREAT')
###
####test_that("enrichment_tables are produced from queried BED files.", {
####          load('beds.Rdata')
####          enrichment_tables <- query_great(beds,
####                                           outdir='.',
####                                           saveTables=TRUE,
####                                           loadTables=FALSE,
####                                           assembly='hg19',
####                                           ontologies=NULL)
####          expect_equal(1, 1)
####})

#test_that("enrichment_tables can be loaded from saved query.", {
#          load('beds.Rdata')
#          enrichment_tables <- query_great(beds,
#                                           outdir='.',
#                                           saveTables=FALSE,
#                                           loadTables=TRUE,
#                                           assembly='hg19',
#                                           ontologies=NULL)
#          expect_equal(length(enrichment_tables$LT$'GO Cellular Component'),
#                       17)
#          save(enrichment_tables, file='enrichment_tables.Rdata')
#})

context("Post-query sample collapse can be done")
#test_that("Post-query sample collapse can be done", {
#          load("enrichment_tables.Rdata")
#          collapseSamples <- 'LB+LT=L,MONO'
#          enrichment_tables_for_collapsed_samples <- collapse_samples(enrichment_tables,
#                                                                      yaml_path='conf.yaml')
#          expect_equal(1,1)
#          save(enrichment_tables_for_collapsed_samples, file='enrichment_tables_for_collapsed_samples.Rdata')
#
#})
#
#
context("Slim ontologies can be added to enrichment_tables")
#
#test_that("Default slim ontologies can be added to enrichment_tables", {
#          load('enrichment_tables.Rdata')
#          enrichment_tables_with_slim <- add_slim_ontologies(enrichment_tables)
#          expect_equal(length(enrichment_tables$LT) + 3, 
#                       length(enrichment_tables_with_slim$LT))
#          #save(enrichment_tables_with_slim, file='enrichment_tables_with_slim.Rdata')
#})
#
#test_that("Slim ontologies from txt file can be added to enrichment_tables", {
#          load('enrichment_tables.Rdata')
#          enrichment_tables_with_slim <- add_slim_ontologies(enrichment_tables,
#                                                             slimList='goslim_generic.txt')
#          expect_equal(length(enrichment_tables$LT) + 3, 
#                       length(enrichment_tables_with_slim$LT))
#          save(enrichment_tables_with_slim, file='enrichment_tables_with_slim.Rdata')
#})
#
#
#
context("Additional metrics can be added to enrichment_tables")
#
#test_that("Additional metrics can be added to enrichment_tables with default settings", {
#          load('enrichment_tables_with_slim.Rdata')
#          enrichment_tables_with_additional_metrics <- add_metrics_to_enrichment_tables(enrichment_tables_with_slim)
#
#          expect_equal(length(enrichment_tables_with_slim$LT$'GO Cellular Component') + 23, 
#                       length(enrichment_tables_with_additional_metrics$LT$'GO Cellular Component'))
#})

#test_that("Additional metrics can be added to enrichment_tables with default settings for problematic atac set", {
#          load('enrichment_tables_atac.Rdata')
#          enrichment_tables_atac <- enrichment_tables
#          enrichment_tables_atac_with_additional_metrics <- add_metrics_to_enrichment_tables(enrichment_tables_atac)
#
#          expect_equal(length(enrichment_tables_atac[[1]][[1]]) + 25, 
#                       length(enrichment_tables_atac_with_additional_metrics[[1]][[1]]))
#          save(enrichment_tables_atac_with_additional_metrics, file='enrichment_tables_atac_with_additional_metrics.Rdata')
#})

#
#test_that("Additional metrics can be added to enrichment_tables with default settings", {
#          load('enrichment_tables_with_slim.Rdata')
#          enrichment_tables_with_additional_metrics <- add_metrics_to_enrichment_tables(enrichment_tables_with_slim,
#                                                                                        filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Fold_Enrichment_Rank'),
#                                                                                        filterGreaterLowerThans=c('greater','lower','lower','lower'),
#                                                                                        filterThresholds=c(1.5,0.05,0.05,30))
#          expect_equal(length(enrichment_tables_with_slim$LT$'GO Cellular Component') + 23, 
#                       length(enrichment_tables_with_additional_metrics$LT$'GO Cellular Component'))
#})

#test_that("Additional metrics can be added to enrichment_tables without post filter metric", {
#          load('enrichment_tables_with_slim.Rdata')
#          enrichment_tables_with_additional_metrics <- add_metrics_to_enrichment_tables(enrichment_tables_with_slim,
#                                                                                        filterMetrics=c('Binom_Fold_Enrichment','Binom_Adjp_BH','Hyper_Adjp_BH'),
#                                                                                        filterGreaterLowerThans=c('greater','lower','lower'),
#                                                                                        filterThresholds=c(2,0.05,0.05))
#          expect_equal(length(enrichment_tables_with_slim$LT$'GO Cellular Component') + 23, 
#                       length(enrichment_tables_with_additional_metrics$LT$'GO Cellular Component'))
#          save(enrichment_tables_with_additional_metrics, file='enrichment_tables_with_additional_metrics.Rdata')
#})
#
#context("Similarity-filtered go ontologies can be added to enrichment tables")
#test_that("similarity-filtered go ontologies can be added to enrichment_tables", {
#           load('enrichment_tables_with_additional_metrics.Rdata')
#           enrichment_tables_with_similarity_filtered_ontologies <- add_similarity_filtered_ontologies(enrichment_tables_with_additional_metrics)
#           save(enrichment_tables_with_similarity_filtered_ontologies, file='enrichment_tables_with_similarity_filtered_ontologies.Rdata')
#           expect_equal(1,1)
#})
#
#
#test_that("similarity-filtered go ontologies can be added to enrichment_tables", {
#           load('enrichment_tables_atac_with_additional_metrics.Rdata')
#           enrichment_tables_atac_with_similarity_filtered_ontologies <- add_similarity_filtered_ontologies(enrichment_tables_atac_with_additional_metrics)
#           save(enrichment_tables_atac_with_similarity_filtered_ontologies, file='enrichment_tables_atac_with_similarity_filtered_ontologies.Rdata')
#           expect_equal(1,1)
#
#})
#
#context("Custom ontologies can be added to enrichment tables")
#test_that("Custom ontologies can be added to enrichment tables", {
#          load('enrichment_tables_with_similarity_filtered_ontologies.Rdata')
#          enrichment_tables_with_custom_ontologies <- add_custom_ontologies(enrichment_tables_with_similarity_filtered_ontologies,
#                                                                            yaml='custom_ontologies.yaml')
#          save(enrichment_tables_with_custom_ontologies, file='enrichment_tables_with_custom_ontologies.Rdata')
#          expect_equal(1,1)
#})
#
###test_that("enrichment_table can be converted to ggplot2 format", {
###          load('enrichment_tables_with_additional_metrics.Rdata')
###          data_for_heatmap <- prepare_data_for_heatmap(enrichmentTables=enrichment_tables_with_additional_metrics,
###                                               showedMetric='Binom_Bonf_PValue',
###                                               transformation=c('mlog10','Zscore','x100','unsignifAsNa'),
###                                               filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),                                                                   
###                                               filterThresholds=c(1.5,0.05,0.05,5),
###                                               filterGreaterLowerThans=c('greater','lower','lower','lower'),
###                                               ontology='MSigDB Pathway',
###                                               orderGOTerms=FALSE,
###                                               goLabels='name')
###          save(data_for_heatmap, file='data_for_heatmap.Rdata')
###})
###
###

context("Conversion to ggplot2 format")

test_that("enrichment_table can be converted to ggplot2 format", {
          load('enrichment_tables_with_custom_ontologies.Rdata')
          data_for_heatmap2 <- prepare_data_for_heatmap2(enrichmentTables = enrichment_tables_with_custom_ontologies,
                                                         goLabels='name')
          save(data_for_heatmap2, file='data_for_heatmap2.Rdata')
          expect_equal(1, 1)
})

#test_that("enrichment_table can be converted to ggplot2 format", {
#          load('enrichment_tables_with_custom_ontologies.Rdata')
#          data_for_heatmap2_hclust_bfe <- prepare_data_for_heatmap2(enrichmentTables = enrichment_tables_with_custom_ontologies,
#                                                         clusterTermsBy="Binom_Fold_Enrichment",
#                                                         outdir='clusterBy_bfe_test',
#                                                         goLabels='name')
#          #save(data_for_heatmap2_hclust_zbfe, file='data_for_heatmap2_hclust_zbfe.Rdata')
#          expect_equal(1, 1)
#})
#
#test_that("enrichment_table can be converted to ggplot2 format", {
#          load('enrichment_tables_with_custom_ontologies.Rdata')
#          data_for_heatmap2_hclust_bfe <- prepare_data_for_heatmap2(enrichmentTables = enrichment_tables_with_custom_ontologies,
#                                                         clusterTermsBy="Binom_Fold_Enrichment",
#                                                         goLabels='name')
#          save(data_for_heatmap2_hclust_bfe, file='data_for_heatmap2_hclust_bfe.Rdata')
#          expect_equal(1, 1)
#})
#
#
#test_that("enrichment_table can be converted to ggplot2 format", {
#          load('enrichment_tables_with_similarity_filtered_ontologies.Rdata')
#          data_for_heatmap2 <- prepare_data_for_heatmap2(enrichmentTables = enrichment_tables_with_similarity_filtered_ontologies,
#                                                         clusterTermsBy=NULL,
#                                                         goLabels='name')
#          save(data_for_heatmap2, file='data_for_heatmap2.Rdata')
#          expect_equal(1, 1)
#})
#
##test_that("melted data can be plotted", {
##          load('data_for_heatmap2.Rdata')
##          #debug(plot_melted_data2)
##          plot_melted_data2(melted = data_for_heatmap2,
##                            outpath = 'heatmap.pdf')
##          expect_equal(1, 1)
##
##})
##

context("Plot heatmaps")

#test_that("all heatmaps can be produced in one call", {
#          load('data_for_heatmap2.Rdata')
#          plot_all_heatmaps(d = data_for_heatmap2,
#                            device = 'pdf')
#          expect_equal(1, 1)
#})

#test_that("all heatmaps can be produced in one call", {
#          load('data_for_heatmap2_hclust_zbfe.Rdata')
#          plot_all_heatmaps(d = data_for_heatmap2_hclust_zbfe,
#                            device = 'pdf')
#          expect_equal(1, 1)
#})
#
#test_that("all heatmaps can be produced in one call", {
#          load('data_for_heatmap2_hclust_bfe.Rdata')
#          plot_all_heatmaps(d = data_for_heatmap2_hclust_bfe,
#                            device = 'pdf')
#          expect_equal(1, 1)
#})
#
#
## Should solve this issue for big files
##cairo error 'invalid value (typically too big) for the size of the input (surface, pattern, etc.)'
##test_query_great.R:114: error: all heatmaps can be produced in one call
##test_that("all heatmaps can be produced in one call", {
##          load('data_for_heatmap2.Rdata')
##          plot_all_heatmaps(d = data_for_heatmap2,
##                            device = 'png')
##          expect_equal(1, 1)
##})
#
context("run complete pipeline")
test_that("run_pipeline", {
          run_pipeline(indir='.',
                       outdir='run_pipeline',
                       assembly='hg19',
                       yaml='conf.yaml')
          expect_equal(1, 1)
})

