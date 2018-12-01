context('Test queries to GREAT')
#
##test_that("enrichment_tables are produced from queried BED files.", {
##          load('beds.Rdata')
##          enrichment_tables <- query_great(beds,
##                                           outdir='.',
##                                           saveTables=TRUE,
##                                           loadTables=FALSE,
##                                           assembly='hg19',
##                                           ontologies=NULL)
##          expect_equal(1, 1)
##})
#
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
#

context("Modifications on enrichment_tables")
## Deprecated: add_similarity_filtered_ontologies should work better for this aim.
#test_that("best_in_group go ontologies can be added to enrichment_tables", {
#          load('enrichment_tables.Rdata')
#          enrichment_tables_with_best_in_go_group <- add_best_in_go_group_ontologies(enrichment_tables)
#          expect_equal(1,1)
#          #expect_equal(length(enrichment_tables$LT) + 3, 
#          #             length(enrichment_tables_with_slim$LT))
#          save(enrichment_tables_with_best_in_go_group, file='enrichment_tables_with_best_in_go_group.Rdata')
#})

context("Slim ontologies can be added to enrichment_tables")

test_that("slim ontologies can be added to enrichment_tables", {
          load('enrichment_tables.Rdata')
          enrichment_tables_with_slim <- add_slim_ontologies(enrichment_tables,
                                                             slimList='goslim_generic.txt')
          expect_equal(length(enrichment_tables$LT) + 3, 
                       length(enrichment_tables_with_slim$LT))
          save(enrichment_tables_with_slim, file='enrichment_tables_with_slim.Rdata')
})

test_that("additional metrics can be added to enrichment_tables", {
          load('enrichment_tables_with_slim.Rdata')
          enrichment_tables_with_additional_metrics <- add_metrics_to_enrichment_tables(enrichment_tables_with_slim,
                                                                                        filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                                                                        filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                                                                        filterThresholds=c('1.5','0.05','0.05','5'))
          expect_equal(length(enrichment_tables_with_slim$LT$'GO Cellular Component') + 19, 
                       length(enrichment_tables_with_additional_metrics$LT$'GO Cellular Component'))
          save(enrichment_tables_with_additional_metrics, file='enrichment_tables_with_additional_metrics.Rdata')
})

test_that("similarity-filtered go ontologies can be added to enrichment_tables", {
           load('enrichment_tables_with_additional_metrics.Rdata')
           enrichment_tables_with_similarity_filtered_ontologies <- add_similarity_filtered_ontologies(enrichment_tables_with_additional_metrics)
           save(enrichment_tables_with_similarity_filtered_ontologies, file='enrichment_tables_with_similarity_filtered_ontologies.Rdata')
           expect_equal(1,1)
})

#test_that("enrichment_table can be converted to ggplot2 format", {
#          load('enrichment_tables_with_additional_metrics.Rdata')
#          data_for_heatmap <- prepare_data_for_heatmap(enrichmentTables=enrichment_tables_with_additional_metrics,
#                                               showedMetric='Binom_Bonf_PValue',
#                                               transformation=c('mlog10','Zscore','x100','unsignifAsNa'),
#                                               filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),                                                                   
#                                               filterThresholds=c(1.5,0.05,0.05,5),
#                                               filterGreaterLowerThans=c('greater','lower','lower','lower'),
#                                               ontology='MSigDB Pathway',
#                                               orderGOTerms=FALSE,
#                                               goLabels='name')
#          save(data_for_heatmap, file='data_for_heatmap.Rdata')
#})
#
#

test_that("enrichment_table can be converted to ggplot2 format", {
          load('enrichment_tables_with_similarity_filtered_ontologies.Rdata')
          data_for_heatmap2 <- prepare_data_for_heatmap2(enrichmentTables = enrichment_tables_with_similarity_filtered_ontologies,
                                                         clusterTermsBy=NULL,
                                                         goLabels='name')
          save(data_for_heatmap2, file='data_for_heatmap2.Rdata')
          expect_equal(1, 1)
})

test_that("melted data can be plotted", {
          load('data_for_heatmap2.Rdata')
          #debug(plot_melted_data2)
          plot_melted_data2(melted = data_for_heatmap2,
                            outpath = 'heatmap.pdf')
          expect_equal(1, 1)

})

test_that("all heatmaps can be produced in one call", {
          load('data_for_heatmap2.Rdata')
          plot_all_heatmaps(d = data_for_heatmap2,
                            device = 'pdf')
          expect_equal(1, 1)
})




