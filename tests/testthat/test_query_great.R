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
#

test_that("enrichment_tables can be loaded from saved query.", {
              data(beds, package="greatr")
              enrichment_tables <- query_great(beds,
                                               outdir='.',
                                               saveTables=FALSE,
                                               loadTables=FALSE,
                                               assembly='hg19',
                                               ontologies=NULL)
              expect_equal(length(enrichment_tables$LT$'GO Cellular Component'), 17)
              save(enrichment_tables, file='enrichment_tables.Rdata')
              enrichment_tables_out <- enrichment_tables
              rm(enrichment_tables)
              data("enrichment_tables", package="greatr")
              expect_identical(enrichment_tables, enrichment_tables_out)
})

