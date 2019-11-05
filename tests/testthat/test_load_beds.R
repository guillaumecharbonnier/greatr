context('Load BED')

test_that("BED files are loaded from a list of files", {
          beds <- load_beds(indir=system.file("extdata", package="greatr"),
                            smartLabels=TRUE,
                            assembly='hg19')
          save(beds, file="beds.Rdata")
          beds_out <- beds
          rm(beds)
          data("beds", package="greatr")
          expect_identical(beds_out, beds)
})


