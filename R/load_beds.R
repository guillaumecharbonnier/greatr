#' Load bed files
#'
#' @param outdir The output directory for plots and data.
#' @param outfile The filename of the heatmap. If NULL, a smart name is given based on chosen settings.
#' @param indir The input directory for bed files, will be combined with ``files'' parameter to get full path.
#' @param files A comma-separated list of filenames
#' @param beds A GRangesList object of the bed files to query to GREAT. If provided, ``indir''  and ``files'' are ignored.
#' @param enrichmentTables A list of queried tables. If provided, ``indir'', ``files'' and ``beds'' are ignored.
#' @param sampleLabels A comma-separated list of labels to use for filenames
#' @param smartLabels If TRUE, labels will be taken from filename without extension and dirname.
#' @param showedMetric Define the metric displayed in the heatmap. Can be modified by transformation parameter.
#' @param subsampleSize Specify the number of features to draw from beds. Keep it to NULL and the function will use the number of features of the sample with lowest number of features.
#' @param subsample_replicates The number of drawing replicates for each sample.
#' @param transformation A vector of transformations that can be done on the showed metric. Note that if "mlog10" and "Zscore" are given, then log transformation is always applied before Zscore.
#' @param filterMetric The metric used to limit the number of rows displayed in the heatmap.
#' @param filterThreshold GO terms where at least one sample has its "filterMetric" value below "filterThreshold" will be taken for display in the heatmap.
#' @param assembly The genome assembly of the provided input files.
#' @param ontology The displayed ontology. 'MSigDB Pathway' by default. Choose only one from available choices. Can accept ontology with '_' in name.
#' @param goLabels Use 'ID', 'name' or the combination of both 'ID-name' from ontologies as labels on heatmaps.
#' @return enrichmentTables that can be used to shortcut bed loading and GREAT querying for subsequent use of ``make_great_heatmap'' for different metric and transformation.
#'
#' @export
load_beds <- function(outdir='.',
                      indir='.',
                      files=NULL,
                      sampleLabels=NULL,
                      smartLabels=TRUE,
                      assembly=c('hg19','mm9','mm10','danRer7')){
    dir.create(outdir,
               recursive = TRUE,
               showWarnings = FALSE)
    if (is.null(files)){
        files <- list.files(path=indir,
                   pattern=".*.bed$",
                   full.names = FALSE)
    } else {
        files <- unlist(strsplit(x=files, split=','))
    }
    paths <- file.path(indir,files)
    if (smartLabels){
        labels <- tools::file_path_sans_ext(basename(paths))
    }
    print(paths)
    names(paths) <- labels
    print(paths)
    beds <- list()
    for (sample in names(paths)){
        beds[sample] <- rtracklayer::import(con=as.character(paths[sample]), format='BED', genome=assembly)
        # GREAT does not like non-integer score
        if (!is.null(beds[[sample]]@elementMetadata$score)){
            beds[[sample]]@elementMetadata$score <- round(beds[[sample]]@elementMetadata$score)
        }
        #beds[sample] <- read.table(as.character(paths[sample]), sep='\t')[,1:3]
    }
    beds <- GenomicRanges::GRangesList(beds)
    save(beds, file=file.path(outdir,'beds.Rdata'))
    return(beds)
}


