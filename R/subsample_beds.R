#' Subsampling to equal number of features
#' 
#' @param beds The GRangesList containing the bed file of samples to consider.
#' @param subsample_size Specify the number of features to draw from beds. Keep it to NULL and the function will use the number of features of the sample with lowest number of features.
#' @param subsample_replicates The number of drawing replicates for each sample.
subsample_beds <- function(beds,
                           subsampleSize=NULL,
                           subsampleReplicates=3){
    if (is.null(subsampleSize)){
        subsampleSize <- min(elementNROWS(beds))
    }
    print(paste('subsampling to', subsampleSize))

    for (sample in names(beds)){
        if (elementNROWS(beds[sample]) > subsampleSize){
            print(paste('  for', sample))
            for (subsampleReplicate in 1:subsampleReplicates){
                print(paste('    replicate',subsampleReplicate))
                subsample_name <- paste0(sample,
                                         '_subsample_size-',
                                         subsampleSize,
                                         '_replicate-',
                                         subsampleReplicate)
                beds[[subsample_name]] <- sample(beds[[sample]],
                                                 size = subsampleSize,
                                                 replace = FALSE)
            }
        } else {
            print(paste('  skipping', sample))
        }
    }
    return(beds)
}

