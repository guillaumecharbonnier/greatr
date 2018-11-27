#' TODO Write function to merge replicates.
#'
#' param data_for_heatmap the matrix
#' subsample_merge function applied to merge replicates.
#' removeReplicates If TRUE, replicates are removed from the matrix.
merge_replicates <- function(data_for_heatmap,
                             subsample_merge=c('median','mean'),
                             removeReplicates=FALSE){
    subsample_merge <- match.arg(subsample_merge)
    pattern <- "(.*?)_replicate-1"
    text <- colnames(data_for_heatmap)
    result <- regmatches(text,
                         regexec(pattern,text))
    sample_subsample_sizes <- sapply(result,`[`,2)
    sample_subsample_sizes <- sample_subsample_sizes[!is.na(sample_subsample_sizes)]
    for (sample_subsample_size in sample_subsample_sizes){
        samples_to_merge <- grep(x=text,
                                 pattern=sample_subsample_size)
        sample_merge_name <- paste(sample_subsample_size,
                                   subsample_merge,
                                   sep='_')
        data_for_heatmap[,sample_merge_name] <- apply(data_for_heatmap[,samples_to_merge],
                                                      1,
                                                      eval(parse(text=subsample_merge)))
        if (removeReplicates){
            data_for_heatmap[,samples_to_merge] <- NULL
        }
    }
    return(data_for_heatmap)
}

