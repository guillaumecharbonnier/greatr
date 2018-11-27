#' Plot heatmap from data produced by "prepare_data_for_heatmap"
#'
#' @param data_for_heatmap produced by "prepare_data_for_heatmap"
#' @param outfile 
#' @param device set the image format. Will be ignored if outfile is given as its extension will be used instead.
#' @param showedMetric only used to define fillLabel if not specified
#' @param transformation only used to define fillLabel if not specified. If it contain "Zscore", it also default to a more appropriate palette.
#' @param fillLabel The label for the color scale. 
#' @param samplesAsRow If TRUE, samples are rows and GO terms are columns. Need to adjust nc and nr because the layout is bugged right now if TRUE.
plot_heatmap <- function(data_for_heatmap,
                         outdir,
                         outfile=NULL,
                         device=c('pdf','svg','png'),
                         showedMetric,
                         transformation,
                         fillLabel=NULL,
                         samplesAsRows=FALSE,
                         skipIfOutfileExists=FALSE){
    device <- match.arg(device)
    data_for_heatmap <- data_for_heatmap$metric
    nr <- nrow(data_for_heatmap)
    nc <- ncol(data_for_heatmap)
    transformation <- paste(as.vector(transformation),
                            collapse='_')
    # Totally empirical threshold to prevent segfault which happen for largest heatmaps
    # e.g 'all' and 'pbinom' threshold segfault
    # whereas 'all and 'rank 10' is ok.
    if(is.null(outfile)){
        outfile <- paste0('showedMetric-',
                          showedMetric,
                          '_transformation-',
                          transformation,
                          '.',
                          device)
    }

    outpath <- file.path(outdir,outfile)
    if (is.null(fillLabel)){
        fillLabel <- paste(showedMetric, transformation, sep="\n")
    }

    cell_limit <- 300*100
    if (nc*nr > cell_limit){
        warning('plotting skipped because matrix is too big and likely to cause segfault')
    } else {
        dir.create(outdir, recursive=T, showWarnings=F)
        # Skipping plot if it already exists.
        #if (! (file.exists(outfile) & skipIfOutfileExists){

            melted <- reshape::melt(data_for_heatmap)
            # loose rownames:
            #dm <- apply(data.matrix(dm), 2, as.numeric)
            height <- min(4+0.1*nr,100)
            width <- min(nc*0.2+0.1*max(nchar(as.character(data_for_heatmap$label))),100)

            p <- ggplot2::ggplot(data = melted, ggplot2::aes(x = variable, y = label))
            p <- p + ggplot2::geom_tile(ggplot2::aes(fill = value))
            p <- p + ggplot2::geom_text(ggplot2::aes(label = round(value, 1)))

            p <- p + ggplot2::scale_fill_gradient(low = "white", high = "darkblue")
            # Does not work:
            # Error: Continuous value supplied to discrete scale
            #p <- p + scale_fill_brewer(palette = palette_name)
            p <- p + ggplot2::theme(axis.text.x  = ggplot2::element_text(angle=320, vjust=1, hjust=0),
                           legend.position = "bottom")
            p <- p + ggplot2::scale_y_discrete(position = "right")

            #p <- p + scale_fill_gradientn(limits=quantile(dm,c(0.01,0.99)))
            p <- p + ggplot2::xlab('Samples')
            p <- p + ggplot2::ylab('Gene Ontology Terms')
            p <- p + ggplot2::labs(fill = fillLabel)
            if (samplesAsRows){
                p <- p + ggplot2::coord_flip()
            }
            ggplot2::ggsave(filename=outpath, plot=p, width=width, height=height, limitsize=F)
        #}
    }
    return(p)
}

