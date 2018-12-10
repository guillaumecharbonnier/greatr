#' Plot melted data. variable should contain samples, label should contain go_labels.
#' value should contain the metric values.
#' Additionaly, a column "metric" can be used to plot in facets (one column for each metric) multiple metrics
#' TODO: A column "ontology" should be added to plot in facets (one row for each ontology).
#' A column "signif_binom" and "signif_hyper" containing TRUE/FALSE can be added to bold/italicize values when sample/goTerm pair pass these tests. 
#' 
#'
#' @export
plot_melted_data2 <- function(melted,
                              outpath){

    n_samples <- length(unique(melted$Sample))
    n_go_terms <- length(unique(melted$label))
    if (is.null(melted$metric)){
        n_metrics <- 1
    } else {
        n_metrics <- length(unique(melted$metric))
    }
    n_char_samples <- max(nchar(as.character(melted$Sample)))
    n_char_go_terms <- max(nchar(as.character(melted$label)))
    n_rows_metrics <- 2 # I may improve this looking at how much '\n' I have in 'metric' and multiply this number by 0.2 + 0.1 for spacing.

    width_ontology_name <- 1
    width_inter_metrics <- 0.1
    width <- (0.2 * n_samples + width_inter_metrics) * n_metrics + width_ontology_name + 0.09 * n_char_go_terms
    height_top <- 0.2 * n_rows_metrics + 0.1
    height_heatmap <- 0.2 * n_go_terms
    height_bottom <- 0.2 + 0.15 * n_char_samples
    height <- height_top + height_heatmap + height_bottom 
    
    p <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Sample, y = label))
    p <- p + ggplot2::geom_tile(ggplot2::aes(fill = scaled))
    p <- p + ggplot2::geom_text(ggplot2::aes(color = scaled,
                           label = round(value, 1),
                           fontface = ifelse(signif_hyper,
                                             ifelse(signif_binom,
                                                    "bold.italic",
                                                    "bold"),
                                             ifelse(signif_binom,
                                                    "italic",
                                                    "plain")
                                             )
                           )
    )
    if ('metric' %in% colnames(melted) & ! 'Ontology' %in% colnames(melted)){
        p <- p + ggplot2::facet_grid(. ~ metric) #cols = ggplot2::vars(metric))
    }
    if ('Ontology' %in% colnames(melted) & ! 'metric' %in% colnames(melted)){
        p <- p + ggplot2::facet_grid(Ontology ~ ., scales="free_y", space="free_y")
    }
    if ('Ontology' %in% colnames(melted) & 'metric' %in% colnames(melted)){
        p <- p + ggplot2::facet_grid(Ontology ~ metric, 
                                     scales="free_y", 
                                     space="free_y",
                                     labeller = ggplot2::label_wrap_gen(22)) #With 16, Panther PATHWAY stay on 1 line. Used to be 12.
    }
    #p <- p + scale_fill_gradient(low = "white", high = "darkblue", guide=FALSE)
    #p <- p + scale_color_gradient2(low = "black", mid="yellow", high = "white", midpoint=0.5, guide=FALSE)
    #p <- p + scale_fill_gradientn(colors=c('#C3C3C5','#ABABDA','#5B5BC8','#2121D8','#0000BE'),guide=FALSE)
    #p <- p + scale_color_gradientn(colors=c('#FFBF00','#FFBF01','#FFD65C','#FFEFBE','#FFFEFB'),guide=FALSE)
    #p <- p + scale_fill_gradientn(colors=c('#CACACA','#667F93','#2A506D','#052B48','#001322'),guide=FALSE)
    #p <- p + scale_color_gradientn(colors=c('#351F00','#704000','#AA7A3A','#E4C499','#FFFFFF'),guide=FALSE)
    p <- p + ggplot2::scale_fill_gradientn(colors=c('#FFFFFF','#C3DAFF','#3986FF','#00307D','#00060E'),guide=FALSE)
    p <- p + ggplot2::scale_color_gradientn(colors=c('#160E00','#160E00','#FFB025','#FFE7BD','#FFFFFF'),guide=FALSE)

    # Does not work:
    # Error: Continuous value supplied to discrete scale
    #p <- p + scale_fill_brewer(palette = palette_name)
    p <- p + ggplot2::theme(axis.text.x  = ggplot2::element_text(angle=320, vjust=1, hjust=0),
                   legend.position = "bottom")
    p <- p + ggplot2::scale_y_discrete(position = "right")

    #p <- p + scale_fill_gradientn(limits=quantile(dm,c(0.01,0.99)))
    p <- p + ggplot2::xlab('Samples')
    p <- p + ggplot2::ylab('Gene Ontology Terms')
    #p <- p + labs(fill = fillLabel)
    print(paste('width height:', width, height))
    print(unique(melted$Ontology))
    print(unique(melted$metric))

    ggplot2::ggsave(filename=outpath,
           plot=p, 
           width=width,
           height=height, 
           limitsize=F)
}



