#' Plot heatmap from data produced by "prepare_data_for_heatmap"
#'
#' @param d melted data produced by "prepare_data_for_heatmap2"
#' @param outdir 
#' @param device set the image format.
#' @export
plot_all_heatmaps <- function(d,
                              outdir='.',
                              device=c('pdf','svg','png')){
    #    # plot all metrics and ontologies heatmap
    #    plot_melted_data2(melted = d,
    #                      outpath = paste0(outdir,'/ont_all__met_all.',device))
    #
    #    # Plot unitary heatmaps
    #    for (metrics in unique(d$metric)){
    #        for (ontologies in unique(d$Ontology)){
    #            print(paste('Plotting ontologies',ontologies,'and metrics',metrics)) 
    #            plot_melted_data2(melted = d[d$metric %in% metrics & d$Ontology %in% ontologies,],
    #                              outpath = file.path(outdir,
    #                                                  make.names(paste0('ont_',
    #                                                                    ontologies,
    #                                                                    '__met_',
    #                                                                    metrics,
    #                                                                    '.',
    #                                                                    device))))
    #        }
    #    }
    #    

    # plot interesting combinations of metrics and ontologies.
    #GO BP with 
    ontologies <- as.character(unique(d$Ontology))
    ontologies_unique <- as.list(ontologies)
    names(ontologies_unique) <- ontologies

    ontologies_groups <- list(all = ontologies,
                              various_go_bp = c('GO Biological Process',
                                                'Slim GO Biological Process',
                                                'Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.5) GO Biological Process',
                                                'Similarity filtered (Best Binom_Raw_PValue; Wang 0.5) GO Biological Process'),
                              various_go_cc = c('GO Cellular Component',
                                                'Slim GO Cellular Component',
                                                'Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.5) GO Cellular Component',
                                                'Similarity filtered (Best Binom_Raw_PValue; Wang 0.5) GO Cellular Component'),
                              various_go_mf = c('GO Molecular Function',
                                                'Slim GO Molecular Function',
                                                'Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.5) GO Molecular Function',
                                                'Similarity filtered (Best Binom_Raw_PValue; Wang 0.5) GO Molecular Function'),
                              bfe_wang1_go = c('Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.1) GO Biological Process',
                                           'Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.1) GO Cellular Component',
                                           'Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.1) GO Molecular Function'),
                              bfe_wang2_go = c('Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.2) GO Biological Process',
                                           'Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.2) GO Cellular Component',
                                           'Similarity filtered (Best Binom_Fold_Enrichment; Wang 0.2) GO Molecular Function'))

    ontologies <- c(ontologies_unique,ontologies_groups)

    metrics <- as.character(unique(d$metric))
    metrics_unique <- as.list(metrics)
    names(metrics_unique) <- metrics

    metrics_groups <- list(all = metrics,
                           filters = c('mlog10 Binom BH PValue',
                                       'mlog10 Hyper BH PValue',
                                       'Binom Fold Enrichment'))
    metrics <- c(metrics_unique, metrics_groups)
    for (metrics_category in names(metrics)){
        for (ontologies_category in names(ontologies)){
            plot_melted_data2(melted = d[d$metric %in% metrics[[metrics_category]] & d$Ontology %in% ontologies[[ontologies_category]],],
                              outpath = file.path(outdir,
                                                  make.names(paste0('ont_',
                                                                    ontologies_category,
                                                                    '__met_',
                                                                    metrics_category,
                                                                    '.',
                                                                    device))))
        }
    }
}

