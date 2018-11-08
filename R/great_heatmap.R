library(tools) # file_path_sans_ext
library(rtracklayer) # Loading BED
library(rGREAT) # Queries
library(ggplot2)
library(ggpubr)
library(reshape) # melt
library(scales)
##library(ggdendro)
##library(mvtnorm) 
#library(gplots) # Need heatmap.2
#library(RColorBrewer) # Need a better palette for the heatmap
##library(klaR) # kmode
#library(cluster) # PAM gower

#' Produce heatmaps of ontology enrichment for a list of bed files, or a list of bed objects
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
make_great_heatmap <- function(outdir='.',
                               outfile=NULL,
                               indir='.',
                               files=NULL,
                               beds=NULL,
                               enrichmentTables=NULL,
                               sampleLabels,
                               smartLabels=TRUE,
                               subsampleReplicates=NULL,
                               subsampleSize=NULL,
                               showedMetric=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Binom_Adjp_BH','Binom_Raw_PValue','Hyper_Bonf_PValue','Hyper_Adjp_BH','Hyper_Fold_Enrichment','Hyper_Raw_PValue','Binom_Rank','Hyper_Rank', 'mlog10_Binom_Bonf_PValue', 'mlog10_Hyper_Bonf_PValue','sum_Norm_mlog10_Binom_Bonf_PValue','sum_Norm_mlog10_Hyper_Bonf_PValue','Post_Filter_Binom_Rank','Post_Filter_Hyper_Rank','sum_Signif_Norm_mlog10_Binom_Bonf_PValue','sum_Signif_Norm_mlog10_Hyper_Bonf_PValue','sum_Displayed_Norm_mlog10_Binom_Bonf_PValue','sum_Displayed_Norm_mlog10_Hyper_Bonf_PValue'),
                               transformation=c('mlog10','Zscore','x100','unsignifAsNa'),
                               filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                               filterGreaterLowerThans=c('greater','lower','lower','lower'),
                               filterThresholds=c('1.5','0.05','0.05','10'),
                               assembly=c('hg19','mm9','mm10'),
                               ontology=c('MSigDB Pathway','GO Molecular Function','GO Biological Process','PANTHER Pathway','Disease Ontology',
                                          'MSigDB_Pathway','GO_Molecular_Function','GO_Biological_Process','PANTHER_Pathway','Disease_Ontology'),
                               goLabels=c('name','ID-name','ID')){
    showedMetric <- match.arg(showedMetric) # Default to the first one.
    #filterMetric <- match.arg(filterMetric)
    assembly <- match.arg(assembly)
    ontology <- match.arg(ontology)
    goLabels <- match.arg(goLabels)

    ontology <- gsub(pattern='_',
                     replacement=' ',
                     x=ontology)
    if (is.null(enrichmentTables)){
        if (is.null(beds)){
            beds <- load_beds(outdir,
                              indir,
                              files,
                              assembly='hg19')
            if (!is.null(subsampleReplicates)){
                beds <- subsample_beds(beds,
                                       subsampleSize=subsampleSize,
                                       subsampleReplicates=subsampleReplicates)
            }
        }
        enrichmentTables <- query_great(outdir,
                                         assembly=assembly,
                                         beds=beds)
        enrichmentTables <- add_metrics_to_enrichment_tables(enrichmentTables,
                                                             filterMetrics=filterMetrics,
                                                             filterGreaterLowerThans=filterGreaterLowerThans,
                                                             filterThresholds=filterThresholds)
    }
    data_for_heatmap <- prepare_data_for_heatmap(enrichmentTables=enrichmentTables,
                                                 showedMetric=showedMetric,
                                                 transformation=transformation,
                                                 filterMetrics=filterMetrics,
                                                 filterGreaterLowerThans=filterGreaterLowerThans,
                                                 filterThresholds=filterThresholds,
                                                 ontology=ontology,
                                                 goLabels=goLabels)
    p <- plot_heatmap(data_for_heatmap,
                 outdir=outdir,
                 showedMetric=showedMetric,
                 transformation=transformation)
    great_heatmap <- list(enrichmentTables=enrichmentTables,
                          data_for_heatmap=data_for_heatmap,
                          p=p)
    return(great_heatmap)
}

load_beds <- function(outdir,
                      indir,
                      files,
                      sampleLabels,
                      smartLabels=TRUE,
                      assembly=c('hg19')){
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
    return(beds)
}

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

#' Query GREAT for the list of loaded bed files.
#'
#' @param beds A list of bed files loaded by load_bed function.
#' @param outdir The output directory for plots and data.
#' @param files A comma-separated list of filenames
#' @param labels A comma-separated list of labels to use for filenames
#' @param assembly The genome assembly of the provided input files.
#' @param ontologies The queried ontologies. All listed ontologies are queried by default.
#' @return enrichment_tables from GREAT.
#'
#' @export
query_great <- function(beds,
                        outdir='out/r/great_heatmap',
                        saveTables=TRUE,
                        loadTables=TRUE,
                        assembly=c('hg19','mm9','mm10'),
                        ontologies=c('GO Molecular Function','GO Biological Process','MSigDB Pathway','PANTHER Pathway','Disease Ontology')){
    assembly <- match.arg(assembly)
    enrichment_tables <- list()
    enrichment_tables_path <- file.path(outdir, 'enrichment_tables.Rdata')
    if (loadTables & file.exists(enrichment_tables_path)){
        print('Loading existing tables to save time. Set loadTables to FALSE if you want to force new queries to GREAT.')
        load(file=enrichment_tables_path)
    }

    for (sample in names(beds)){
        # Avoid redoing analysis if already done and loaded from save file.
        if (!sample %in% names(enrichment_tables)){
            print(paste0('Running GREAT analysis for ',sample,'.'))
            job = rGREAT::submitGreatJob(beds[[sample]], species = assembly, request_interval=5)
            enrichment_tables[[sample]] = rGREAT::getEnrichmentTables(job, ontology = ontologies)
        }
    }
    if (saveTables){
        print('Saving tables to save time in future uses. Set saveTables to FALSE to disable.')
        save(enrichment_tables, file=enrichment_tables_path)
    }
    return(enrichment_tables)
}

add_metrics_to_enrichment_tables <- function(enrichment_tables,
                                             filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                             filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                             filterThresholds=c('1.5','0.05','0.05','5')){
    for (sample_label in names(enrichment_tables)){
        for (ontology in names(enrichment_tables[[sample_label]])){
            enrichment_tables[[sample_label]][[ontology]] <- compute_additional_metrics(enrichment_tables[[sample_label]][[ontology]],
                                                                                        filterMetrics=filterMetrics,
                                                                                        filterGreaterLowerThans=filterGreaterLowerThans,
                                                                                        filterThresholds=filterThresholds)
        }
    }
    return(enrichment_tables)
}

#' Add various additional metrics to the provided enrichment table.
compute_additional_metrics <- function(enrichment_table,
                                       filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                       filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                       filterThresholds=c('1.5','0.05','0.05','5')){
    enrichment_table$Binom_Rank <- rank(enrichment_table$Binom_Raw_PValue)
    enrichment_table$Hyper_Rank <- rank(enrichment_table$Hyper_Raw_PValue)
    enrichment_table$Binom_Bonf_PValue <- p.adjust(enrichment_table$Binom_Raw_PValue, method='bonferroni')
    enrichment_table$Hyper_Bonf_PValue <- p.adjust(enrichment_table$Hyper_Raw_PValue, method='bonferroni')
    enrichment_table$mlog10_Binom_Bonf_PValue <- -log10(enrichment_table$Binom_Bonf_PValue)
    enrichment_table$mlog10_Hyper_Bonf_PValue <- -log10(enrichment_table$Hyper_Bonf_PValue)

    filterPostFilter <- grepl(pattern='^Post_Filter', x=filterMetrics)
    filterPreFilter <- !grepl(pattern='^Post_Filter', x=filterMetrics)
    filterPreFilterMetrics <- filterMetrics[filterPreFilter]
    filterPreFilterThresholds <- filterThresholds[filterPreFilter]
    filterPreFilterGreaterLowerThans <- filterGreaterLowerThans[filterPreFilter]
    filterPostFilterMetrics <- filterMetrics[filterPostFilter]
    filterPostFilterThresholds <- filterThresholds[filterPostFilter]
    filterPostFilterGreaterLowerThans <- filterGreaterLowerThans[filterPostFilter]

    # These two tests should be used to set font to bold and or italic in plots:
    enrichment_table$signif_binom <- pass_tests(enrichmentTable=enrichment_table,
                                                filterMetrics='Binom_Adjp_BH',
                                                filterThresholds='0.05',
                                                filterGreaterLowerThans='lower')
    enrichment_table$signif_hyper <- pass_tests(enrichmentTable=enrichment_table,
                                                filterMetrics='Hyper_Adjp_BH',
                                                filterThresholds='0.05',
                                                filterGreaterLowerThans='lower')

    print(paste("debug enrichment_table$signif_binom", str(enrichment_table$signif_binom)))

    enrichment_table$pass_signif_tests <- pass_tests(enrichmentTable=enrichment_table,
                                                            filterMetrics=filterPreFilterMetrics,
                                                            filterThresholds=filterPreFilterThresholds,
                                                            filterGreaterLowerThans=filterPreFilterGreaterLowerThans)
    enrichment_table[enrichment_table$pass_signif_tests, 'Post_Filter_Binom_Rank'] <- rank(enrichment_table[enrichment_table$pass_signif_tests,'Binom_Raw_PValue'])
    enrichment_table[enrichment_table$pass_signif_tests, 'Post_Filter_Hyper_Rank'] <- rank(enrichment_table[enrichment_table$pass_signif_tests,'Hyper_Raw_PValue'])
#    enrichment_table$Post_Filter_Binom_Rank <- ifelse(enrichment_table$pass_signif_tests,
#                                                      rank(enrichment_table[enrichment_table$pass_signif_tests,'Binom_Raw_PValue']),
#                                                      NA)
#    enrichment_table$Post_Filter_Hyper_Rank <- ifelse(enrichment_table$pass_signif_tests,
#                                                      rank(enrichment_table[enrichment_table$pass_signif_tests,'Hyper_Raw_PValue']),
#                                                      NA)
    enrichment_table$sum_Norm_mlog10_Binom_Bonf_PValue <- enrichment_table$mlog10_Binom_Bonf_PValue / sum(enrichment_table$mlog10_Binom_Bonf_PValue, na.rm=T)
    enrichment_table$sum_Norm_mlog10_Hyper_Bonf_PValue <- enrichment_table$mlog10_Hyper_Bonf_PValue / sum(enrichment_table$mlog10_Hyper_Bonf_PValue, na.rm=T)

    enrichment_table$sum_Signif_Norm_mlog10_Binom_Bonf_PValue <- ifelse(enrichment_table$pass_signif_tests,
                                                                        enrichment_table$mlog10_Binom_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_signif_tests,'mlog10_Binom_Bonf_PValue']),
                                                                        NA)
    enrichment_table$sum_Signif_Norm_mlog10_Hyper_Bonf_PValue <- ifelse(enrichment_table$pass_signif_tests,
                                                                        enrichment_table$mlog10_Hyper_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_signif_tests,'mlog10_Hyper_Bonf_PValue']),
                                                                        NA)
    # Set post_filter as true in case there is none of them.
    #enrichment_table$pass_post_filter_tests <- TRUE
    if (length(filterPostFilterMetrics) > 0){
        enrichment_table$pass_post_filter_tests <- pass_tests(enrichmentTable=enrichment_table,
                                                              filterMetrics=filterPostFilterMetrics,
                                                              filterThresholds=filterPostFilterThresholds,
                                                              filterGreaterLowerThans=filterPostFilterGreaterLowerThans)
    } else {
        enrichment_table$pass_post_filter_tests <- TRUE
    }

    enrichment_table$sum_Displayed_Norm_mlog10_Binom_Bonf_PValue <- ifelse(enrichment_table$pass_post_filter_tests,
                                                                           enrichment_table$mlog10_Binom_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_post_filter_tests,'mlog10_Binom_Bonf_PValue']),
                                                                           NA)
    enrichment_table$sum_Displayed_Norm_mlog10_Hyper_Bonf_PValue <- ifelse(enrichment_table$pass_post_filter_tests,
                                                                           enrichment_table$mlog10_Hyper_Bonf_PValue / sum(enrichment_table[enrichment_table$pass_post_filter_tests,'mlog10_Hyper_Bonf_PValue']),
                                                                           NA)


    return(enrichment_table)
}

#' Look for the filterMetrics and return a boolean vector flagging GO failing to pass significance test.
#' @return Could return any metric column from enrichmentTable. I should rewrite this function to work directly on the vector. Later...
pass_tests <- function(enrichmentTable,
                              filterMetrics=c('Binom_Fold_Enrichment','Binom_Adjp_BH','Hyper_Adjp_BH'),
                              filterThresholds=c(2,0.05,0.05),
                              filterGreaterLowerThans=c('greater','lower','lower')){
    filterThresholds <- as.numeric(filterThresholds)
    filterMetricsDf <- enrichmentTable[, filterMetrics, drop=FALSE]
    filterMetricsBool <- filterMetricsDf
    for (filterNumber in 1:length(filterMetrics)){
        threshold <- filterThresholds[filterNumber]
        if (filterGreaterLowerThans[filterNumber] == 'lower'){
            filterMetricsBool[,filterNumber] <- filterMetricsDf[,filterNumber] <= threshold
        } else if (filterGreaterLowerThans[filterNumber] == 'greater'){
            filterMetricsBool[,filterNumber] <- filterMetricsDf[,filterNumber] >= threshold
        }
    }
    d <- apply(filterMetricsBool,1,all)
    d[is.na(d)] <- FALSE
    return(d)
}

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

prepare_data_for_heatmap <- function(enrichmentTables,
                                     showedMetric='Binom_Bonf_PValue',
                                     transformation=c('mlog10','Zscore','x100','unsignifAsNa'),
                                     filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                     filterThresholds=c(1.5,0.05,0.05,5),
                                     filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                     ontology='MSigDB Pathway',
                                     orderGOTerms=FALSE,
                                     goLabels=c('name','ID','ID-name')){
    filterThresholds <- as.numeric(filterThresholds)
    goLabels <- match.arg(goLabels) # Default to the first one.
    data_for_heatmap <- enrichmentTables[[1]][[ontology]]

    print(paste("debug data_for_heatmap", str(data_for_heatmap)))
    col_to_remove_after_merge <- colnames(data_for_heatmap)

    #col_to_remove_after_merge <- col_to_remove_after_merge[! col_to_remove_after_merge %in% c('ID','name')]
    for (sample_label in names(enrichmentTables)){
        to_merge <- enrichmentTables[[sample_label]][[ontology]]
        if ('unsignifAsNa' %in% transformation){
            to_merge$unsignifAsNaShowedMetric <- to_merge[,showedMetric]
            to_merge$unsignifAsNaShowedMetric[!to_merge$pass_signif_tests] <- NA
        }
        data_for_heatmap <- merge(data_for_heatmap,
                                  to_merge,
                                  by=c('ID','name'),
                                  suffix=c('',
                                           paste0('.',
                                                  sample_label)))
    }
    if (goLabels == 'ID-name'){
        labels_for_heatmap <- paste(data_for_heatmap[,'ID'],
                                    data_for_heatmap[,'name'],
                                    sep=' - ')
    } else {
        labels_for_heatmap <- data_for_heatmap[,goLabels]
    }
    data_for_heatmap[,col_to_remove_after_merge] <- NULL
    # Deduplicating labels with suffix number.
    # It occurs for some databases when using only "name".
    labels_for_heatmap <- ifelse(duplicated(labels_for_heatmap) | duplicated(labels_for_heatmap, fromLast=TRUE), 
                                 paste(labels_for_heatmap,
                                       ave(labels_for_heatmap, 
                                           labels_for_heatmap, 
                                           FUN=seq_along),
                                       sep='_'), 
                                 labels_for_heatmap)
    rownames(data_for_heatmap) <- labels_for_heatmap
    print(paste('Number of GO terms before filtering: ',nrow(data_for_heatmap)))
    data_for_heatmap$any_sample_pass_signif_tests <- apply(data_for_heatmap[,grep(pattern="pass_signif_tests",
                                                                                  x=colnames(data_for_heatmap))],
                                                           1,
                                                           any)
    data_for_heatmap$any_sample_pass_post_filter_tests <- apply(data_for_heatmap[,grep(pattern="pass_post_filter_tests",
                                                                                       x=colnames(data_for_heatmap))],
                                                                1,
                                                                any)
    data_for_heatmap <- data_for_heatmap[data_for_heatmap$any_sample_pass_signif_tests & data_for_heatmap$any_sample_pass_post_filter_tests,]
    #for (filterNumber in 1:length(filterMetrics)){
    #    pattern <- paste0('^',
    #                      filterMetrics[filterNumber])
    #    threshold <- filterThresholds[filterNumber]
    #    print(paste('filterMetric pattern:',
    #                  pattern))
    #    print(paste('threshold:',
    #                threshold))

    #    if (filterGreaterLowerThans[filterNumber] == 'lower'){
    #        minFilterMetric <- apply(X=data_for_heatmap[,grep(x=colnames(data_for_heatmap),
    #                                                          pattern=pattern)],
    #                                 MARGIN=1,
    #                                 FUN=min)
    #        data_for_heatmap <- data_for_heatmap[minFilterMetric <= threshold,]
    #        #if (pattern == '^Post_Filter_Binom_Rank'){
    #        #    browser()
    #        #}
    #    } else if (filterGreaterLowerThans[filterNumber] == 'greater'){
    #        maxFilterMetric <- apply(X=data_for_heatmap[,grep(x=colnames(data_for_heatmap),
    #                                                          pattern=pattern)],
    #                                 MARGIN=1,
    #                                 FUN=max)
    #        data_for_heatmap <- data_for_heatmap[maxFilterMetric >= threshold,]
    #    }
    #    print(paste('Number of GO terms after filtering for',
    #                filterMetrics[filterNumber],
    #                filterGreaterLowerThans[filterNumber],
    #                'than',
    #                filterThresholds[filterNumber],
    #                ':',
    #                nrow(data_for_heatmap)))

    #    #print(apply(data_for_heatmap[,grep(pattern='Post_Filter_Binom_Rank', x=colnames(data_for_heatmap))],1,min))
    #    # Computing post filter ranks to be able to select only top N terms if Post_Filter_(Binom|Hyper)_Rank is used in next iteration.
    #    post_filter_ranks <- apply(data_for_heatmap[,grep(pattern='^(Binom|Hyper)_Raw_PValue',
    #                                                           x=colnames(data_for_heatmap))], 
    #                                    MARGIN=2, 
    #                                    FUN=rank)
    #    colnames(post_filter_ranks) <- gsub(pattern='^(Binom|Hyper)_Raw_PValue',
    #                                             replacement="Post_Filter_\\1_Rank",
    #                                             x=colnames(post_filter_ranks))
    #    data_for_heatmap[,grep(pattern='Post_Filter_(Binom|Hyper)_Rank',
    #                           x=colnames(data_for_heatmap))] <- NULL
    #    data_for_heatmap <- cbind(data_for_heatmap,
    #                              post_filter_ranks)
    #}
    showedMetricPattern <- paste0('^',
                                  showedMetric)
    if ('unsignifAsNa' %in% transformation){
        showedMetricPattern <- paste0('^',
                                      'unsignifAsNaShowedMetric')
    }
    print(paste("debug data_for_heatmap",str(data_for_heatmap)))
    signif_binom <- data_for_heatmap[,grep(x=colnames(data_for_heatmap),
                                           pattern="^signif_binom")]
    colnames(signif_binom) <- gsub(pattern=paste0("^signif_binom",'.'),
                                  replacement='',
                                  x=colnames(signif_binom))
    signif_hyper <- data_for_heatmap[,grep(x=colnames(data_for_heatmap),
                                           pattern="^signif_hyper")]
    colnames(signif_hyper) <- gsub(pattern=paste0("^signif_hyper",'.'),
                                  replacement='',
                                  x=colnames(signif_hyper))
    data_for_heatmap <- data_for_heatmap[,grep(x=colnames(data_for_heatmap),
                                               pattern=showedMetricPattern)]
    colnames(data_for_heatmap) <- gsub(pattern=paste0(showedMetricPattern,'.'),
                                       replacement='',
                                       x=colnames(data_for_heatmap))

    data_for_heatmap <- merge_replicates(data_for_heatmap)

    if ('mlog10' %in% transformation){
        data_for_heatmap <- -log10(data_for_heatmap)
        data_for_heatmap[data_for_heatmap == Inf] <- 333
    }
    if ('pseudolog2' %in% transformation){
        data_for_heatmap <- log2(data_for_heatmap + 1)  
        #data_for_heatmap[data_for_heatmap == Inf] <- 333
    }
    if ('Zscore' %in% transformation){
        data_for_heatmap <- t(scale(t(data_for_heatmap)))
        attr(data_for_heatmap,"scaled:center") <- NULL
        attr(data_for_heatmap,"scaled:scale") <- NULL 
        data_for_heatmap[is.nan(data_for_heatmap)] <- 0
        data_for_heatmap <- as.data.frame(data_for_heatmap)
    }
    if ('x100' %in% transformation){
        data_for_heatmap <- 100 * data_for_heatmap
    }
    labels_for_heatmap <- rownames(data_for_heatmap)
    if(orderGOTerms){
        hclust_GOterms <- hclust(dist(data_for_heatmap,
                                      method = "euclidean"),
                                 method = "ward.D" )
        labels_for_heatmap <- factor(x=labels_for_heatmap,
                                     levels=labels_for_heatmap[hclust_GOterms$order])
    }
    metric <- cbind(label=labels_for_heatmap,
                                  data_for_heatmap)
    signif_binom <- cbind(label=labels_for_heatmap,
                          signif_binom)
    signif_hyper <- cbind(label=labels_for_heatmap,
                          signif_hyper)
    data_for_heatmap <- list(metric=metric,
                             signif_binom=signif_binom,
                             signif_hyper=signif_hyper)
    return(data_for_heatmap)
}

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

plot_facet_heatmaps <- function(indir='.',
                                 outdir='.',
                                 files=NULL,
                                 filterMetrics=c('Binom_Fold_Enrichment','Binom_Adjp_BH','Hyper_Adjp_BH','Post_Filter_Binom_Rank'),
                                 filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                 filterThresholds=c('2','0.05','0.05','10'),
                                 subsampleReplicates=NULL,
                                 assembly=c('hg19','mm9','mm10'),
                                 ontology='GO Biological Process'){
    if (ontology=='GOBP'){
        ontology <- 'GO Biological Process'
    }
    if (ontology=='MSDBP'){
        ontology <- 'MSigDB Pathway'
    }
    great_heatmap <- make_great_heatmap(indir=indir,
                                           outdir=outdir,
                                           ontology=ontology,
                                           files=files,
                                           filterMetrics=filterMetrics,
                                           filterGreaterLowerThans=filterGreaterLowerThans,
                                           filterThresholds=filterThresholds,
                                           subsampleReplicates=subsampleReplicates,
                                           showedMetric='Binom_Fold_Enrichment',
                                           transformation='')
    melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted$metric <- 'Binom Fold\nEnrichment'
    melted$scaled <- scales::rescale(melted$value)

    #After the first heatmap we can get an idea of the dimension of the final plot
    nr <- nrow(great_heatmap$data_for_heatmap)
    nc <- ncol(great_heatmap$data_for_heatmap)
    height_heatmap <- 0.2 * nr
    height_top <- 0.5 # I may improve this looking at how much '\n' I have in 'metric' and multiply this number by 0.2
    height_bottom <- 0.2 + 0.15 * max(nchar(colnames(great_heatmap$data_for_heatmap)))
    height <- height_top + height_bottom + height_heatmap
    width_heatmap <- nc*0.25
    width_labels <- 0.15* max(nchar(as.character(great_heatmap$data_for_heatmap$label)))
    # Manually set n_heatmaps here but can be improved later
    n_heatmaps <- 7
    width <- n_heatmaps * width_heatmap + width_labels

    #great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
    #                   outdir=outdir,
    #                   ontology=ontology,
    #                   filterMetrics=filterMetrics,
    #                   filterGreaterLowerThans=filterGreaterLowerThans,
    #                   filterThresholds=filterThresholds,
    #                   showedMetric='Binom_Fold_Enrichment',
    #                   transformation='unsignifAsNa')
    #plots[['Binom_Fold_Enrichment']] <- plots[['Binom_Fold_Enrichment']] + theme(axis.text.y = element_blank())

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Fold_Enrichment',
                       transformation='Zscore')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'Zscore by GO term\nBinom Fold Enrich'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Rank',
                       transformation='')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'Binom Rank'
    append_to_melted$scaled <- scales::rescale(-append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Post_Filter_Binom_Rank',
                       transformation='')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'Post-Filter\nBinom Rank'
    append_to_melted$scaled <- scales::rescale(-append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Adjp_BH',
                       transformation='mlog10')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'mlog10 Binom\nAdjp_BH'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Hyper_Adjp_BH',
                       transformation='mlog10')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- 'mlog10 Hyper\nAdjp_BH'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- '% of sum of mlog10\nBonf PValue in sample'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Signif_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- '% of sum of signif mlog10\nBonf PValue in sample'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Displayed_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    append_to_melted <- reshape::melt(great_heatmap$data_for_heatmap$metric)
    append_to_melted$metric <- '% of sum of displayed mlog10\nBonf PValue in sample'
    append_to_melted$scaled <- scales::rescale(append_to_melted$value)
    append_to_melted$signif_binom <- reshape::melt(great_heatmap$data_for_heatmap$signif_binom)$value
    append_to_melted$signif_hyper <- reshape::melt(great_heatmap$data_for_heatmap$signif_hyper)$value
    melted <- rbind(melted,append_to_melted)

    n_heatmaps <- 4
    print(unique(melted$metric))
    melted_filters <- melted[melted$metric %in% c('Post-Filter\nBinom Rank',
                                                   'mlog10 Binom\nAdjp_BH',
                                                   'mlog10 Hyper\nAdjp_BH',
                                                   'Binom Fold\nEnrichment'),]
    outpath <- file.path(outdir,'facet_filters.pdf')
    plot_melted_data(melted_filters,
                     outpath)
    
    melted_filters <- melted[melted$metric %in% c('Zscore by GO term\nBinom Fold Enrich',
                                                  'Binom Fold\nEnrichment'),]
    outpath <- file.path(outdir,'facet_fold_enrich.pdf')
    plot_melted_data(melted_filters,
                     outpath)

    melted_filters <- melted[melted$metric %in% c('% of sum of displayed mlog10\nBonf PValue in sample',
                                                  '% of sum of signif mlog10\nBonf PValue in sample',
                                                  '% of sum of mlog10\nBonf PValue in sample',
                                                  'Binom Fold\nEnrichment'),]
    outpath <- file.path(outdir,'facet_sum_norm_mlog10.pdf')
    plot_melted_data(melted_filters,
                     outpath)

    outpath <- file.path(outdir,'facet_all.pdf')
    plot_melted_data(melted,
                     outpath)
}

#' Plot melted data. variable should contain samples, label should contain go_labels.
#' value should contain the metric values.
#' Additionaly, a column "metric" can be used to plot in facets (one column for each metric) multiple metrics
#' TODO: A column "ontology" should be addede to plot in facets (one row for each ontology).
#' A column "signif_binom" and "signif_hyper" containing TRUE/FALSE can be added to bold/italicize values when sample/goTerm pair pass these tests. 
#'
#' @export
plot_melted_data <- function(melted,
                             outpath){
    n_samples <- length(unique(melted$variable))
    n_go_terms <- length(unique(melted$label))
    if (is.null(melted$metric)){
        n_metrics <- 1
    } else {
        n_metrics <- length(unique(melted$metric))
    }
    n_char_samples <- max(nchar(as.character(melted$variable)))
    n_char_go_terms <- max(nchar(as.character(melted$label)))
    n_rows_metrics <- 2 # I may improve this looking at how much '\n' I have in 'metric' and multiply this number by 0.2 + 0.1 for spacing.

    width <- 0.25 * n_samples * n_metrics + 0.15 * n_char_go_terms
    height_top <- 0.2 * n_rows_metrics + 0.1
    height_heatmap <- 0.2 * n_go_terms
    height_bottom <- 0.2 + 0.15 * n_char_samples
    height <- height_top + height_heatmap + height_bottom 
    
    p <- ggplot2::ggplot(data = melted, ggplot2::aes(x = variable, y = label))
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
    if ('metric' %in% colnames(melted)){
        p <- p + ggplot2::facet_grid(cols = ggplot2::vars(metric))
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
    ggplot2::ggsave(filename=outpath,
           plot=p, 
           width=width,
           height=height, 
           limitsize=F)
}



#for (palette_name in 'gradient-white-darkblue'){
#    print(paste('Palette:', palette_name))
#    for (metric in names(heatmap_matrix)){
#        print(paste('Metric:', metric))
#        for (filter_threshold in names(heatmap_matrix[[metric]])){
#            print(paste('Filter threshold:', filter_threshold))
#            for (ontology in names(heatmap_matrix[[metric]][[filter_threshold]])){
#                print(paste('Ontology:',ontology))
#                for (sample_category in names(heatmap_matrix[[metric]][[filter_threshold]][[ontology]])){
#                    print(paste('Sample category:',sample_category))
#
#                    path_prefix <- paste0(outdir,
#                                          '/heatmap_metric-',
#                                          metric,
#                                          '/filter-treshold-',
#                                          filter_threshold,
#                                          '/ontology-',
#                                          ontology,
#                                          '/palette-',
#                                          palette_name,
#                                          '_sample_category-',
#                                          sample_category)
#
#
#                    dir.create(dirname(path_prefix), recursive=T)
#                    # Skipping plot if it already exists.
#                    if (!file.exists(paste0(path_prefix,'.pdf'))){
#                        dm <- heatmap_matrix[[metric]][[filter_threshold]][[ontology]][[sample_category]]
#                        nr <- dim(dm)[1]
#                        nc <- dim(dm)[2]
#                        
    #heatmap_matrix <- list()
##heatmap_matrix[[metric]][[filter]][[ontology]]
#for (ontology in ontologies){
#    print(paste('Preparing heatmap matrix for',ontology))
#    # I initialize the dataframe with the whole first sample so merge 'suffixes' are used
#    # for all samples and then I remove the columns for the first samples added here without suffixes.
#    for (sample_category in names(sample_categories)){
#        print(paste('Category:', sample_category))
#        #for (sample in names(states)){
#        merged <- enrichment_tables[[1]][[ontology]]
#        merged <- merged[order(merged$ID),]
#        rownames_to_add_after_merge <- paste(merged$ID, merged$name)
#        col_to_remove_after_merge <- colnames(merged)
#
#
#        for (sample in sample_categories[[sample_category]]){
#            print(paste('sample :',sample))
#
#            sample_to_merge <- enrichment_tables[[sample]][[ontology]]
#            merged <- merge(merged,
#                            sample_to_merge, 
#                            by='ID', 
#                            suffix=c('', 
#                                     paste0('.', sample)))
#        }
#        rownames(merged) <- rownames_to_add_after_merge
#        merged[,col_to_remove_after_merge] <- NULL
#
#        # Gathering indices for filtering
#        filters <- c('Binom_Rank',
#                     'Hyper_Rank',
#                     'Binom_Bonf_PValue',
#                     'Hyper_Bonf_PValue')
#
#        for (filter in filters){
#            merged[,paste0('min_',filter)] <- apply(X=merged[,grep(colnames(merged), pattern=filter)],
#                                                    MARGIN=1,
#                                                    FUN=min)
#        }
#
#        metrics <- c('Binom_Bonf_PValue', 
#                     'Hyper_Bonf_PValue',
#                     'Binom_Fold_Enrichment',
#                     'Hyper_Fold_Enrichment',
#                     'Binom_Rank',
#                     'Hyper_Rank')
#
#        for (metric in metrics){
#            for (filter in filters){
#                if (grepl(filter, pattern='PValue')){
#                    thresholds <- c(0.05, 0.001, 0.000001)
#                } else if (grepl(filter, pattern='Rank')){
#                    thresholds <- c(10,5)
#                }
#
#                for (threshold in thresholds){
#                    indices_selected_terms <- merged[,paste0('min_',filter)] <= threshold
#                    dm <-  merged[indices_selected_terms,
#                                  grep(colnames(merged),
#                                       pattern=metric)]
#                    colnames(dm) <- sub(pattern=paste0(metric,'.'),
#                                        replacement='',
#                                        x=colnames(dm))
#
#                    dm[,paste0('min_',filter)] <- NULL
#                    dm[,paste0('min_',metric)] <- NULL
#
#                    filter_threshold <- paste(filter,
#                                              threshold,
#                                              sep="_")
#
#
#                    if (grepl(metric, pattern='PValue')){
#                        
#                        # This scaling can lead to -Inf values.
#                        # It is useless anyway because Z-score looks better on log transformation.
#                        #heatmap_matrix[[paste0('Zscore_by_sample_',metric)]][[filter_threshold]][[ontology]][[sample_category]] <- scale(dm)
#
#                        #zscore_by_go_term_dm <- t(scale(t(dm)))
#                        # When scale is done on line with the same number, it produces NaN.
#                        # This happen quite a lot when P-val = 0 so I am just manually changing these values to 0.
#                        #zscore_by_go_term_dm[is.nan(zscore_by_go_term_dm)] <- 0
#                        #heatmap_matrix[[paste0('Zscore_by_GOterm_',metric)]][[filter_threshold]][[ontology]][[sample_category]] <- zscore_by_go_term_dm
#                        minuslog10dm <- -log10(dm)
#                        minuslog10dm[minuslog10dm == Inf]  <- 333
#                        heatmap_matrix[[paste0('-log10_',metric)]][[filter_threshold]][[ontology]][[sample_category]] <- minuslog10dm
#
#                        # Z-score looks better on log-transformed data.
#                        heatmap_matrix[[paste0('Zscore_by_sample_-log10_',metric)]][[filter_threshold]][[ontology]][[sample_category]] <- scale(minuslog10dm)
#                        zscore_by_go_term_dm <- t(scale(t(minuslog10dm)))
#                        zscore_by_go_term_dm[is.nan(zscore_by_go_term_dm)] <- 0
#                        heatmap_matrix[[paste0('Zscore_by_GOterm_-log10_',metric)]][[filter_threshold]][[ontology]][[sample_category]] <- zscore_by_go_term_dm
#
#                    } else {
#                        print(paste('Metric should not need log or Z-score transformations', metric))
#                        heatmap_matrix[[metric]][[filter_threshold]][[ontology]][[sample_category]] <- dm
#                    }
#                }
#            }
#        }
#    }
#}
#```
#
### Creating subcategories
#```{r, creating-subcategories, cache=T}
#sample_categories <- list()
#sample_categories[['raw']] <- names(paths)
#
## Removing this sample it has weird GO outputs
### Note afterward on 2018-07-31 22:33:30:
### It is likely that this sample is perfectly fine and that 
### weird ouptuts are actually the one for B_VB_1 but I have
### a bug in the code right now that lead to columns shifted
### the the right.
#to_remove <- 'B_VB_2'
#vec <- sample_categories[['raw']]
#vec <- vec[! vec %in% to_remove]
#sample_categories[['raw_no_B_VB_2']] <- vec
#
## Removing this sample it has low number of active enhancer and may hinder subsampling strategy below.
## This assumption has yet to be checked.
#to_remove <- 'B_VB_1'
#vec <- sample_categories[['raw']]
#vec <- vec[! vec %in% to_remove]
#sample_categories[['raw_no_B_VB_1']] <- vec
#
## Removing both of them because the other two B_VB samples look great anyway.
#outliers_to_remove <- c('B_VB_1','B_VB_2')
#vec <- sample_categories[['raw']]
#vec <- vec[! vec %in% outliers_to_remove]
#sample_categories[['raw_no_B_VB_outliers']] <- vec
#
## Removing these samples because Salva is less interested in them for the article
#uninteresting_to_remove <- c("pCD4_mem_1", 
#               "pCD4_mem_2",
#               "pCD8_naive_mem",
#               "AAM")
#vec <- sample_categories[['raw']]
#vec <- vec[! vec %in% uninteresting_to_remove]
#sample_categories[['raw_no_uninteresting']] <- vec
#
## Removing both uninteresting and outliers
#to_remove <- c(outliers_to_remove,uninteresting_to_remove)
#vec <- sample_categories[['raw']]
#vec <- vec[! vec %in% to_remove]
#sample_categories[['raw_no_uninteresting_and_outliers']] <- vec
#```
#



### Median-merging samples by stage
#```{r median-merging-samples-by-stage, cache=T}
## Fake median for pCD4_naive_1 which has 3439 features:
#enrichment_tables[['pCD4_naive_1_subsample_size-3439_median']] <- enrichment_tables[['pCD4_naive_1']]
#
#patterns <- list()
#patterns[['HSC_subsample_size-3439_median']] <- 'HSC_[0-9]_subsample_size-3439_median'
##patterns[['HSC']] <- 'HSC_[0-9]_subsample_size-3439_median'
##patterns[['CD34_subsample_size-3439_median']] <- 'CD34_subsample_size-3439_median'
##patterns[['EC_subsample_size-3439_median']] <- 'EC_subsample_size-3439_median'
##patterns[['LC_subsample_size-3439_median']] <- 'LC_subsample_size-3439_median'
##patterns[['SP4_subsample_size-3439_median']] <- 'SP4_subsample_size-3439_median'
##patterns[['SP8_subsample_size-3439_median']] <- 'SP8_subsample_size-3439_median'
#patterns[['pCD4_subsample_size-3439_median']] <- 'pCD4_naive_[0-9]_subsample_size-3439_median'
#patterns[['pCD8_subsample_size-3439_median']] <- 'pCD8_naive_[0-9]_subsample_size-3439_median'
#patterns[['B_VB_subsample_size-3439_median']] <- 'B_VB_[0-9]_subsample_size-3439_median'
#patterns[['B_TGC_subsample_size-3439_median']] <- 'B_TGC_[0-9]_subsample_size-3439_median'
#patterns[['NE_subsample_size-3439_median']] <- 'NE_[P|V]B_[0-9]_subsample_size-3439_median'
#patterns[['MONO_subsample_size-3439_median']] <- 'MONO_.*_subsample_size-3439_median'
#patterns[['IM_subsample_size-3439_median']] <- 'IM_[0-9]_subsample_size-3439_median'
#
##First grep to get all first replicate
##pattern="(.*?)_subsample_size-([0-9]+)_replicate-([0-9]+)"
#text <- names(enrichment_tables)
#for (pattern_name in names(patterns)){
#    print(pattern_name)
#    samples_to_merge <- grep(x=text, pattern=patterns[[pattern_name]])
#    enrichment_tables_to_merge <- enrichment_tables[samples_to_merge]
#    for (ontology in names(enrichment_tables[[1]])){
#        print(ontology)
#        #enrichment_tables_to_merge_for_one_ontology <- lapply(enrichment_tables_to_merge,
#        #                                                      `[[`,
#        #                                                      ontology)
#
#        merged <- enrichment_tables_to_merge[[1]][[ontology]]
#        merged <- merged[order(merged$ID),]
#        merged_ID_name <- merged[,c('ID','name')]
#        rownames_to_add_after_merge <- paste(merged$ID, merged$name)
#        col_to_remove_after_merge <- colnames(merged)
#
#
#        for (sample in names(enrichment_tables_to_merge)){
#            merged <- merge(merged,
#                            enrichment_tables_to_merge[[sample]][[ontology]],
#                            by='ID',
#                            suffix=c('',
#                                     paste0('.', sample)))
#        }
#        rownames(merged) <- rownames_to_add_after_merge
#        merged[,col_to_remove_after_merge] <- NULL
#
#        enrichment_tables[[pattern_name]][[ontology]] <- as.data.frame(matrix(0,
#                                                                              nrow=nrow(enrichment_tables[[1]][[ontology]]),
#                                                                              ncol=ncol(enrichment_tables[[1]][[ontology]])))
#        colnames(enrichment_tables[[pattern_name]][[ontology]]) <- colnames(enrichment_tables[[1]][[ontology]])
#        for (column in (names(enrichment_tables[[1]][[ontology]]))){
#            if (column %in% c('ID','name')){
#                enrichment_tables[[pattern_name]][[ontology]][,column] <- merged_ID_name[,column]
#            } else {
#                enrichment_tables[[pattern_name]][[ontology]][,column] <- apply(merged[,grep(colnames(merged),pattern=column)],1,median)
#            }
#        }
#    }
#}
#
#sample_categories[['median_no_uninteresting_and_outliers']] <- c('HSC_subsample_size-3439_median',
#                                                              'CD34_subsample_size-3439_median',
#                                                              'EC_subsample_size-3439_median',
#                                                              'LC_subsample_size-3439_median',
#                                                              'SP4_subsample_size-3439_median',
#                                                              'SP8_subsample_size-3439_median',
#                                                              'pCD4_subsample_size-3439_median',
#                                                              'pCD8_subsample_size-3439_median',
#                                                              'B_VB_subsample_size-3439_median',
#                                                              'B_TGC_subsample_size-3439_median',
#                                                              'NE_subsample_size-3439_median',
#                                                              'MONO_subsample_size-3439_median',
#                                                              'IM_subsample_size-3439_median')
#
#```
#

### Plotting heatmaps
#```{r plot-heatmaps, cache=T}
##sf <- 2
##n_col <- 36
##palette <- list()
##palette[['Blues']] <- colorRampPalette(brewer.pal(9,"Blues"))(n_col)
## I do not like the palette below after testing:
##palette[['BuPu']] <- colorRampPalette(brewer.pal(9,"BuPu"))(n_col)
##
##gg_color_hue <- function(n) {
##    hues = seq(15, 375, length = n + 1)
##    hcl(h = hues, l = 65, c = 100)[1:n]
##}
##
##palette[['gg_default']] <- gg_color_hue(n_col)
#devices <- c('pdf','png')
#for (palette_name in 'gradient-white-darkblue'){
#    print(paste('Palette:', palette_name))
#    for (metric in names(heatmap_matrix)){
#        print(paste('Metric:', metric))
#        for (filter_threshold in names(heatmap_matrix[[metric]])){
#            print(paste('Filter threshold:', filter_threshold))
#            for (ontology in names(heatmap_matrix[[metric]][[filter_threshold]])){
#                print(paste('Ontology:',ontology))
#                for (sample_category in names(heatmap_matrix[[metric]][[filter_threshold]][[ontology]])){
#                    print(paste('Sample category:',sample_category))
#
#                    path_prefix <- paste0(outdir,
#                                          '/heatmap_metric-',
#                                          metric,
#                                          '/filter-treshold-',
#                                          filter_threshold,
#                                          '/ontology-',
#                                          ontology,
#                                          '/palette-',
#                                          palette_name,
#                                          '_sample_category-',
#                                          sample_category)
#
#
#                    dir.create(dirname(path_prefix), recursive=T)
#                    # Skipping plot if it already exists.
#                    if (!file.exists(paste0(path_prefix,'.pdf'))){
#                        dm <- heatmap_matrix[[metric]][[filter_threshold]][[ontology]][[sample_category]]
#                        nr <- dim(dm)[1]
#                        nc <- dim(dm)[2]
#                        
#                        # Totally empirical threshold to prevent segfault which happen for largest heatmaps
#                        # e.g 'all' and 'pbinom' threshold segfault
#                        # whereas 'all and 'rank 10' is ok.
#
#                        cell_limit <- 300*100
#                        if (nc*nr < cell_limit){
#
#                            hclust_GOterms <- hclust( dist(dm, method = "euclidean"), method = "ward.D" )
#
#                            dm <- data.frame(dm)
#                            dm$GOterm <- factor(hclust_GOterms$labels,
#                                                levels=hclust_GOterms$labels[hclust_GOterms$order])
#
#
#
#                            melted <- reshape::melt(dm)
#                            # loose rownames:
#                            #dm <- apply(data.matrix(dm), 2, as.numeric)
#
#                            height <- min(4+0.1*nr,100)
#                            width <- min(nc*0.2+0.1*max(nchar(hclust_GOterms$labels)),100)
#
#
#
#                            p <- ggplot(data = melted, aes(x = variable, y = GOterm))
#                            p <- p + geom_tile(aes(fill = value))
#                            p <- p + scale_fill_gradient(low = "white", high = "darkblue")
#                            # Does not work:
#                            # Error: Continuous value supplied to discrete scale
#                            #p <- p + scale_fill_brewer(palette = palette_name)
#                            p <- p + theme(axis.text.x  = element_text(angle=320, vjust=1, hjust=0),
#                                           legend.position = "bottom")
#                            p <- p + scale_y_discrete(position = "right")
#
#                            #p <- p + scale_fill_gradientn(limits=quantile(dm,c(0.01,0.99)))
#                            p <- p + xlab('Samples')
#                            p <- p + ylab('Gene Ontology Terms')
#                            p <- p + labs(fill = metric)
#
#                            if (sample_category == 'median_no_uninteresting_and_outliers'){
#                                # If labels are updated, use this line to remove plots and allow for update:
#                                # find out/r -name '*median*' -delete
#                                labels <- c('HSC',
#                                            'CD34',
#                                            'EC',
#                                            'LC',
#                                            'SP4',
#                                            'SP8',
#                                            'pCD4',
#                                            'pCD8',
#                                            'B VB',
#                                            'B TGC',
#                                            'NE',
#                                            'MONO',
#                                            'IM')
#                                p <- p + scale_x_discrete(labels=labels)
#                            }
#
#                            for (device in devices){
#                                path <- paste0(path_prefix,'.',device)
#                                ggsave(filename=path, plot=p, width=width, height=height, limitsize=F)
#                            }
#                        } else {
#                            warning('plotting skipped because matrix is too big and likely to cause segfault')
#                        }
#                    }
#                }
#            }
#        }
#    }
#}
#
#```


#' Make multiple heatmaps with standard metrics and their transformation.
#'
#' @export
make_preset_heatmaps <- function(indir='.',
                                 outdir='.',
                                 files=NULL,
                                 filterMetrics=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Hyper_Bonf_PValue','Post_Filter_Binom_Rank'),
                                 filterGreaterLowerThans=c('greater','lower','lower','lower'),
                                 filterThresholds=c('1.5','0.05','0.05','10'),
                                 subsampleReplicates=NULL,
                                 assembly=c('hg19','mm9','mm10'),
                                 ontology='GO Biological Process'){
    if (ontology=='GOBP'){
        ontology <- 'GO Biological Process'
    }
    if (ontology=='MSDBP'){
        ontology <- 'MSigDB Pathway'
    }
    great_heatmap <- make_great_heatmap(indir=indir,
                                           outdir=outdir,
                                           ontology=ontology,
                                           files=files,
                                           filterMetrics=filterMetrics,
                                           filterGreaterLowerThans=filterGreaterLowerThans,
                                           filterThresholds=filterThresholds,
                                           subsampleReplicates=subsampleReplicates,
                                           showedMetric='Binom_Fold_Enrichment',
                                           transformation='')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Fold_Enrichment',
                       transformation='unsignifAsNa')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir, 
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Post_Filter_Binom_Rank',
                       transformation='')
    ## Post_Filter_Binom_Rank is computed after NA transformation.
    ## Thus it can not work in the current implementation.
    ## It is not very useful so I am not going to make effort to make it works.
    #    enrichmentTables <- make_great_heatmap(enrichmentTables=enrichmentTables, 
    #                       outdir=outdir, 
    #                       ontology=ontology,
    #                       filterMetrics=filterMetrics,
    #                       filterGreaterLowerThans=filterGreaterLowerThans,
    #                       filterThresholds=filterThresholds,
    #                       showedMetric='Post_Filter_Binom_Rank',
    #                       transformation='unsignifAsNa')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Bonf_PValue',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Hyper_Bonf_PValue',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Adjp_BH',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Hyper_Adjp_BH',
                       transformation='mlog10')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Bonf_PValue',
                       transformation=c('mlog10','unsignifAsNa'))
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Bonf_PValue',
                       transformation=c('Zscore','mlog10'))
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='Binom_Fold_Enrichment',
                       transformation='Zscore')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue',
                       transformation=c('x100','unsignifAsNa'))
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Signif_Norm_mlog10_Binom_Bonf_PValue',
                       transformation='x100')
    great_heatmap <- make_great_heatmap(enrichmentTables=great_heatmap$enrichmentTables, 
                       outdir=outdir,
                       ontology=ontology,
                       filterMetrics=filterMetrics,
                       filterGreaterLowerThans=filterGreaterLowerThans,
                       filterThresholds=filterThresholds,
                       showedMetric='sum_Signif_Norm_mlog10_Binom_Bonf_PValue',
                       transformation=c('x100','unsignifAsNa'))
}
#make_preset_heatmaps(indir='inp/bed/alex/peaks',
#                     outdir='out/r/great_heatmap/alex2')
#    beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alex', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation='x100')
#    beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alex', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation='pseudolog2')
#    beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alex', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation='Zscore')
#    beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alex', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation=c('pseudolog2','Zscore'))
#
#
#
#
##if(!exists('first_compilation')){
#    beds <- make_great_heatmap(indir='inp/bed/alex/peaks', outdir='out/r/great_heatmap/alex', ontology='GO Biological Process')
##}
##
#first_compilation <- FALSE
#
#
#enrichment_tables <- make_great_heatmap(indir='inp/bed/alex/peaks', outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process')
#beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process',showedMetric='Post_Filter_Binom_Rank', transformation='')
#beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process',showedMetric='Binom_Fold_Enrichment', transformation='')
#beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process',showedMetric='Binom_Bonf_PValue', transformation='mlog10')
#beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation='')
#beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation='pseudolog2')
#beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation='Zscore')
#beds <- make_great_heatmap(beds=beds, outdir='out/r/great_heatmap/alexNosubsampling', ontology='GO Biological Process',showedMetric='sum_Norm_mlog10_Binom_Bonf_PValue', transformation=c('pseudolog2','Zscore'))
#

#(outdir='.',
# outfile=NULL,
# indir='.',
# files=NULL,
# sampleLabels,
# smartLabels=TRUE,
# showedMetric=c('Binom_Fold_Enrichment','Binom_Bonf_PValue','Binom_Raw_PValue','Hyper_Bonf_PValue','Hyper_Fold_Enrichment','Hyper_Raw_PValue','Binom_Rank','Hyper_Rank', 'mlog10_Binom_Bonf_PValue', 'mlog10_Hyper_Bonf_PValue','sum_Norm_mlog10_Binom_Bonf_PValue','sum_Norm_mlog10_Hyper_Bonf_PValue'),
# transformation=c('mlog10','Zscore'),
# filterMetrics=c('Hyper_Fold_Enrichment','Binom_Bonf_PValue','Post_Filter_Binom_Rank'),
# filterGreaterLowerThans=c('greater','lower','lower'),
# filterThresholds=c('2','0.05','5'),
# assembly=c('hg19','mm9','mm10'),
# ontology=c('MSigDB Pathway','GO Molecular Function','GO Biological Process','PANTHER Pathway','Disease Ontology'),
# goLabels=c('name','ID-name','ID')){

