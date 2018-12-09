#' slimList <- 'out/sed/extract_go_id_from_obo/wget/http/www.geneontology.org/ontology/subsets/goslim_generic.txt'
#'
#' @docType data
#'
#' @usage data(goslim_generic)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references Moore et al. (2013) Genetics 195:1077-1086
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23979570}{PubMed})
#'
#' @source \href{http://www.geneontology.org/ontology/subsets/goslim_generic.obo}{Source}
#' then converted to one column format using "sed -n 's/^.*id: GO/GO/p' {input} > {output}"
#'
#' @format A vector of characters containing GO slim generic IDs
"goslim_generic"
