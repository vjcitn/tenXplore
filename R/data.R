#' cellTypes: data.frame with ids and terms
#' @import methods
#' @importFrom methods slot
#' @importFrom utils data
#' @docType data
#' @format TermSet instance
#' @source efo.owl, August 2017, subclasses of \url{http://www.ebi.ac.uk/efo/EFO_0000324}
#' @examples
#' data(CellTypes)
#' head(slot(CellTypes, "cleanFrame"))
"CellTypes"
#' tenx500: serialized full SummarizedExperiment for demonstration
#' @docType data
#' @format SummarizedExperiment instance
#' @source restfulSE se1.3M pared down to 500 samples, assay materialized and assigned
#' @examples
#' data(tenx500)
#' tenx500
"tenx500"
