
apd = function() dirname(dir(system.file("app", package = "tenXplore"), full = TRUE))[1]

#' basic shiny interface to 10x data with ontological setup for cell selection
#' @import ontoProc
#' @import shiny
#' @import AnnotationDbi
#' @import org.Mm.eg.db
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom restfulSE se1.3M
#' @importFrom stats prcomp biplot na.omit
#' @importFrom matrixStats rowSds
#' @return shiny app invocation
#' @note Starts slowly as it sets up connection to HDF Server.
#' @examples
#' tenXplore
#' @export
tenXplore = function() {
  shiny::runApp(apd())
}

