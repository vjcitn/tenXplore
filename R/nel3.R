
nestedL = function(vec, clsupp, inSE, inmem) { 
maxnsamp = ncol(inSE)

 ui = fluidPage(
  tags$head(
    tags$style(HTML("
      .shiny-output-error-validation {
        color: green;
      }
    "))
  ),
  sidebarLayout(
   sidebarPanel(
    helpText(h3("Ontology-based gene filtering for 10x 1.3 million neuron dataset")),
    helpText(strong("Cell type (top level) selection"),
      "Use command key to add selections individually.  You can also delete selections from the box by reclicking selected items."),
    selectInput("topLevel", "Type", choices=vec, selected="neuron",
         multiple=TRUE, selectize=FALSE, size=4),
    helpText(strong("Cell subtype selection"),
      "Cell subtypes obtained using subclasses of top level selections."),
    uiOutput("secLevUI"),
    selectInput("trans", "Transformation", choices=c("ident.", "log(x+1)"),
          selected="log(x+1)", selectize=FALSE),
    numericInput("nsamp", "# 10x samples", value=400, min=100, max=maxnsamp, step=100),
    numericInput("pc1", "pc axis 1", 1, min=1, max=10, step=1),
    numericInput("pc2", "pc axis 2", 2, min=2, max=10, step=1),
    numericInput("expa", "biplot expand setting", .9, min=.5, max=1.1, step=.1)
    ),
  mainPanel(
   tabsetPanel(
    tabPanel("PCA", 
      plotOutput("pcs"),
#      plotlyOutput("pcs2"),
      tableOutput("def")
        ),
    tabPanel("Cell:~:GO",
     dataTableOutput("godt")
     ),
    tabPanel("Subtype",
     helpText(sprintf("%d top level cats\n", length(vec))),
     textOutput("getSecLNum"),
     dataTableOutput("thedt")
     ),
    tabPanel("Counts",
     dataTableOutput("counts")
     ),
    tabPanel("Comments",
     helpText("Source for this package is at the ", a("tenXplore",
      href="http://github.com/vjcitn/tenXplore"), " repository."),
     helpText("10x data are obtained using the ", a("restfulSE", href="http://github.com/vjcitn/restfulSE")," package.  The HDF5 server technology is used to provide data on request; server hosted using resources provided by National Cancer Institute grant U01 CA214846-01"),
     helpText("OWL documents for Cell Ontology and Experimental Factor Ontology are included with tenXplore package; processing using the", a("redland", href="https://cran.r-project.org/package=redland"), "package.")
     )
    )
   )
  )
 )
 server = function(input, output) {
  CellTypes = NULL
  allGOterms = NULL
  data("CellTypes")
  data("allGOterms", package="ontoProc")
  output$def = renderTable( {
     validate(need(input$secLevel, ""))
     goinds = lapply(c(input$secLevel), function(x) agrep( x, allGOterms[,2] ) ) 
     ll = sapply(goinds, length)
     validate(need(sum(unlist(ll))>0, "no matches, please try another term"))
     tabs = lapply(goinds, function(x) allGOterms[x,])
     data.frame(subtype=input$secLevel , nterms=sapply(tabs,nrow))
     } )
  output$getSecLNum = renderText({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen( input$topLevel, clsupp)
   validate(need(!is.null(ss), "no subtypes; please choose another cell type above"))
   allchoices = ss@cleanFrame[,"clean"]
   paste(as.character(length(allchoices)), " second level options")
   })
  buildGOTab = reactive({
   validate(need(input$secLevel, "select second order term"))
   curinds = unlist(lapply(c(input$secLevel), function(x) agrep( x, allGOterms[,2] ) ) )
     ll = sapply(curinds, length)
     validate(need(sum(unlist(ll))>0, "no matches, please try another term"))
   allGOterms[curinds,]
   })
  godtSetup = reactive({
   tab = buildGOTab()
   gotags = tab[,1]
   oktags = keys(org.Mm.eg.db, keytype="GO")
   gotags = intersect(gotags, oktags)
   gtab = AnnotationDbi::select(org.Mm.eg.db, keys=gotags, keytype="GO", columns=c("ENSEMBL", "SYMBOL"))
   allg = unique(gtab[,"SYMBOL"])
   list(dataframe=gtab, genes=allg)
  })
  getData = reactive({
   allg = godtSetup()$genes
   featinds = match(allg, rowData(inSE)$symbol)
   validate(need(length(featinds)>0, "no expression data for this subtype, please revise"))
   if (!inmem) 
       showNotification(paste("acquiring ", input$nsamp, 
            " records from HDF5 server"), id="acqnote")
   dat=t(assay(finse <- inSE[na.omit(featinds),1:input$nsamp]))
   if (!inmem) removeNotification(id="acqnote")
   list(finse=finse, data=dat)
   })
  output$counts = renderDataTable({
   strs = getData()
   dat = t(strs$data)
   syms = make.names(rowData(strs$finse)$symbol, unique=TRUE)
   dat = data.frame(gene=syms, rowmn=rowMeans(dat), rowsd=rowSds(dat), dat)
   colnames(dat) = c("gene", "rowmn", "rowsd", as.character(1:input$nsamp))
   if (ncol(dat)>200) dat= dat[,1:3]
   dat
   })
  output$pcs = renderPlot({
   strs = getData()
   data = strs$data
   finse = strs$finse
   if (input$trans == "log(x+1)") data = log(data+1)
   syms = make.names(rowData(finse)$symbol, unique=TRUE)
   colnames(data) = syms
   pcs = prcomp(data)
   suppressWarnings(  # arrow warnings
     biplot(pcs, xlabs=rep(".", nrow(data)), choices=c(input$pc1, input$pc2),
       expand=input$expa, 
       main=paste("# genes = ", nrow(finse), ", # cells = ", ncol(finse)))
     )
   })
#  output$pcs2 = renderPlotly({
#   strs = getData()
#   data = strs$data
#   finse = strs$finse
#   if (input$trans == "log(x+1)") data = log(data+1)
#   syms = make.names(rowData(finse)$symbol, unique=TRUE)
#   colnames(data) = syms
#   pcs = prcomp(data)
##varname.size=6, alpha=.2, labels.size=1, var.scale=1
#   suppressWarnings(  # arrow warnings
#     ggplotly(ggbiplot(pcs, xlabs=rep(".", nrow(data)), choices=c(input$pc1, input$pc2),
#        varname.size=6, alpha=.2, labels.size=1, var.scale=1) + ggtitle(
#       paste("# genes = ", nrow(finse), ", # cells = ", ncol(finse))))
#     )
#   })
  output$godt = renderDataTable( godtSetup()$dataframe )
  output$thedt = renderDataTable({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(input$topLevel, clsupp)
   allchoices = ss@cleanFrame[,"clean"]
   validate(need(input$secLevel, "select second order term"))
   as.data.frame(ss@cleanFrame)
   })
  output$secLevUI = renderUI({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(input$topLevel, clsupp)
   validate(need(!is.null(ss), "no subtypes; please choose another cell type above"))
   allchoices = unname(ss@cleanFrame[,"clean"])
   selectInput("secLevel", "Subtype", choices=allchoices, size=4, multiple=TRUE,
      selectize=FALSE, selected=c("dopaminergic neuron", "GABAergic neuron"))
   })
 }
app = shinyApp(ui, server)
runApp(app)
}

#' basic shiny interface to 10x data with ontological setup for cell selection
#' @import ontoProc
#' @import shiny
#' @import AnnotationDbi
#' @import org.Mm.eg.db
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom restfulSE RESTfulSummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @importFrom restfulSE se1.3M
#' @importFrom stats prcomp biplot na.omit
#' @importFrom matrixStats rowSds
#' @param remouse optional SummarizedExperiment instance, assay data in memory
#' @return shiny app invocation
#' @examples
#' tenXplore
#' @export
tenXplore = function(remouse) {
  inmem = TRUE
  if (missing(remouse)) {
    remouse = se1.3M()
    inmem = FALSE
    }
  data("allGOterms", package="ontoProc")
  data("CellTypes")
  clsupp = getCellOnto()
  nestedL(slot(CellTypes, "cleanFrame")$clean, clsupp, remouse, inmem)
}
