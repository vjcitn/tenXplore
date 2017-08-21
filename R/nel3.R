#' @exportClass rrdfSupport
setClass("rrdfSupport", representation(model="ANY", world="ANY"))

#' read and model the Experimental Factor Ontology as shipped in OWL and parsed in redland
#' @export
buildEFOOntSupport = function() {
 EFworld <- new("World")
 EFstorage <- new("Storage", EFworld, "hashes", name="", options="hash-type='memory'")
 EFmodel <- new("Model", world=EFworld, EFstorage, options="")
 EFparser <- new("Parser", EFworld)
 parseFileIntoModel(EFparser, EFworld, system.file("owl/efo.owl", package="tenXplore"), EFmodel)
 new("rrdfSupport", model=EFmodel, world=EFworld)
}

#' read and model the Cell Ontology as shipped in OWL and parsed in redland
#' @export
buildCellOntSupport = function() {
 CLworld <- new("World")
 CLstorage <- new("Storage", CLworld, "hashes", name="", options="hash-type='memory'")
 CLmodel <- new("Model", world=CLworld, CLstorage, options="")
 CLparser <- new("Parser", CLworld)
 parseFileIntoModel(CLparser, CLworld, system.file("owl/cl.owl", package="tenXplore"), CLmodel)
 new("rrdfSupport", model=CLmodel, world=CLworld)
}

secLevGen = function( TermSet, choices, rrdfsupp ) {
   inds = match(choices, TermSet@cleanFrame$clean)
   set = lapply(inds, function(x) 
     try(children_URL(TermSet@cleanFrame$url[x], model=rrdfsupp@model, world=rrdfsupp@world), silent=TRUE))
   chk = sapply(set, function(x) inherits(x, "try-error"))
   if (any(chk)) set = set[-which(chk)]
   if (length(set)==1) return(set[[1]])
   do.call(c, set)
}

nestedL = function(vec, clsupp, inSE) { 
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText(paste("Cell type (top level) selection -- on click a dropdown is produced and multiple",
      "cell type names from Cell Ontology are available.  Use command key to add selections individually.  You can also delete selections from the box by reclicking selected items.")),
    selectInput("topLevel", "Type", choices=vec, selected="neuron",
         multiple=TRUE, selectize=FALSE, size=4),
    helpText(paste("Cell subtype -- on click a dropdown is produced and multiple",
      "cell type names from Cell Ontology are available, annotated as subclasses of top level selections.")),
    uiOutput("secLevUI"),
    numericInput("pc1", "pc axis 1", 1, min=1, max=10, step=1),
    numericInput("pc2", "pc axis 2", 1, min=2, max=10, step=1)
    ),
  mainPanel(
   tabsetPanel(
    tabPanel("PCA", plotOutput("pcs")),
    tabPanel("agrep",
     dataTableOutput("godt")
     ),
    tabPanel("basics",
     helpText(sprintf("%d top level cats\n", length(vec))),
     textOutput("getSecLNum"),
     dataTableOutput("thedt")
     )
    )
   )
  )
 )
 server = function(input, output) {
  output$getSecLNum = renderText({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
   allchoices = ss@cleanFrame[,"clean"]
   paste(as.character(length(allchoices)), " second level options")
   })
  buildGOTab = reactive({
   validate(need(input$secLevel, "select second order term"))
   curinds = unlist(lapply(input$secLevel, function(x) agrep( x, allGOterms[,2] ) ) )
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
  output$pcs = renderPlot({
   allg = godtSetup()$genes
   data = t(assay(finse <- inSE[na.omit(match(allg, rowData(inSE)$symbol)),1:30]))
   syms = make.names(rowData(finse)$symbol, unique=TRUE)
   colnames(data) = syms
   pcs = prcomp(data)
   biplot(pcs, xlabs=rep("x", nrow(data)), choices=c(input$pc1, input$pc2))
   })
  output$godt = renderDataTable( godtSetup()$dataframe )
  output$thedt = renderDataTable({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
   allchoices = ss@cleanFrame[,"clean"]
   validate(need(input$secLevel, "select second order term"))
   cbind(sec=input$secLevel,as.data.frame(ss@cleanFrame))
   })
  output$secLevUI = renderUI({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
   allchoices = ss@cleanFrame[,"clean"]
   selectInput("secLevel", "Subtype", choices=allchoices, size=4, multiple=TRUE,
      selectize=FALSE, selected=c("dopaminergic neuron", "GABAergic neuron"))
   })
 }
app = shinyApp(ui, server)
runApp(app)
}

#' basic shiny interface to 10x data with ontological setup for cell selection
#' @export
tenXplore = function() {
if (!exists("remouse")) remouse=se1.3M()
data("allGOterms")
data("CellTypes")
clsupp = buildCellOntSupport()
nestedL(CellTypes@cleanFrame$clean, clsupp, remouse)
}
