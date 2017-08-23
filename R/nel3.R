#' @exportClass rrdfSupport
setClass("rrdfSupport", representation(model="ANY", world="ANY"))
setMethod("show", "rrdfSupport", function(object) {
  cat("rrdfSupport instance from tenXplore\n")
  cat(sprintf("there are %d class statements\n", countClasses(object)))
})

#' accessors for rrdfSupport
#' @export
getModel = function(x) x@model
#' @rdname getModel
#' @aliases getWorld
#' @export
getWorld = function(x) x@world

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
     try(children_URL(TermSet@cleanFrame$url[x], model=getModel(rrdfsupp), world=getWorld(rrdfsupp)), silent=TRUE))
   chk = sapply(set, function(x) inherits(x, "try-error"))
   if (any(chk)) set = set[-which(chk)]
   if (length(set)==0) return(NULL)
   if (length(set)==1) return(set[[1]])
   do.call(c, set)
}

nestedL = function(vec, clsupp, inSE) { 
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
    helpText(h3("Ontology-based gene filtering for 10x 1.3 million dataset")),
    helpText(strong("Cell type (top level) selection")," -- on click a dropdown is produced and multiple",
      "cell type names from Cell Ontology are available.  Use command key to add selections individually.  You can also delete selections from the box by reclicking selected items."),
    selectInput("topLevel", "Type", choices=vec, selected="neuron",
         multiple=TRUE, selectize=FALSE, size=4),
    helpText(strong("Cell subtype selection")," -- on click a dropdown is produced and multiple",
      "cell type names from Cell Ontology are available, annotated as subclasses of top level selections."),
    uiOutput("secLevUI"),
    selectInput("trans", "Transformation", choices=c("ident.", "log(x+1)"),
          selected="log(x+1)", selectize=FALSE),
    numericInput("nsamp", "# 10x samples", value=400, min=100, max=3000, step=50),
    numericInput("pc1", "pc axis 1", 1, min=1, max=10, step=1),
    numericInput("pc2", "pc axis 2", 2, min=2, max=10, step=1),
    numericInput("expa", "biplot expand setting", .9, min=.5, max=1.1, step=.1)
    ),
  mainPanel(
   tabsetPanel(
    tabPanel("PCA", 
      plotOutput("pcs"),
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
     helpText(a("foo", href="bar"),p("i don't know"))
     )
    )
   )
  )
 )
 server = function(input, output) {
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
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
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
   dat=t(assay(finse <- inSE[na.omit(featinds),1:input$nsamp]))
   list(finse=finse, data=dat)
   })
  output$counts = renderDataTable({
   strs = getData()
   dat = t(strs$data)
   syms = make.names(rowData(strs$finse)$symbol, unique=TRUE)
   dat = data.frame(gene=syms, rowmn=rowMeans(dat), rowsd=rowSds(dat), dat)
   colnames(dat) = c("gene", "rowmn", "rowsd", as.character(1:input$nsamp))
   dat
   })
  output$pcs = renderPlot({
#   allg = godtSetup()$genes
#   featinds = match(allg, rowData(inSE)$symbol)
#   validate(need(length(featinds)>0, "no expression data for this subtype, please revise"))
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
  output$godt = renderDataTable( godtSetup()$dataframe )
  output$thedt = renderDataTable({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
   allchoices = ss@cleanFrame[,"clean"]
   validate(need(input$secLevel, "select second order term"))
   as.data.frame(ss@cleanFrame)
   })
  output$secLevUI = renderUI({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
   validate(need(!is.null(ss), "no subtypes; please choose another cell type above"))
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
if (!exists("remouse")) remouse <<- se1.3M()
data("allGOterms")
data("CellTypes")
clsupp = buildCellOntSupport()
nestedL(CellTypes@cleanFrame$clean, clsupp, remouse)
}
