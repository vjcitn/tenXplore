#secLevGen = function( TermSet, choices, rrdfsupp ) {
#   inds = match(choices, TermSet@cleanFrame$clean)
#   set = lapply(inds, function(x) 
#     try(children_URL(TermSet@cleanFrame$url[x], model=getModel(rrdfsupp), world=getWorld(rrdfsupp)), silent=TRUE))
#   chk = sapply(set, function(x) inherits(x, "try-error"))
#   if (any(chk)) set = set[-which(chk)]
#   if (length(set)==0) return(NULL)
#   if (length(set)==1) return(set[[1]])
#   do.call(c, set)
#}

library(shiny)
library(ontoProc)
library(tenXplore)
data("allGOterms", package="ontoProc")
data("CellTypes", package="tenXplore")


vec = slot(CellTypes, "cleanFrame")$clean
clsupp = getOnto("cellOnto")

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
    numericInput("nsamp", "# 10x samples", value=400, min=100, max=1306127, step=100),
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
