library(restfulSE)
library(tenXplore)

data("allGOterms")
data("CellTypes")
clsupp = buildCellOntSupport()

vec = CellTypes@cleanFrame$clean

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
    helpText(strong("Cell type (top level) selection"),
      "Use command key to add selections individually.  You can also delete selections from the box by reclicking selected items."),
    selectInput("topLevel", "Type", choices=vec, selected="neuron",
         multiple=TRUE, selectize=FALSE, size=4),
    helpText(strong("Cell subtype selection"),
      "Cell subtypes obtained using subclasses of top level selections."),
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
