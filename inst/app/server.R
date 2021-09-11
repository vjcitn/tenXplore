library(restfulSE)
library(org.Mm.eg.db)
library(tenXplore)
library(AnnotationDbi)
library(SummarizedExperiment)

# following reflects fact that se1.3M ExperimentHub element
# lacks rowData!

inSE = se1.3M()[,seq_len(100000)]  # emulate se100k
#okSE = se100k()
rowData(inSE) = rowData(inSE)

data("allGOterms", package="ontoProc")
data("CellTypes")
clsupp = getCellOnto()

 server = function(input, output) {
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
   ss = secLevGen(input$topLevel, clsupp)
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
   showNotification(paste("acquiring ", input$nsamp, " records from HDF5 server"), id="acqnote")
   featinds = sort(as.numeric(na.omit(featinds)))  # new, seems to avoid 'bad slice' at HSDS
   dat=t(as.matrix(assay(finse <- inSE[featinds,1:input$nsamp])))
   removeNotification(id="acqnote")
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
  output$godt = renderDataTable( godtSetup()$dataframe )
  output$thedt = renderDataTable({
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(input$topLevel, clsupp)
   allchoices = ss@cleanFrame[,"clean"]
   validate(need(input$secLevel, "select second order term"))
   as.data.frame(ss@cleanFrame)
   })
  output$secLevUI = renderUI({
   showNotification(paste("setting up"), id="setupnote")
   validate(need(input$topLevel, "select top level term"))
   ss = secLevGen(input$topLevel, clsupp)
   removeNotification(id="setupnote")
   validate(need(!is.null(ss), "no subtypes; please choose another cell type above"))
   allchoices = ss@cleanFrame[,"clean"]
   selectInput("secLevel", "Subtype", choices=unname(allchoices), size=4, multiple=TRUE,
      selectize=FALSE, selected=c("dopaminergic neuron", "GABAergic neuron"))
   })
 }
