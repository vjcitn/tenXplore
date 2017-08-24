
library(restfulSE)
library(tenXplore)
data("allGOterms")
data("CellTypes")
clsupp = buildCellOntSupport()

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

 server = function(input, output) {
 showNotification(paste("establishing HDF5 server interface"), id="intfnote")
 if (!exists("inSE")) inSE <<- se1.3M()
 removeNotification(id="intfnote")
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
   showNotification(paste("acquiring ", input$nsamp, " records from HDF5 server"), id="acqnote")
   dat=t(assay(finse <- inSE[na.omit(featinds),1:input$nsamp]))
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
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
   allchoices = ss@cleanFrame[,"clean"]
   validate(need(input$secLevel, "select second order term"))
   as.data.frame(ss@cleanFrame)
   })
  output$secLevUI = renderUI({
   validate(need(input$topLevel, "select top level term"))
   showNotification("acquiring subtype terms", id="subtnote")
   ss = secLevGen(CellTypes, input$topLevel, clsupp)
   removeNotification(id="subtnote")
   validate(need(!is.null(ss), "no subtypes; please choose another cell type above"))
   allchoices = ss@cleanFrame[,"clean"]
   selectInput("secLevel", "Subtype", choices=allchoices, size=4, multiple=TRUE,
      selectize=FALSE, selected=c("dopaminergic neuron", "GABAergic neuron"))
   })
 }
