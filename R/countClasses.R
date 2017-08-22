
countClasses = function(co) {
 biq = new("Query", co@world, 
   sprintf("SELECT ?a WHERE {?a %s %s}",
     "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>",
     "<http://www.w3.org/2002/07/owl#Class>"),
       base_uri=NULL, query_language="sparql", query_uri=NULL)
   xx = executeQuery(biq, co@model)
   nc = 0
   while (!is.null(zz <- getNextResult(xx))) {
     if (substr(zz, 1, 3) != "_:r") nc=nc+1
     }
   nc
}

