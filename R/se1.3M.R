#' add/retrieve HSDS-based SE to/from cache
#' @import BiocFileCache SummarizedExperiment
#' @param cache BiocFileCache-like cache
#' @export
se1.3M = function(cache = BiocFileCache::BiocFileCache()) {
  requireNamespace("SummarizedExperiment", quietly=TRUE)
  q = bfcquery(cache, "BiocRestful/tenx1.3M.rds")
  nans = nrow(q)
  if (nans >= 1) return(readRDS(q$rpath[nans]))
  ans = bfcadd(cache, rname="https://mghp.osn.xsede.org/bir190004-bucket01/BiocRestful/tenx1.3M.rds",
     action="copy", download=TRUE)
  readRDS(ans)
}

