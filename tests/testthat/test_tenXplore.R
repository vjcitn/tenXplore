
library(tenXplore)

context("ontology processing")

test_that("children compute", {
  library(ontoProc)
  efoOnto = getEFOOnto()
  cc = children_TAG("EFO:0000324", efoOnto) 
  ss = secLevGen("B cell", efoOnto)
  cl = slot(ss, "cleanFrame")
  expect_true(class(ss) == "TermSet")
  expect_true("mature B cell" %in% cl[,1])
})


