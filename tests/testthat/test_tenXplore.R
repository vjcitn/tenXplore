
library(tenXplore)

context("ontology processing")

test_that("children compute", {
  .efosupp = buildEFOOntSupport()
  cc = children_URL("<http://www.ebi.ac.uk/efo/EFO_0000324>", 
      .efosupp@model, .efosupp@world)
  ss = secLevGen(cc, "B cell", .efosupp)
  cl = slot(ss, "cleanFrame")
  expect_true(class(ss) == "TermSet")
  expect_true("mature B cell" %in% cl[,1])
})


