library(testthat)
library(PharmacoGx)
library(PharmacoExploreR)

# Load small PSet for testing
pset <- loadPSet("NCI60_2021")


test_that("correlateExpressionAUC returns a data.frame with expected columns", {
  drug <- drugNames(pset)[1]
  res <- correlateExpressionAUC(pset, sensitivity.measure = "aac_recomputed", drug = drug)
  
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene", "correlation", "p.value") %in% colnames(res)))
})


test_that("plotExprAUC returns a ggplot object", {
  gene <- res$gene[1]
  p <- plotExprAUC(pset = pset, corResults = res, sensitivity.measure = "aac_recomputed",
                   gene = gene, drug = drug)
  expect_s3_class(p, "ggplot")
})


test_that("defineResponseGroups returns a named vector", {
  groups <- defineResponseGroups(pset, drug, sensitivity.measure = "aac_recomputed", method = "median")
  expect_type(groups, "character")
  expect_equal(length(groups), length(cellNames(pset)))
})

test_that("runDiffExpr returns a data.frame with p-values and adjusted p-values", {
  diff_res <- runDiffExpr(pset, groupLabels = groups, mDataType = "rna", method = "t.test")
  expect_s3_class(diff_res, "data.frame")
  expect_true(all(c("gene", "pval", "adj_pval") %in% colnames(diff_res)))
})

test_that("plotGeneBoxplot returns a ggplot object", {
  p <- plotGeneBoxplot(pset, gene_of_interest = "TP53", response_groups = groups, 
                       mDataType = "rna", plotType = "boxplot")
  expect_s3_class(p, "ggplot")
})

test_that("plotDoseResponse returns a ggplot object", {
  cellLines <- rownames(pset@treatmentResponse$raw[, , "Viability"])[1:2]
  p <- plotDoseResponse(pset, drug = drug, sensitivity.measure = "aac_recomputed", cell.lines = cellLines)
  expect_s3_class(p, "ggplot")
})



