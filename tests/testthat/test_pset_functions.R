library(testthat)
library(PharmacoGx)
library(PharmacoExploreR)

# Skip tests if PharmacoGx data can't be downloaded
skip_if_no_pharmacogx <- function() {
  skip_if_not_installed("PharmacoGx")
  skip_on_cran()
}

# Create a shared fixture for tests (runs once)
# Note: In practice, you might want to use a smaller/cached dataset
# or create mock data to speed up tests
setup_test_pset <- function() {
  skip_if_no_pharmacogx()
  
  # Load PSet (this will be slow the first time)
  # Consider caching this or using a smaller subset
  pset <- tryCatch({
    downloadPSet("NCI60_2021")
  }, error = function(e) {
    skip("Could not download PSet for testing")
  })
  
  return(pset)
}

# Test correlateExpressionAUC
test_that("correlateExpressionAUC returns a data.frame with expected columns", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  res <- correlateExpressionAUC(
    pset, 
    drug = drug,
    mDataType = "rna",
    sensitivity.measure = "aac_recomputed"
  )
  
  # Check it's a data frame
  expect_s3_class(res, "data.frame")
  
  # Check expected columns (match your actual output!)
  expect_true(all(c("gene", "cor", "pval", "adj_pval") %in% colnames(res)))
  
  # Check dimensions
  expect_gt(nrow(res), 0)
  expect_equal(ncol(res), 4)
  
  # Check data types
  expect_type(res$gene, "character")
  expect_type(res$cor, "double")
  expect_type(res$pval, "double")
  expect_type(res$adj_pval, "double")
  
  # Check correlation values are in valid range
  expect_true(all(abs(res$cor[!is.na(res$cor)]) <= 1))
  
  # Check p-values are in valid range
  expect_true(all(res$pval[!is.na(res$pval)] >= 0 & 
                    res$pval[!is.na(res$pval)] <= 1))
})

# Test plotExprAUC
test_that("plotExprAUC returns a ggplot object", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  # Get correlation results first
  res <- correlateExpressionAUC(
    pset, 
    drug = drug,
    sensitivity.measure = "aac_recomputed"
  )
  
  # Select a gene that has non-NA correlation
  valid_genes <- res$gene[!is.na(res$cor)]
  gene <- valid_genes[1]
  
  # Create plot
  p <- plotExprAUC(
    pset = pset, 
    corResults = res, 
    gene = gene, 
    drug = drug,
    mDataType = "rna",
    sensitivity.measure = "aac_recomputed"
  )
  
  expect_s3_class(p, "ggplot")
  
  # Check plot has required components
  expect_true(!is.null(p$layers))
  expect_true(length(p$layers) > 0)
})

# Test volcanoAUC
test_that("volcanoAUC returns a ggplot object", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  # Get correlation results
  res <- correlateExpressionAUC(
    pset, 
    drug = drug,
    sensitivity.measure = "aac_recomputed"
  )
  
  # Create volcano plot
  p <- volcanoAUC(
    corResults = res,
    corThreshold = 0.3,
    pvalThreshold = 0.05
  )
  
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$layers))
})

# Test defineResponseGroups
test_that("defineResponseGroups returns a factor with correct levels", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  # Test median method
  groups <- defineResponseGroups(
    pset, 
    drug = drug,
    sensitivity.measure = "aac_recomputed", 
    method = "median"
  )
  
  # Check it's a factor
  expect_s3_class(groups, "factor")
  
  # Check levels
  expect_true(all(c("sensitive", "resistant") %in% levels(groups)))
  
  # Check it's named
  expect_true(!is.null(names(groups)))
  
  # Check length (should have entries for samples with non-NA values)
  expect_gt(length(groups), 0)
  
  # Check both groups are present
  expect_true("sensitive" %in% groups)
  expect_true("resistant" %in% groups)
})

# Test defineResponseGroups with quantile method
test_that("defineResponseGroups works with quantile method", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  groups <- defineResponseGroups(
    pset, 
    drug = drug,
    sensitivity.measure = "aac_recomputed", 
    method = "quantile",
    quantiles = c(0.25, 0.75)
  )
  
  expect_s3_class(groups, "factor")
  
  # With quantile method, some samples may be NA
  expect_true(any(!is.na(groups)))
})

# Test runDiffExpr
test_that("runDiffExpr returns a data.frame with expected columns", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  # Define groups first
  groups <- defineResponseGroups(
    pset, 
    drug = drug,
    sensitivity.measure = "aac_recomputed", 
    method = "median"
  )
  
  # Run differential expression
  diff_res <- runDiffExpr(
    pset, 
    groupLabels = groups, 
    mDataType = "rna", 
    method = "t.test"
  )
  
  # Check structure
  expect_s3_class(diff_res, "data.frame")
  
  # Check columns (match your actual output!)
  expect_true(all(c("gene", "logFC", "pval", "adj_pval") %in% colnames(diff_res)))
  
  # Check dimensions
  expect_gt(nrow(diff_res), 0)
  
  # Check p-values are in valid range (where not NA)
  valid_pvals <- diff_res$pval[!is.na(diff_res$pval)]
  expect_true(all(valid_pvals >= 0 & valid_pvals <= 1))
})

# Test plotGeneBoxplot
test_that("plotGeneBoxplot returns a ggplot object", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  # Define groups
  groups <- defineResponseGroups(
    pset, 
    drug = drug,
    sensitivity.measure = "aac_recomputed", 
    method = "median"
  )
  
  # Get a valid gene name
  expr <- summarizeMolecularProfiles(pset, mDataType = "rna", summary.stat = "mean")
  if (inherits(expr, "SummarizedExperiment")) {
    expr <- SummarizedExperiment::assay(expr)
  }
  gene <- rownames(expr)[1]
  
  # Create boxplot
  p <- plotGeneBoxplot(
    pset, 
    gene = gene,  # Parameter name is 'gene' not 'gene_of_interest'
    groupLabels = groups,  # Parameter name is 'groupLabels' not 'response_groups'
    mDataType = "rna", 
    plotType = "boxplot"
  )
  
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$layers))
})

# Test plotGeneBoxplot with violin plot
test_that("plotGeneBoxplot works with violin plot type", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  drug <- drugNames(pset)[1]
  
  groups <- defineResponseGroups(
    pset, 
    drug = drug,
    sensitivity.measure = "aac_recomputed", 
    method = "median"
  )
  
  expr <- summarizeMolecularProfiles(pset, mDataType = "rna", summary.stat = "mean")
  if (inherits(expr, "SummarizedExperiment")) {
    expr <- SummarizedExperiment::assay(expr)
  }
  gene <- rownames(expr)[1]
  
  p <- plotGeneBoxplot(
    pset, 
    gene = gene,
    groupLabels = groups,
    mDataType = "rna", 
    plotType = "violin"
  )
  
  expect_s3_class(p, "ggplot")
})

# Test plotDoseResponse
test_that("plotDoseResponse returns a ggplot object", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  
  # Get valid cell lines
  valid_cell_lines <- rownames(pset@treatmentResponse$raw[, , "Viability"])
  
  # Select first 2-3 cell lines
  cell_lines <- valid_cell_lines[1:min(3, length(valid_cell_lines))]
  
  # Create dose-response plot
  # Note: plotDoseResponse doesn't take 'drug' or 'sensitivity.measure' parameters
  # based on your function definition
  p <- plotDoseResponse(
    pset, 
    cell.lines = cell_lines,
    sensitivity.measure = "Viability"
  )
  
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(p$layers))
})

# Test error handling - invalid drug name
test_that("correlateExpressionAUC handles invalid drug name", {
  skip_if_no_pharmacogx()
  
  pset <- setup_test_pset()
  
  expect_error(
    correlateExpressionAUC(
      pset, 
      drug = "NonExistentDrug123",
      sensitivity.measure = "aac_recomputed"
    ),
    "not found"
  )
})

# Test error handling - invalid PSet
test_that("correlateExpressionAUC handles invalid PSet", {
  expect_error(
    correlateExpressionAUC(
      pset = list(),  # Not a valid PSet
      drug = "TestDrug",
      sensitivity.measure = "aac_recomputed"
    ),
    "PharmacoSet"
  )
})


