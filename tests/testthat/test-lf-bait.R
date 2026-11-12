lf_bait_fixture <- function() {
  samples <- data.frame(RUN = 1:4, originalRUN = paste0("run", 1:4),
                        GROUP = rep(c("Control", "Treatment"), each = 2),
                        SUBJECT = paste0("bio", 1:4), LABEL = "L")
  list(ProteinLevelData = rbind(
    transform(samples, Protein = "bait|id", LogIntensities = c(20, 22, 21, 24), Variance = .09),
    transform(samples, Protein = "client", LogIntensities = c(18, 20, 20, 23), Variance = .04)),
    FeatureLevelData = samples)
}

test_that("label-free QC includes absent, non-finite and fully lost sample summaries", {
  x <- lf_bait_fixture()
  p <- x$ProteinLevelData[-2, ]
  p$LogIntensities[p$Protein == "bait|id" & p$RUN == 3] <- Inf
  p <- p[p$RUN != 4, ]
  qc <- bait_qc_data(p, "bait|id", x$FeatureLevelData)
  expect_equal(qc$Status[match(1:4, qc$RUN)], c("Measured", "Absent", "Non-finite", "Absent"))
  traces <- plotly::plotly_build(bait_qc_plot(qc))$x$data
  expect_equal(sum(vapply(traces, function(x) length(x$x), integer(1))), 4L)
  expect_match(traces[[1]]$text, "Biological replicate bio1 | Run run1", fixed = TRUE)
})

test_that("label-free settings restrict the workflow and default to unweighted opt-in", {
  expect_true(is_bait_protein(list(BIO = "Protein", DDA_DIA = "LType")))
  expect_false(is_bait_protein(list(BIO = "PTM", DDA_DIA = "LType")))
  expect_false(is_bait_protein(list(BIO = "Peptide", DDA_DIA = "LType")))
  expect_false(tmt_bait_settings(list())$enabled)
  expect_identical(tmt_bait_settings(list())$variance_method, "ignore")
})

test_that("label-free adjustment delegates variance, missing samples and audit to MSstats", {
  skip_if_not("normalizeToBait" %in% getNamespaceExports("MSstats"))
  x <- lf_bait_fixture()
  adjusted <- bait_adjust_summary(x, "bait|id", "error", label_free = TRUE)
  expect_null(adjusted$error)
  expect_equal(adjusted$data$ProteinLevelData$LogIntensities, c(-2, -2, -1, -1))
  expect_true(all(is.na(adjusted$data$ProteinLevelData$Variance)))
  propagated <- bait_adjust_summary(x, "bait|id", "error", "propagate", TRUE)
  expect_equal(propagated$data$ProteinLevelData$Variance, rep(.13, 4))
  x$ProteinLevelData <- x$ProteinLevelData[-1, ]
  expect_match(bait_adjust_summary(x, "bait|id", "error", label_free = TRUE)$error, "absent")
  dropped <- bait_adjust_summary(x, "bait|id", "drop", label_free = TRUE)
  expect_match(dropped$warnings, "Dropped 1 client row", fixed = TRUE)
  expect_equal(attr(dropped$data, "BaitNormalization")$dropped_samples$originalRUN, "run1")
  coverage <- bait_remaining_coverage(dropped$data$ProteinLevelData, c("Control", "Treatment", "Lost"))
  expect_equal(coverage$Samples, c(1L, 0L, 2L))
  expect_equal(x$ProteinLevelData$LogIntensities[x$ProteinLevelData$Protein == "client"], c(18, 20, 20, 23))
})

test_that("label-free contrasts guard conditions lost by dropping samples", {
  x <- lf_bait_fixture()
  x$ProteinLevelData <- x$ProteinLevelData[x$ProteinLevelData$GROUP == "Control", ]
  contrast <- matrix(c(-1, 1), 1, dimnames = list("T-C", c("Control", "Treatment")))
  expect_error(bait_contrast_matrix(x, contrast), "contrast conditions: Treatment")
})

test_that("label-free exports retain exact bait, variance policy and original summary", {
  contrast <- matrix(c(-1, 1), 1, dimnames = list("T-C", c("Control", "Treatment")))
  bait <- list(enabled = TRUE, protein = 'bait|"id"', policy = "drop", variance_method = "propagate")
  code <- generate_analysis_code(list(log = "2", norm = "FALSE"),
    list(BIO = "Protein", DDA_DIA = "LType", filetype = "phil"), contrast, list(), bait = bait)
  expect_no_error(parse(text = code))
  expect_match(code, paste("bait_protein =", deparse1(bait$protein)), fixed = TRUE)
  expect_match(code, 'variance_method = "propagate"', fixed = TRUE)
  expect_match(code, "comparison_data = MSstats::normalizeToBait(", fixed = TRUE)
  expect_match(code, "data = comparison_data", fixed = TRUE)
  bait$enabled <- FALSE
  off <- generate_lf_protein_code(list(), list(BIO = "Protein", DDA_DIA = "LType", filetype = "phil"), contrast, bait)
  expect_match(off, "comparison_data = summarized", fixed = TRUE)
  expect_false(grepl("MSstats::normalizeToBait(", off, fixed = TRUE))
})

test_that("missing fork functions and unsupported label-free settings fail clearly", {
  mockery::stub(bait_adjust_summary, "getNamespaceExports", character())
  expect_match(bait_adjust_summary(lf_bait_fixture(), "bait|id", "error", label_free = TRUE)$error,
               "MSstats build that exports normalizeToBait()", fixed = TRUE)
})

test_that("label-free normalization preserves log2 and variance safeguards", {
  skip_if_not("normalizeToBait" %in% getNamespaceExports("MSstats"))
  x <- lf_bait_fixture()
  attr(x$ProteinLevelData, "MSstatsAbundanceScale") <- "log10"
  expect_match(bait_adjust_summary(x, "bait|id", "error", label_free = TRUE)$error, "log2 input")
  attr(x$ProteinLevelData, "MSstatsAbundanceScale") <- "log2"
  x$ProteinLevelData$Variance <- NA_real_
  expect_match(bait_adjust_summary(x, "bait|id", "error", "propagate", TRUE)$error,
               "Invalid or unavailable")
  expect_null(bait_adjust_summary(x, "bait|id", "error", "ignore", TRUE)$error)
})
