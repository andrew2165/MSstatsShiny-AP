test_that("bait QC distinguishes absent rows from non-finite abundances", {
  channels <- data.frame(
    Mixture = "Plex", TechRepMixture = 1, Run = "Run01",
    Channel = c("126", "127", "128", "129", "130"),
    Condition = c("Control", "Control", "Treatment", "Treatment", "Norm"),
    BioReplicate = c("C1", "C2", "T1", "T2", "Reference")
  )
  proteins <- rbind(
    transform(channels[c(1, 2, 5), ], Protein = "sp|001|BAIT", Abundance = c(12, Inf, 11)),
    transform(channels[1:3, ], Protein = "CLIENT", Abundance = 10)
  )

  qc <- bait_qc_data(proteins, "sp|001|BAIT", channels)
  expect_equal(qc$Status[match(channels$Channel[1:4], qc$Channel)],
               c("Measured", "Non-finite", "Absent", "Absent"))
  expect_equal(nrow(qc), 4L)
  expect_equal(qc$Abundance[qc$Channel == "126"], 12)
  traces <- plotly::plotly_build(bait_qc_plot(qc))$x$data
  expect_equal(vapply(traces, function(x) length(x$x), integer(1)), c(1L, 3L))
})

test_that("bait adjustment is opt-in and reports the fork's missing-channel audit", {
  expect_identical(tmt_bait_settings(list())$enabled, FALSE)
  expect_identical(tmt_bait_settings(list())$policy, "error")
  skip_if_not("normalizeToBaitTMT" %in% getNamespaceExports("MSstatsTMTAP"))

  channels <- data.frame(
    Mixture = "Plex", TechRepMixture = 1, Run = "Run01",
    Channel = c("126", "127"), Condition = c("Control", "Treatment"),
    BioReplicate = c("C1", "T1")
  )
  proteins <- rbind(
    transform(channels[1, ], Protein = "BAIT", Abundance = 12),
    transform(channels, Protein = "CLIENT", Abundance = c(10, 14))
  )
  original <- list(ProteinLevelData = proteins, FeatureLevelData = data.frame())
  err <- bait_adjust_summary(original, "BAIT", "error")
  expect_match(err$error, "absent")
  dropped <- bait_adjust_summary(original, "BAIT", "drop")
  expect_match(dropped$warnings, "Dropped 1 client row", fixed = TRUE)
  expect_equal(dropped$data$ProteinLevelData$Abundance, -2)
  expect_equal(attr(dropped$data, "BaitNormalization")$dropped_samples$Channel, "127")
  expect_identical(original$ProteinLevelData, proteins)
})

test_that("TMT protein annotation requires a mixture ID for each row", {
  input <- list(BIO = "Protein", DDA_DIA = "TMT")
  for (mixtures in list(c("PlexA", "PlexB"), c(1L, 2L))) {
    annotation <- data.frame(Mixture = mixtures)
    expect_identical(check_tmt_annotation_mixture(annotation, input), annotation)
  }
  for (mixture in c(NA_character_, "", "  ")) {
    expect_error(check_tmt_annotation_mixture(data.frame(Mixture = mixture), input),
                 "Mixture column in the sample annotation CSV")
  }
  expect_error(check_tmt_annotation_mixture(data.frame(Run = "Run01"), input),
               "sample annotation CSV needs a Mixture column")
})

test_that("FragPipe upload reports a blank mixture before conversion", {
  report <- tempfile(fileext = ".csv")
  annotation <- tempfile(fileext = ".csv")
  on.exit(unlink(c(report, annotation)))
  writeLines(c("Spectrum.File", "Run01.mzML"), report)
  writeLines(c("Run,Mixture,Channel", "Run01,,126"), annotation)
  input <- list(BIO = "Protein", DDA_DIA = "TMT", filetype = "phil", sep_data = ",",
                data = list(datapath = report), annot = list(datapath = annotation))

  mockery::stub(getData, "show_modal_spinner", NULL)
  expect_error(getData(input), "Mixture column in the sample annotation CSV")
})

test_that("TMT export quotes the exact bait and keeps the original summary", {
  bait <- 'sp|001|BAIT "identifier"'
  settings <- list(summarization = "Median", global_norm = FALSE, reference_norm = TRUE)
  upload <- list(BIO = "Protein", DDA_DIA = "TMT", filetype = "phil", sep_data = ",")
  contrast <- matrix(c(-1, 1), nrow = 1,
                     dimnames = list("Treatment minus Control", c("Control", "Treatment")))
  code <- generate_tmt_protein_code(settings, upload, contrast, list(),
                                    list(enabled = TRUE, protein = bait, policy = "drop"))
  expect_no_error(parse(text = code))
  expect_match(code, paste("bait_protein =", deparse1(bait)), fixed = TRUE)
  expect_match(code, 'missing_bait = "drop"', fixed = TRUE)
  expect_match(code, "comparison_data = MSstatsTMTAP::normalizeToBaitTMT(", fixed = TRUE)
  expect_match(code, "comparison_data, contrast.matrix = contrast.matrix", fixed = TRUE)
})

test_that("bait comparison rejects a contrast that uses an entirely lost condition", {
  adjusted <- list(ProteinLevelData = data.frame(
    Condition = c("Control", "Treatment"), Abundance = c(0, 1)
  ))
  contrast <- matrix(c(-1, 1, 0), nrow = 1,
                     dimnames = list("Treatment minus Control", c("Control", "Treatment", "Lost")))
  expect_identical(bait_contrast_matrix(adjusted, contrast), contrast[, 1:2, drop = FALSE])
  contrast[1, "Lost"] <- 1
  expect_error(bait_contrast_matrix(adjusted, contrast), "contrast conditions: Lost")
})

test_that("coverage includes lost conditions and counts distinct biological replicates", {
  proteins <- data.frame(
    Mixture = "Plex", TechRepMixture = c(1, 1, 2, 2), Run = c("R1", "R1", "R2", "R2"),
    Channel = c("126", "126", "126", "127"), BioReplicate = c("C1", "C1", "C1", "C2"),
    Condition = factor("Control"), Abundance = c(0, 1, 0, NA)
  )
  coverage <- bait_remaining_coverage(proteins, factor(c("Treatment", "Control")))
  expect_identical(coverage$Condition, c("Control", "Treatment"))
  expect_identical(coverage$Channels, c(2L, 0L))
  expect_identical(coverage[["Biological replicates"]], c(1L, 0L))
})
