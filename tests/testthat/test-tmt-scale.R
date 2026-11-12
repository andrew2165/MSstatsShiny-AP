tmt_scale_settings <- function(log = "2") {
  list(summarization = "msstats", global_norm = FALSE, reference_norm = FALSE,
       remove_norm_channel = FALSE, null = TRUE, null1 = TRUE, norm = "FALSE",
       log = log, features_used = "all", n_feat = 3, summaryMethod = "TMP",
       MBi = FALSE, censInt = "NA", remove50 = FALSE)
}

test_that("TMT summarization ignores hidden label-free log base without mutating it", {
  data <- as.data.frame(MSstatsTMTAP::input.pd)
  data <- data[data$Run == data$Run[1], ]
  settings <- tmt_scale_settings("10")
  before <- settings
  expected <- tmt_summarization_loop(data, tmt_scale_settings(), list(DDA_DIA = "TMT"), FALSE)
  actual <- tmt_summarization_loop(data, settings, list(DDA_DIA = "TMT"), FALSE)
  expect_equal(actual, expected)
  expect_identical(settings, before)
  expect_identical(attr(actual$ProteinLevelData, "MSstatsTMTAbundanceScale"), "log2")
})

test_that("TMT exported settings record the effective log2 scale", {
  settings <- tmt_scale_settings("10")
  code <- generate_tmt_protein_code(settings,
    list(BIO = "Protein", DDA_DIA = "TMT", filetype = "phil", sep_data = ","),
    matrix(c(-1, 1), 1, dimnames = list("T-C", c("Control", "Treatment"))),
    list(), list(enabled = TRUE, protein = "bait", policy = "error"))
  statements <- as.list(parse(text = code))
  assignment <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("=")) &&
                         identical(x[[2]], as.name("qc_settings")), statements)[[1]]
  env <- new.env()
  eval(assignment, env)
  expect_identical(env$qc_settings$log, "2")
  expect_identical(settings$log, "10")
})

test_that("TMT bait QC precedes normalization and rejects unverified summaries", {
  channels <- data.frame(Mixture = "Plex", TechRepMixture = 1, Run = "R1",
    Channel = c("126", "127"), Condition = c("Control", "Treatment"),
    BioReplicate = c("C1", "T1"))
  original <- list(FeatureLevelData = data.frame(), ProteinLevelData = rbind(
    transform(channels, Protein = "bait", Abundance = c(20, 22)),
    transform(channels, Protein = "client", Abundance = c(18, 21))))
  attr(original$ProteinLevelData, "MSstatsTMTAbundanceScale") <- "log2"
  shiny::testServer(function(input, output, session) {
    summary <- shiny::reactiveVal(list(data = original,
      upload = list(BIO = "Protein", DDA_DIA = "TMT", proceed1 = 1L),
      qc = list(log = "2"), channels = NULL))
    adjustment <- bait_panel_server(input, output, session,
      shiny::reactive(summary()$upload), summary)
  }, {
    session$setInputs(bait_enabled = FALSE, bait_protein = "bait", bait_missing = "error")
    expect_null(adjustment())
    expect_true(nchar(output$bait_qc_plot) > 0)
    session$setInputs(bait_enabled = TRUE)
    expect_null(adjustment()$error)
    expect_equal(adjustment()$data$ProteinLevelData$Abundance, c(-2, -1))
    for (scale in list(NULL, "log10", "log2_client_bait")) {
      invalid <- summary()
      attr(invalid$data$ProteinLevelData, "MSstatsTMTAbundanceScale") <- scale
      summary(invalid)
      session$flushReact()
      expect_match(adjustment()$error, "Rerun protein summarization")
      expect_error(output$bait_qc_plot, "Rerun protein summarization")
    }
  })
})

test_that("switching from label-free to TMT snapshots log2 without resetting the UI", {
  server <- qcServer
  original <- MSstatsTMTAP::quant.pd.msstats
  attr(original$ProteinLevelData, "MSstatsTMTAbundanceScale") <- "log2"
  mockery::stub(server, "preprocessData", original)
  shiny::testServer(function(input, output, session) {
    upload <- shiny::reactiveValues(BIO = "Protein", DDA_DIA = "LType", proceed1 = 1L)
    qc <- server(input, output, session, session, shiny::reactive(upload),
      shiny::reactive(original$ProteinLevelData),
      app_template = shiny::reactive(TEMPLATES$default))
  }, {
    do.call(session$setInputs, tmt_scale_settings("10"))
    session$setInputs(run = 1L)
    expect_identical(qc$tmtSummary()$qc$log, "10")
    upload$DDA_DIA <- "TMT"
    upload$proceed1 <- 2L
    session$setInputs(run = 2L)
    expect_identical(qc$tmtSummary()$qc$log, "2")
    expect_identical(input$log, "10")
  })
})
