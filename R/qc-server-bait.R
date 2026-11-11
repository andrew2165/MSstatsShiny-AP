is_tmt_protein <- function(loadpage_input) {
  identical(loadpage_input$BIO, "Protein") && identical(loadpage_input$DDA_DIA, "TMT")
}

tmt_bait_settings <- function(input) {
  list(enabled = isTRUE(input$bait_enabled), protein = input$bait_protein,
       policy = if (identical(input$bait_missing, "drop")) "drop" else "error")
}

bait_panel_ui <- function(ns) {
  conditionalPanel(
    "input['loadpage-BIO'] === 'Protein' && input['loadpage-DDA_DIA'] === 'TMT'",
    wellPanel(
      h4("Bait normalization and QC"),
      p("After summarization, select the exact bait protein ID to inspect its log2 abundance."),
      selectizeInput(ns("bait_protein"), "Bait protein", choices = c("")),
      p("With normalization enabled, comparisons use client log2 abundance minus bait log2 abundance within each channel."),
      shinyjs::disabled(
        selectInput(ns("bait_missing"), "When a biological channel has no usable bait",
                    choices = c("Stop the comparison" = "error",
                                "Drop affected channels" = "drop"), selected = "error",
                    selectize = FALSE)
      ),
      conditionalPanel("input.bait_protein", ns = ns,
                       plotlyOutput(ns("bait_qc_plot"), height = "340px")),
      uiOutput(ns("bait_status"))
    )
  )
}

bait_qc_data <- function(proteins, bait, channels) {
  keys <- c("Mixture", "TechRepMixture", "Run", "Channel")
  sample_columns <- c(keys, "Condition", "BioReplicate")
  proteins <- as.data.frame(proteins)
  samples <- unique(rbind(proteins[sample_columns], channels))
  samples <- samples[!samples$Condition %in% c("Norm", "Empty"), , drop = FALSE]
  bait_rows <- proteins[which(proteins$Protein == bait), c(keys, "Abundance")]
  bait_rows$Status <- rep("Measured", nrow(bait_rows))
  qc <- merge(samples, bait_rows, by = keys, all.x = TRUE, sort = FALSE)
  qc$Status <- ifelse(is.na(qc$Status), "Absent",
                      ifelse(is.finite(qc$Abundance), "Measured", "Non-finite"))
  qc
}

bait_qc_plot <- function(qc) {
  groups <- sort(unique(as.character(qc$Condition)))
  # Give each channel a distinct horizontal position even for identical abundances.
  qc$position <- match(qc$Condition, groups) +
    stats::ave(seq_len(nrow(qc)), qc$Condition, FUN = function(i) {
      if (length(i) == 1L) 0 else seq(-0.32, 0.32, length.out = length(i))
    })
  qc$hover <- htmltools::htmlEscape(with(qc, paste0(
    Condition, " | Biological replicate ", BioReplicate, " | Run ", Run,
    " | Channel ", Channel, " | Mixture ", Mixture, " | Technical replicate ",
    TechRepMixture, " | ", Status,
    ifelse(is.finite(Abundance), paste0(" | Log2 abundance ", signif(Abundance, 5)), "")
  )))
  axis <- list(title = "Condition", tickvals = seq_along(groups), ticktext = groups,
               range = c(0.5, length(groups) + 0.5))
  measured <- qc[qc$Status == "Measured", , drop = FALSE]
  missing <- qc[qc$Status != "Measured", , drop = FALSE]
  abundance <- plot_ly(measured, x = ~position, y = ~Abundance, text = ~hover,
                       type = "scatter", mode = "markers", hoverinfo = "text",
                       showlegend = FALSE, marker = list(color = "#24678d"))
  abundance <- layout(abundance, xaxis = axis,
                      yaxis = list(title = "Bait log2 abundance"))
  if (!nrow(missing)) return(abundance)

  flags <- plot_ly(missing, x = ~position, y = ~Status, text = ~hover,
                   type = "scatter", mode = "markers", hoverinfo = "text",
                   showlegend = FALSE, marker = list(symbol = "x", color = "#b53322"))
  flags <- layout(flags, xaxis = axis, yaxis = list(title = "Bait QC"))
  if (!nrow(measured)) return(flags)
  plotly::subplot(abundance, flags, nrows = 2, heights = c(0.7, 0.3),
                  shareX = TRUE, titleX = TRUE, titleY = TRUE)
}

bait_adjust_summary <- function(summary, bait, policy) {
  warnings <- character()
  adjusted <- tryCatch(
    withCallingHandlers({
      if (is.null(bait) || !nzchar(bait)) {
        stop("Select a bait protein before enabling bait normalization.", call. = FALSE)
      }
      if (!"normalizeToBaitTMT" %in% getNamespaceExports("MSstatsTMT")) {
        stop("Bait normalization requires an MSstatsTMT build that exports normalizeToBaitTMT().",
             call. = FALSE)
      }
      MSstatsTMT::normalizeToBaitTMT(summary, bait = bait, missing_bait = policy)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = identity
  )
  if (inherits(adjusted, "error")) return(list(error = conditionMessage(adjusted)))
  list(data = adjusted, warnings = warnings)
}

bait_remaining_coverage <- function(proteins, conditions) {
  retained <- as.data.frame(proteins)
  retained <- retained[is.finite(retained$Abundance), , drop = FALSE]
  keys <- c("Mixture", "TechRepMixture", "Run", "Channel")
  counts <- lapply(sort(unique(as.character(conditions))), function(condition) {
    rows <- retained[retained$Condition == condition, , drop = FALSE]
    data.frame(Condition = condition,
               Channels = nrow(unique(rows[keys])),
               `Biological replicates` = length(unique(rows$BioReplicate)), check.names = FALSE)
  })
  do.call(rbind, counts)
}

bait_panel_server <- function(input, output, session, loadpage_input, summarization) {
  observe({
    shinyjs::toggleState("bait_missing", condition = isTRUE(input$bait_enabled))
  })
  choices <- reactive({
    req(is_tmt_protein(loadpage_input()))
    sort(unique(as.character(summarization()$data$ProteinLevelData$Protein)))
  })
  observeEvent(choices(), {
    selected <- intersect(input$bait_protein, choices())
    updateSelectizeInput(session, "bait_protein", choices = c("", choices()),
                         selected = if (length(selected)) selected else "", server = TRUE)
  })
  observeEvent(loadpage_input()$proceed1, {
    shiny::updateCheckboxInput(session, "bait_enabled", value = FALSE)
    updateSelectInput(session, "bait_missing", selected = "error")
    updateSelectizeInput(session, "bait_protein", selected = "")
  }, ignoreInit = TRUE)

  qc <- reactive({
    req(input$bait_protein %in% choices())
    result <- summarization()
    bait_qc_data(result$data$ProteinLevelData, input$bait_protein, result$channels)
  })
  output$bait_qc_plot <- renderPlotly({
    req(nrow(qc()) > 0L)
    bait_qc_plot(qc())
  })

  adjustment <- reactive({
    if (!is_tmt_protein(loadpage_input())) return(NULL)
    settings <- tmt_bait_settings(input)
    if (!settings$enabled) return(NULL)
    bait_adjust_summary(summarization()$data, settings$protein, settings$policy)
  })
  dropped <- reactive({
    attr(adjustment()$data, "BaitNormalization", exact = TRUE)$dropped_samples
  })
  output$bait_status <- renderUI({
    result <- adjustment()
    if (is.null(result)) return(NULL)
    if (!is.null(result$error)) {
      return(div(class = "alert alert-danger", role = "alert", result$error))
    }
    if (NROW(dropped())) {
      tagList(
        div(class = "alert alert-warning", role = "alert", result$warnings),
        h5("Dropped channels"),
        div(style = "overflow-x: auto", tableOutput(session$ns("bait_dropped_channels"))),
        h5("Remaining coverage for group comparison"),
        tableOutput(session$ns("bait_coverage"))
      )
    }
  })
  output$bait_dropped_channels <- renderTable({
    audit <- dropped()
    req(NROW(audit) > 0L)
    first <- c("Condition", "BioReplicate", "Run", "Channel", "Reason", "ClientRows")
    audit[c(first, setdiff(names(audit), first))]
  })
  output$bait_coverage <- renderTable({
    req(NROW(dropped()) > 0L)
    bait_remaining_coverage(adjustment()$data$ProteinLevelData, qc()$Condition)
  })
  adjustment
}
