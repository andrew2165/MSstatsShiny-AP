generate_analysis_code = function(qc_input, loadpage_input, comp_mat, input, app_template = TEMPLATES$default, bait = NULL) {
  if (!is.null(bait)) {
    return(generate_tmt_protein_code(qc_input, loadpage_input, comp_mat, input, bait))
  }
  codes = preprocessDataCode(qc_input, loadpage_input)

  # Check if this is a response curve analysis
  is_response_curve = !is.null(input[[NAMESPACE_STATMODEL$comparison_mode]]) &&
    input[[NAMESPACE_STATMODEL$comparison_mode]] == CONSTANTS_STATMODEL$comparison_mode_response_curve

  if (is_response_curve) {
    increasing = input[[NAMESPACE_STATMODEL$modeling_response_curve_increasing_trend]]
    transform_dose = TRUE
    ratio_response = isTRUE(input[[NAMESPACE_STATMODEL$visualization_response_curve_ratio_scale]])

    codes = paste(codes, "\n# Set up dose response analysis\n", sep = "")
    codes = paste(codes, "library(MSstatsResponse)\n", sep = "")

    # Serialize the contrast matrix as a data frame
    codes = paste(codes, "group_metadata = data.frame(\n", sep = "")
    codes = paste(codes, "  GROUP = c(\"", paste(comp_mat$GROUP, collapse = "\",\""), "\"),\n", sep = "")
    for (col_name in setdiff(colnames(comp_mat), "GROUP")) {
      col_values = comp_mat[[col_name]]
      if (is.numeric(col_values)) {
        codes = paste(codes, "  ", col_name, " = c(", paste(col_values, collapse = ","), "),\n", sep = "")
      } else {
        codes = paste(codes, "  ", col_name, " = c(\"", paste(col_values, collapse = "\",\""), "\"),\n", sep = "")
      }
    }
    codes = paste(codes, "  stringsAsFactors = FALSE\n)\n", sep = "")

    # Merge and map columns explicitly (prepare_dose_response_fit is internal to MSstatsShiny)
    codes = paste(codes, "\n# Merge metadata with protein-level data\n", sep = "")
    codes = paste(codes, "protein_level_data = merge(summarized$ProteinLevelData, group_metadata, by = \"GROUP\")\n", sep = "")

    codes = paste(codes, "\n# Map columns to MSstatsResponse format\n", sep = "")
    codes = paste(codes, "protein_level_data$protein = protein_level_data$Protein\n", sep = "")
    codes = paste(codes, "protein_level_data$response = protein_level_data$LogIntensities\n", sep = "")

    # Determine the intervention column based on what's in the contrast matrix
    if ("drug" %in% colnames(comp_mat)) {
      codes = paste(codes, "protein_level_data$drug = protein_level_data$drug\n", sep = "")
      codes = paste(codes, "protein_level_data$dose = protein_level_data$dose_value\n", sep = "")
    } else {
      intervention_cols = grep("time|temperature|treatment", colnames(comp_mat),
                               ignore.case = TRUE, value = TRUE)
      if (length(intervention_cols) > 0) {
        intervention_type = sub("_.*", "", intervention_cols[1])
        value_col = paste0(intervention_type, "_value")
        codes = paste(codes, "protein_level_data$drug = \"", intervention_type, "\"\n", sep = "")
        codes = paste(codes, "protein_level_data$dose = protein_level_data$", value_col, "\n", sep = "")
      }
    }

    codes = paste(codes, "prepared_data = protein_level_data[, c(\"protein\", \"drug\", \"dose\", \"response\")]\n", sep = "")

    codes = paste(codes, "\n# Fit dose response curves\n", sep = "")
    codes = paste(codes, "response_results = doseResponseFit(\n",
                  "  data = prepared_data,\n",
                  "  increasing = ", increasing, ",\n",
                  "  transform_dose = ", transform_dose, ",\n",
                  "  ratio_response = FALSE\n)\n", sep = "")

    codes = paste(codes, "\n# Visualize response curves\n", sep = "")
    if (app_template == TEMPLATES$protein_turnover) {
      codes = paste(codes, "visualizeResponseProtein(\n",
                    "  data = prepared_data,\n",
                    "  protein_name = \"Enter protein name here\",\n",
                    "  drug_name = \"time\",\n",
                    "  ratio_response = FALSE,\n",
                    "  show_ic50 = TRUE,\n",
                    "  add_ci = FALSE,\n",
                    "  transform_dose = FALSE,\n",
                    "  n_samples = 1000,\n",
                    "  increasing = ", increasing, ",\n",
                    "  precalculated_ratios = TRUE,\n",
                    "  color_by = \"BaseSequence\",\n",
                    "  target_response = 0.5\n)\n", sep = "")
    } else {
      codes = paste(codes, "visualizeResponseProtein(\n",
                    "  data = prepared_data,\n",
                    "  protein_name = \"Enter protein name here\",\n",
                    "  drug_name = \"Enter drug name here\",\n",
                    "  ratio_response = ", ratio_response, ",\n",
                    "  show_ic50 = TRUE,\n",
                    "  add_ci = TRUE,\n",
                    "  transform_dose = ", transform_dose, ",\n",
                    "  n_samples = 1000,\n",
                    "  increasing = ", increasing, "\n)\n", sep = "")
    }

    return(codes)
  }

  # Standard analysis - build contrast matrix
  codes = paste(codes, "\n# Create the contrast matrix\n", sep = "")
  codes = paste(codes, "contrast.matrix = NULL\n", sep = "")
  
  for (i in seq_len(nrow(comp_mat))) {
    codes = paste(codes, "comparison = matrix(c(", 
                  toString(comp_mat[i,]), "),nrow=1)\n", sep = "")
    codes = paste(codes, "contrast.matrix = rbind(contrast.matrix, comparison)\n", sep = "")
  }
  
  codes = paste(codes, "row.names(contrast.matrix)=c(\"", 
                paste(row.names(comp_mat), collapse = '","'), "\")\n", sep = "")
  codes = paste(codes, "colnames(contrast.matrix)=c(\"", 
                paste(colnames(comp_mat), collapse = '","'), "\")\n", sep = "")
  
  if (loadpage_input$DDA_DIA == "TMT") {
    codes = paste(codes, "\n# Model-based comparison\n", sep = "")
    codes = paste(codes, "model = MSstatsTMT::groupComparisonTMT(summarized,
                       contrast.matrix = contrast.matrix,
                       moderated = ", input[[NAMESPACE_STATMODEL$modeling_tmt_moderation]], ",\t\t\t\t
                       adj.method = \"BH\",
                       remove_norm_channel = TRUE,
                       remove_empty_channel = TRUE
                       )\n", sep = "")
  } else if (loadpage_input$BIO == "PTM") {
    dt = if ((loadpage_input$BIO == "PTM" & loadpage_input$DDA_DIA == "TMT") | 
             loadpage_input$filetype == 'phil') "TMT" else "LabelFree"
    
    codes = paste(codes, "\n# Model-based comparison\n", sep = "")
    codes = paste(codes, "model = MSstatsPTM::groupComparisonPTM(summarized, '",
                  dt, "', \t\t\t\t
                      contrast.matrix = contrast.matrix)\n", sep = "")
  } else {
    codes = paste(codes, "\n# Model-based comparison\n", sep = "")
    codes = paste(codes, "model = MSstats::groupComparison(contrast.matrix, summarized)\n", sep = "")
  }
  
  if (loadpage_input$BIO == "PTM") {
    codes = paste(codes, "groupComparisonPlotsPTM(data=model,
                               type=\"Enter VolcanoPlot, Heatmap, or ComparisonPlot\",
                               which.Comparison=\"all\",
                               which.PTM=\"all\",
                               address=\"\")\n", sep = "")
  } else {
    codes = paste(codes, "groupComparisonPlots(data=model$ComparisonResult,
                               type=\"Enter VolcanoPlot, Heatmap, or ComparisonPlot\",
                               which.Comparison=\"all\",
                               which.Protein=\"all\",isPlotly=FALSE,
                               address=\"\")\n", sep = "")
  }
  
  return(codes)
}

generate_tmt_protein_code <- function(qc_input, loadpage_input, comp_mat, input, bait) {
  literal <- function(value) paste(deparse(value), collapse = "\n")
  qc_fields <- c("summarization", "global_norm", "reference_norm", "remove_norm_channel",
                 "null", "null1", "maxQC1", "norm", "log", "names", "features_used",
                 "n_feat", "summaryMethod", "MBi", "censInt", "remove50")
  moderation <- NAMESPACE_STATMODEL$modeling_tmt_moderation
  moderated <- isTRUE(as.logical(input[[moderation]]))
  sha <- utils::packageDescription("MSstatsTMT")$RemoteSha

  codes <- paste0(
    getDataCode(loadpage_input),
    "\n# Summarize with the same settings and function used in MSstatsShiny.\n",
    "qc_settings = ", literal(qc_input[intersect(qc_fields, names(qc_input))]), "\n",
    "summarized = MSstatsShiny::tmt_summarization_loop(\n",
    "  data, qc_settings, list(DDA_DIA = \"TMT\"), busy_indicator = FALSE\n)\n",
    "\n# Optional original protein intensity plot\n",
    "# MSstatsTMT::dataProcessPlotsTMT(summarized, type = \"QCPlot\")\n",
    "\ncontrast.matrix = ", literal(comp_mat), "\n",
    "\nbait_protein = ", literal(bait$protein), "\n",
    "missing_bait = ", literal(bait$policy), "\n"
  )
  if (isTRUE(bait$enabled)) {
    return(paste0(
      codes,
      if (length(sha)) paste0("# MSstatsTMT source SHA ", sha, "\n"),
      "comparison_data = MSstatsTMT::normalizeToBaitTMT(\n",
      "  summarized, bait = bait_protein, missing_bait = missing_bait\n)\n",
      "bait_audit = attr(comparison_data, \"BaitNormalization\")\n",
      "model = MSstatsTMT::groupComparisonTMT(\n",
      "  comparison_data, contrast.matrix = contrast.matrix, moderated = ", moderated, ",\n",
      "  save_fitted_models = TRUE, use_log_file = FALSE, verbose = FALSE\n)\n"
    ))
  }
  paste0(codes, "comparison_data = summarized\n",
         "model = MSstatsShiny::tmt_model(\n",
         "  comparison_data, ", literal(stats::setNames(list(moderated), moderation)),
         ", contrast.matrix, busy_indicator = FALSE\n)\n")
}
