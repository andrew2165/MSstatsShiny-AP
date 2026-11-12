test_that("exported readers preserve all user-selected separators", {
  cases <- list(
    c("msstats", "Protein", "LType", "sep_msstatsdata"),
    c("msstats", "Protein", "TMT", "sep_msstatsdata"),
    c("msstats", "PTM", "LType", "sep_msstatsptmdata"),
    c("10col", "Protein", "LType", "sep_data"),
    c("sky", "Protein", "LType", "sep_skylinedata"),
    c("prog", "Protein", "LType", "sep_data"),
    c("PD", "Protein", "LType", "sep_data"),
    c("PD", "Protein", "TMT", "sep_data"),
    c("spec", "Protein", "LType", "sep_specdata"),
    c("open", "Protein", "LType", "sep_data"),
    c("openms", "Protein", "LType", "sep_data"),
    c("openms", "Protein", "TMT", "sep_data"),
    c("spmin", "Protein", "TMT", "sep_data"),
    c("phil", "Protein", "TMT", "sep_data")
  )
  path <- tempfile()
  on.exit(unlink(path))
  expected <- data.frame(Protein = c("bait", "client"), Intensity = c(100L, 200L))
  for (case in cases) {
    for (separator in c(",", ";", "\t", "|")) {
      input <- list(filetype = case[1], BIO = case[2], DDA_DIA = case[3],
                    remove = TRUE, which.proteinid = "Protein")
      input[[case[4]]] <- separator
      code <- getDataCode(input)
      expressions <- parse(text = code)
      # Execute the generated reader itself, without invoking downstream converters.
      reads <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("=")) &&
                        is.call(x[[3]]) && is.symbol(x[[3]][[1]]) &&
                        as.character(x[[3]][[1]]) %in%
                        c("read.csv", "read.delim"), as.list(expressions))
      reader <- reads[[1]][[3]]
      expect_identical(reader$sep, separator, info = paste(case[1:3], collapse = "/"))
      write.table(expected, path, sep = separator, row.names = FALSE, quote = TRUE)
      reader[[2]] <- path
      expect_equal(eval(reader), expected)
    }
  }
})
