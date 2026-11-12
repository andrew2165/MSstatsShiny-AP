.onLoad <- function(libname, pkgname) {
  shiny::addResourcePath(
    prefix = "msstats-ap-assets",
    directoryPath = system.file(
      "assets",
      package = "MSstatsShinyAP"
    )
  )
}

.onUnload <- function(libname, pkgname) {
  shiny::removeResourcePath("msstats-ap-assets")
}