#' @title Run the SynergyLMM Shiny application
#' @description
#' `runSynergyLMMApp` launches the interactive 'SynergyLMM' web application locally. The app provides
#' a point-and-click front-end for the whole workflow: uploading tumor growth data, fitting the
#' mixed-effects model with [lmmModel()], evaluating synergy with [lmmSynergy()], inspecting model
#' diagnostics, and running the _post hoc_ and _a priori_ power analyses.
#' @param ... Additional arguments passed to [shiny::runApp()], for example `launch.browser` or
#' `port`.
#' @details
#' The application is shipped with the package, but the packages it needs to run are only listed in
#' `Suggests`, so that they are not a hard dependency of 'SynergyLMM'. If any of them is missing,
#' the function stops and names every missing package in a single message.
#'
#' The app is also available online at <https://synergylmm.uiocloud.no/>.
#' @returns This function is called for its side effect of starting the Shiny application and does
#' not return a meaningful value. It blocks the R session until the app is stopped.
#' @seealso [lmmModel()], [lmmSynergy()]
#' @examples
#' \dontrun{
#' runSynergyLMMApp()
#' }
#' @export

runSynergyLMMApp <- function(...) {
  pkgs <- c("shiny", "bslib", "shinyjs", "shinyWidgets", "shinyhelper",
            "DT", "plotly", "readxl", "openxlsx")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop("The SynergyLMM app requires these packages: ",
         paste(missing, collapse = ", "),
         "\nInstall them with install.packages(c(",
         paste0('"', missing, '"', collapse = ", "),
         "))",
         call. = FALSE)
  }

  app_dir <- system.file("shiny", package = "SynergyLMM")

  if (app_dir == "") {
    stop("Could not find the app directory. Try reinstalling SynergyLMM.", call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
