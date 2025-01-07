# Remove warning message from marginaleffects package
.onLoad <- function(libname, pkgname) {
  # Set default options for your package
  if (is.null(getOption("marginaleffects_safe"))) {
    options(marginaleffects_safe = FALSE)
  }
}