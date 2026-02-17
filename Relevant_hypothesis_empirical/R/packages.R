# Minimal package bootstrapper

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

ensure_packages(c(
  "data.table",
  "stringi",
  "httr",
  "jsonlite",
  "lubridate",
  "ggplot2"
))
