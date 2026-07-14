################################################################################
# INDONESIA INTEGRATED NCD MODEL — PACKAGE DEPENDENCIES
# R/packages.R
# ─────────────────────────────────────────────────────────────────────────────
# Checks that all required packages are installed and stops with a clear
# install instruction if any are missing. Does NOT install packages silently.
#
# Usage: source(here::here("R", "packages.R")) at the top of any run script.
#
# To install all dependencies at once:
#   install.packages(c("here","dplyr","tidyr","readr","purrr","tibble",
#                      "ggplot2","scales","abind","stringr","readxl",
#                      "cowplot","data.table","countrycode","rlang"))
#   devtools::install_github("PPgp/wpp2024")
#
# For reproducible environments use renv:
#   renv::init()   — creates renv.lock from the current library
#   renv::restore()— restores the locked library on a new machine
################################################################################

if (getRversion() < "4.1.0") {
  stop(
    "This pipeline requires R >= 4.1.0 because it uses the native pipe (|>) and anonymous functions.",
    call. = FALSE
  )
}

REQUIRED_PKGS <- c(
  "here", "dplyr", "tidyr", "readr", "purrr", "tibble",
  "ggplot2", "scales", "abind", "stringr", "readxl",
  "cowplot", "data.table", "countrycode", "rlang"
)

missing_pkgs <- REQUIRED_PKGS[!vapply(REQUIRED_PKGS, requireNamespace,
                                       quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  stop(
    "The following packages are required but not installed:\n  ",
    paste(missing_pkgs, collapse = ", "), "\n\n",
    "Install them with:\n",
    "  install.packages(c(",
    paste0('"', missing_pkgs, '"', collapse = ", "),
    "))\n\n",
    "For the wpp2024 demographic package:\n",
    "  devtools::install_github('PPgp/wpp2024')\n",
    call. = FALSE
  )
}

# wpp2024 is only needed for script 02. Check separately so other scripts
# can still run if wpp2024 is absent.
if (!requireNamespace("wpp2024", quietly = TRUE)) {
  message("Note: wpp2024 not installed — required by scripts/02_build_demography.R only.")
  message("  Install: devtools::install_github('PPgp/wpp2024')")
}
