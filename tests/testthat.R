#!/usr/bin/env Rscript
################################################################################
# tests/testthat.R
#
# Description: runs the unit tests in tests/testthat/ against the helpers in R/.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
# Pseudo code: source every R/*.R file, then testthat::test_dir(); exit 1 on
#              any failure.
# Inputs:  none (run from the repository root: Rscript tests/testthat.R)
# Outputs: test report on stdout
################################################################################
suppressPackageStartupMessages({
  library(testthat)
  library(data.table)
})
helper_env <- new.env()
for (helper in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  sys.source(helper, envir = helper_env)
}
results <- test_dir("tests/testthat", env = helper_env, reporter = "summary", stop_on_failure = TRUE)
