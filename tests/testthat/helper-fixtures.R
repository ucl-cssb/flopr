# Shared test fixture helpers.
#
# Several parser functions write their output next to their input file, and
# some (until backlog 001/002/003 were fixed) even overwrote the input file
# itself - tests must never run against tests/testthat/fixtures/ in place,
# or a bug regression could silently corrupt the repo's own committed
# fixtures on every test run. local_fixture_copy() copies the named module's
# fixtures into a fresh temp directory and cleans it up when the calling
# test exits.

#' Copy a fixture module's files into a fresh temp directory
#'
#' @param module subdirectory name under tests/testthat/fixtures/
#' @param env environment whose exit triggers cleanup (defaults to the
#'   caller, i.e. the running test)
#'
#' @return path to the temp directory containing the copied fixtures
local_fixture_copy <- function(module, env = parent.frame()) {
  src_dir <- testthat::test_path("fixtures", module)
  dest_dir <- withr::local_tempdir(pattern = paste0("flopr-", module, "-"), .local_envir = env)
  files <- list.files(src_dir, full.names = TRUE)
  file.copy(files, dest_dir)
  dest_dir
}

#' A compact, order-insensitive numeric fingerprint of a data.frame
#'
#' Full data.frame snapshots get large and can be noisy across platforms
#' (row order, floating point representation). This instead snapshots shape
#' (rows/columns), column names, and per-numeric-column rounded sums/NA
#' counts - sensitive to real regressions in the computed values, stable
#' across row reordering.
#'
#' @param df a data.frame
#' @param digits rounding applied to numeric sums before snapshotting
#'
#' @return a list suitable for testthat::expect_snapshot_value()
fingerprint <- function(df, digits = 3) {
  is_num <- vapply(df, is.numeric, logical(1))
  list(
    dim = dim(df),
    names = names(df),
    numeric_sums = round(vapply(df[is_num], sum, numeric(1), na.rm = TRUE), digits),
    numeric_na_counts = vapply(df[is_num], function(x) sum(is.na(x)), integer(1))
  )
}
