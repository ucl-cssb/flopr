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
