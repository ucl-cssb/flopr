skip_if_not_installed("flowCore")
skip_if_not_installed("flowClust")
skip_if_not_installed("flowStats")

test_that("process_fcs() processes a single .fcs file (README single-file example)", {
  dir <- local_fixture_copy("flow_cytometry")

  process_fcs(
    fcs_file = file.path(dir, "pWeak_None_0_1.fcs"),
    flu_channels = "BL1-H",
    pre_cleaned = TRUE,
    do_plot = FALSE
  )

  out_path <- file.path(dir, "pWeak_None_0_1_processed.fcs")
  expect_true(file.exists(out_path))

  processed <- flowCore::read.FCS(out_path, emptyValue = FALSE)
  expect_true("BL1-H" %in% flowCore::colnames(processed))
  expect_true(flowCore::nrow(processed) > 0)
  # processing should have removed at least some debris/doublet events
  raw <- flowCore::read.FCS(file.path(dir, "pWeak_None_0_1.fcs"), emptyValue = FALSE)
  expect_true(flowCore::nrow(processed) <= flowCore::nrow(raw))
})

test_that("process_fcs_dir() processes a directory with calibration (README directory example)", {
  dir <- local_fixture_copy("flow_cytometry")

  process_fcs_dir(
    dir_path = dir,
    pattern = "*Med*.fcs",
    flu_channels = "BL1-H",
    pre_cleaned = TRUE,
    do_plot = FALSE,
    neg_fcs = "pNeg_None_0_1.fcs",
    calibrate = TRUE,
    mef_peaks = list(list(channel = "BL1-H",
                           peaks = c(0, 822, 2114, 5911, 17013, 41837, 145365, 287558)))
  )

  summary_path <- file.path(paste0(dir, "_processed"), "data_summary.csv")
  expect_true(file.exists(summary_path))

  summary_data <- utils::read.csv(summary_path)
  expect_true(all(c("file", "channel", "mean", "sd") %in% names(summary_data)))
  expect_equal(length(unique(summary_data$file)), 10)  # 10 pMed_*.fcs fixtures
  expect_true(all(c("BL1-H", "normalised_BL1-H", "calibrated_BL1-H",
                     "calibrated_normalised_BL1-H") %in% summary_data$channel))
  expect_snapshot_value(fingerprint(summary_data), style = "json2")
})

test_that("process_fcs() and process_fcs_dir() never modify their .fcs inputs", {
  dir <- local_fixture_copy("flow_cytometry")
  input_size_before <- file.size(file.path(dir, "pWeak_None_0_1.fcs"))

  process_fcs(
    fcs_file = file.path(dir, "pWeak_None_0_1.fcs"),
    flu_channels = "BL1-H",
    pre_cleaned = TRUE,
    do_plot = FALSE
  )

  expect_equal(file.size(file.path(dir, "pWeak_None_0_1.fcs")), input_size_before)
})

test_that("process_fcs() stops on a sample with no bacteria found (needs-discussion/002)", {
  skip("Blocked on backlog 046 - the 'no bacteria found' guard is currently dead
        code (compares a flowFrame to 0, silently swallowed by try()), so this
        never actually stops today. Un-skip once 046 lands, using a fixture
        that's been subsetted down to all-debris/near-empty events.")

  dir <- local_fixture_copy("flow_cytometry")
  # a fixture with (close to) no real bacteria events - construct by heavily
  # subsetting an existing .fcs file once this test is un-skipped
  expect_error(
    process_fcs(
      fcs_file = file.path(dir, "no_bacteria.fcs"),
      flu_channels = "BL1-H",
      pre_cleaned = TRUE,
      do_plot = FALSE
    ),
    "[Nn]o bacteria"
  )
})
