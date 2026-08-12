# Both example fixtures are single-endpoint calibration reads, not kinetic
# timeseries data - they don't contain an "End Kinetic" marker, so despite
# cytation_parse()'s default being timeseries=TRUE, both need
# timeseries=FALSE explicitly (confirmed by inspecting the raw fixture
# content directly; matches examples/plate_reader/biotek_cytation/
# calibration_example.R, which already does this for the membrane file).

test_that("cytation_parse() parses lid calibration data (timeseries=FALSE)", {
  dir <- local_fixture_copy("biotek_cytation")

  out <- cytation_parse(
    data_file = file.path(dir, "231031_pr_calibration_lid.xlsx"),
    layout_csv = file.path(dir, "plate_layout.csv"),
    timeseries = FALSE
  )

  expect_true(file.exists(file.path(dir, "231031_pr_calibration_lid_parsed.csv")))
  expect_true(all(c("well", "row", "column") %in% names(out)))
  expect_snapshot_value(fingerprint(out), style = "json2")
})

test_that("cytation_parse() parses membrane calibration data (timeseries=FALSE)", {
  dir <- local_fixture_copy("biotek_cytation")

  out <- cytation_parse(
    data_file = file.path(dir, "231031_pr_calibration_membrane.xlsx"),
    layout_csv = file.path(dir, "plate_layout.csv"),
    timeseries = FALSE
  )

  expect_true(file.exists(file.path(dir, "231031_pr_calibration_membrane_parsed.csv")))
  expect_true(all(c("well", "row", "column") %in% names(out)))
  expect_snapshot_value(fingerprint(out), style = "json2")
})

test_that("cytation_parse() -> generate_cfs() produces conversion factors", {
  dir <- local_fixture_copy("biotek_cytation")

  cytation_parse(
    data_file = file.path(dir, "231031_pr_calibration_membrane.xlsx"),
    layout_csv = file.path(dir, "plate_layout.csv"),
    timeseries = FALSE
  )
  cfs <- generate_cfs(file.path(dir, "231031_pr_calibration_membrane_parsed.csv"))

  expect_true(all(c("cf", "beta", "calibrant", "fluorophore", "measure") %in% names(cfs)))
  expect_true(nrow(cfs) > 0)
  expect_snapshot_value(fingerprint(cfs), style = "json2")
})

test_that("cytation_parse() never overwrites its .xlsx input with parsed output (backlog 001)", {
  dir <- local_fixture_copy("biotek_cytation")
  input_size_before <- file.size(file.path(dir, "231031_pr_calibration_lid.xlsx"))

  cytation_parse(
    data_file = file.path(dir, "231031_pr_calibration_lid.xlsx"),
    layout_csv = file.path(dir, "plate_layout.csv"),
    timeseries = FALSE
  )

  input_size_after <- file.size(file.path(dir, "231031_pr_calibration_lid.xlsx"))
  expect_equal(input_size_before, input_size_after)
})

test_that("cytation_parse()'s deprecated data_csv argument still works (backlog 048)", {
  dir <- local_fixture_copy("biotek_cytation")

  expect_warning(
    out <- cytation_parse(
      data_csv = file.path(dir, "231031_pr_calibration_lid.xlsx"),
      layout_csv = file.path(dir, "plate_layout.csv"),
      timeseries = FALSE
    ),
    "deprecated"
  )

  expect_true(file.exists(file.path(dir, "231031_pr_calibration_lid_parsed.csv")))
  expect_true(all(c("well", "row", "column") %in% names(out)))
})
