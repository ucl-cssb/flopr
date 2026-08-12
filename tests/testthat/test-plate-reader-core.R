test_that("generate_cfs() produces conversion factors from Tecan Spark calibration data", {
  dir <- local_fixture_copy("tecan_spark")
  spark_parse(
    data_csv = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )

  cfs <- generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))

  expect_true(all(c("cf", "beta", "calibrant", "fluorophore", "measure") %in% names(cfs)))
  expect_true(nrow(cfs) > 0)
  expect_snapshot_value(fingerprint(cfs), style = "json2")
})

test_that("process_plate() handles multi-fluorophore + to_MEFL=TRUE (backlog 004, the README's own example)", {
  dir <- local_fixture_copy("tecan_spark")

  spark_parse(
    data_csv = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )
  generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))
  spark_parse(
    data_csv = file.path(dir, "200228_example_data.csv"),
    layout_csv = file.path(dir, "200228_example_layout.csv"),
    timeseries = TRUE
  )

  # this is the exact call from the README - used to crash before backlog 004
  out <- process_plate(
    data_csv = file.path(dir, "200228_example_data_parsed.csv"),
    blank_well = c("C12", "D12"),
    neg_well = c("C6", "D6", "E6"),
    od_name = "OD700",
    flu_names = c("GFP", "mCherry"),
    af_model = "spline",
    to_MEFL = TRUE,
    flu_gains = 135,
    conversion_factors_csv = file.path(dir, "191219_calibration_membrane_parsed_cfs.csv")
  )

  expect_true("calibrated_GFP" %in% names(out))
  expect_true("normalised_mCherry" %in% names(out))
  expect_true(file.exists(file.path(dir, "200228_example_data_parsed_processed.csv")))
  expect_snapshot_value(fingerprint(out), style = "json2")
})

test_that("calibrate_flu() uses an exact gain match when available (needs-discussion/001)", {
  skip("Blocked on backlog 045 - calibrate_flu()'s tryCatch(finally=...) bug always
        discards the exact match today. Un-skip once 045 lands.")

  dir <- local_fixture_copy("tecan_spark")
  spark_parse(
    data_csv = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )
  cfs <- generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))

  # "GFP 40" (gain 40) is measured directly in the calibration data, so
  # calibrate_flu(flu_gain = 40) should use that exact cf, not an
  # interpolated model fit - this is exactly the path needs-discussion/001
  # found permanently unreachable.
  exact_row <- cfs[cfs$measure == "GFP 40", ]
  expect_equal(nrow(exact_row), 1)

  pr_data <- data.frame(well = "A1", normalised_GFP = 1000)
  out <- calibrate_flu(pr_data, "GFP", flu_gain = 40, "OD700",
                        file.path(dir, "191219_calibration_membrane_parsed_cfs.csv"))
  expect_equal(out$calibrated_GFP, 1000 / exact_row$cf)
})
