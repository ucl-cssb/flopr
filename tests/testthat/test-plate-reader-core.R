test_that("generate_cfs() produces conversion factors from Tecan Spark calibration data", {
  dir <- local_fixture_copy("tecan_spark")
  spark_parse(
    data_file = file.path(dir, "191219_calibration_membrane.csv"),
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
    data_file = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )
  generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))
  spark_parse(
    data_file = file.path(dir, "200228_example_data.csv"),
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
  dir <- local_fixture_copy("tecan_spark")
  spark_parse(
    data_file = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )
  cfs <- generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))

  # "GFP 40" (gain 40) is measured directly in the calibration data, so
  # calibrate_flu(flu_gain = 40) should use that exact cf, not an
  # interpolated model fit - this is exactly the path needs-discussion/001
  # found permanently unreachable. Filter by fluorophore first, matching
  # calibrate_flu()'s own logic - "GFP 40" is also a valid measure column
  # for the microspheres calibrant (fluorophore == ""), so a bare
  # measure-only filter picks up both.
  exact_row <- cfs[cfs$measure == "GFP 40" & cfs$fluorophore == "GFP", ]
  expect_equal(nrow(exact_row), 1)

  pr_data <- data.frame(well = "A1", normalised_GFP = 1000)
  out <- calibrate_flu(pr_data, "GFP", flu_gain = 40, "OD700",
                        file.path(dir, "191219_calibration_membrane_parsed_cfs.csv"))
  expect_equal(out$calibrated_GFP, 1000 / exact_row$cf)
})

test_that("calibrate_od() stops with a clear error when od_calib_name has no matching conversion factor", {
  dir <- local_fixture_copy("tecan_spark")
  spark_parse(
    data_file = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )
  generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))

  # od_calib_name has to match the calibration file's "measure" column
  # exactly - a calibration_csv generated with a different raw-data column
  # naming scheme (e.g. a different plate reader protocol) has no matching
  # row, which used to crash deep inside a data.frame assignment
  # ("replacement has 0 rows, data has N") instead of failing clearly.
  pr_data <- data.frame(well = "A1", "normalised_absorbance:600" = 0.5,
                         check.names = FALSE)
  expect_error(
    calibrate_od(pr_data, "absorbance:600",
                 file.path(dir, "191219_calibration_membrane_parsed_cfs.csv")),
    "No conversion factor found"
  )
})

test_that("calibrate_od()/calibrate_flu() succeed with a decoupled calib name (backlog: multicolour experiment naming mismatch)", {
  dir <- local_fixture_copy("tecan_spark")
  spark_parse(
    data_file = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )
  generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))

  # od_name/flu_name are this experiment's own raw column names, which don't
  # have to match the calibration file's measure/fluorophore values as long
  # as od_calib_name/flu_calib_name are given explicitly
  pr_data <- data.frame(well = "A1", "normalised_absorbance:600" = 0.5,
                         "normalised_GFP:488,530" = 1000, check.names = FALSE)
  out <- calibrate_od(pr_data, "absorbance:600",
                       file.path(dir, "191219_calibration_membrane_parsed_cfs.csv"),
                       od_calib_name = "OD700")
  expect_true("calibrated_absorbance:600" %in% names(out))

  out <- calibrate_flu(pr_data, "GFP:488,530", flu_gain = 40, od_name = "absorbance:600",
                        conversion_factors_csv = file.path(dir, "191219_calibration_membrane_parsed_cfs.csv"),
                        flu_calib_name = "GFP")
  expect_true("calibrated_GFP:488,530" %in% names(out))
})

test_that("calibrate_flu() keeps top-read and bottom-read calibration data separate", {
  # calibrate_flu() only reads conversion_factors_csv, never writes to it, so
  # it's safe to reference examples/ directly here rather than copying a
  # fixture - this specific real calibration file is what has both top- and
  # bottom-read rows for the same fluorophore (no equivalent exists yet
  # among the smaller committed fixtures)
  cfs_path <- testthat::test_path(
    "..", "..", "examples", "plate_reader", "biotek_cytation",
    "231031_pr_calibration_membrane_parsed_cfs.csv"
  )
  skip_if_not(file.exists(cfs_path), "examples/ calibration fixture not found")
  cfs <- utils::read.csv(cfs_path)

  # "GFP 40" and "GFP 40 TOP" are both real, distinct calibration rows for
  # the same fluorophore, differing only in whether fluorescence was read
  # from the top or bottom of the well - these must not get blended together
  expect_true(any(cfs$measure == "GFP 40"))
  expect_true(any(cfs$measure == "GFP 40 TOP"))

  pr_data <- data.frame(well = "A1", normalised_GFP = 1000)
  out_bottom <- calibrate_flu(pr_data, "GFP", flu_gain = 40, "OD",
                               cfs_path, flu_calib_name = "GFP")
  out_top <- calibrate_flu(pr_data, "GFP", flu_gain = 40, "OD",
                            cfs_path, flu_calib_name = "GFP TOP")

  expect_equal(out_bottom$calibrated_GFP, 1000 / cfs$cf[cfs$measure == "GFP 40"])
  expect_equal(out_top$calibrated_GFP, 1000 / cfs$cf[cfs$measure == "GFP 40 TOP"])
  expect_false(isTRUE(all.equal(out_bottom$calibrated_GFP, out_top$calibrated_GFP)))
})

test_that("process_plate() normalises/calibrates multiple OD readings without affecting fluorescence", {
  dir <- local_fixture_copy("tecan_spark")
  spark_parse(
    data_file = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )
  generate_cfs(file.path(dir, "191219_calibration_membrane_parsed.csv"))
  spark_parse(
    data_file = file.path(dir, "200228_example_data.csv"),
    layout_csv = file.path(dir, "200228_example_layout.csv"),
    timeseries = TRUE
  )

  args <- list(
    blank_well = c("C12", "D12"), neg_well = c("C6", "D6", "E6"),
    flu_names = c("GFP", "mCherry"), af_model = "spline", to_MEFL = TRUE,
    flu_gains = 135,
    conversion_factors_csv = file.path(dir, "191219_calibration_membrane_parsed_cfs.csv")
  )

  out_single <- do.call(process_plate, c(
    list(data_csv = file.path(dir, "200228_example_data_parsed.csv"), od_name = "OD700"),
    args
  ))

  dir2 <- local_fixture_copy("tecan_spark")
  file.copy(file.path(dir, "200228_example_data_parsed.csv"),
            file.path(dir2, "200228_example_data_parsed.csv"), overwrite = TRUE)
  out_multi <- do.call(process_plate, c(
    list(data_csv = file.path(dir2, "200228_example_data_parsed.csv"),
         od_name = c("OD700", "OD600")),
    args
  ))

  expect_true(all(c("normalised_OD700", "calibrated_OD700",
                     "normalised_OD600", "calibrated_OD600") %in% names(out_multi)))
  # the first od_name is primary - fluorescence normalisation/calibration
  # must be identical whether or not a second OD is also requested
  expect_equal(out_single$normalised_GFP, out_multi$normalised_GFP)
  expect_equal(out_single$calibrated_GFP, out_multi$calibrated_GFP)
})
