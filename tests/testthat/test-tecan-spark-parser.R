test_that("spark_parse() parses timeseries=FALSE calibration data", {
  dir <- local_fixture_copy("tecan_spark")

  out <- spark_parse(
    data_csv = file.path(dir, "191219_calibration_membrane.csv"),
    layout_csv = file.path(dir, "calibration_plate_layout.csv"),
    timeseries = FALSE
  )

  expect_true(file.exists(file.path(dir, "191219_calibration_membrane_parsed.csv")))
  expect_true(all(c("well", "row", "column", "calibrant", "fluorophore") %in% names(out)))
  expect_snapshot_value(fingerprint(out), style = "json2")
})

test_that("spark_parse() parses timeseries=TRUE sample data", {
  dir <- local_fixture_copy("tecan_spark")

  out <- spark_parse(
    data_csv = file.path(dir, "200228_example_data.csv"),
    layout_csv = file.path(dir, "200228_example_layout.csv"),
    timeseries = TRUE
  )

  expect_true(file.exists(file.path(dir, "200228_example_data_parsed.csv")))
  expect_true(all(c("well", "row", "column", "time") %in% names(out)))
  expect_snapshot_value(fingerprint(out), style = "json2")
})

test_that("spark_parse() input file is never overwritten (backlog 002)", {
  dir <- local_fixture_copy("tecan_spark")
  input_before <- readLines(file.path(dir, "200228_example_data.csv"))

  spark_parse(
    data_csv = file.path(dir, "200228_example_data.csv"),
    layout_csv = file.path(dir, "200228_example_layout.csv"),
    timeseries = TRUE
  )

  input_after <- readLines(file.path(dir, "200228_example_data.csv"))
  expect_identical(input_before, input_after)
})

test_that("spark_parse() errors clearly on a file with only one measurement block (backlog 007)", {
  dir <- local_fixture_copy("tecan_spark")

  # a minimal single-block, non-timeseries file: one "Name" row (the plate-
  # type entry) is not itself a real measurement block, so a file with only
  # that one "Name" row has zero real blocks
  single_block <- data.frame(V1 = "Name", V2 = "Plate info", stringsAsFactors = FALSE)
  single_block_path <- file.path(dir, "single_block.csv")
  utils::write.table(single_block, single_block_path, sep = ",",
                      row.names = FALSE, col.names = FALSE)

  expect_error(
    spark_parse(
      data_csv = single_block_path,
      layout_csv = file.path(dir, "200228_example_layout.csv"),
      timeseries = FALSE
    ),
    "No measurement blocks found"
  )
})
