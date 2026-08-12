test_that("biotek_parse() parses example data (timeseries=TRUE default)", {
  dir <- local_fixture_copy("biotek_neo2")

  out <- biotek_parse(
    data_xl = file.path(dir, "biotek-neo2_example_data.xlsx"),
    layout_csv = file.path(dir, "layout.csv")
  )

  expect_true(file.exists(file.path(dir, "biotek-neo2_example_data_parsed.csv")))
  expect_true(all(c("well", "row", "column") %in% names(out)))
  expect_snapshot_value(fingerprint(out), style = "json2")
})

test_that("biotek_parse(timeseries=FALSE) raises a clear error rather than a silent no-op (backlog 008 precedent)", {
  dir <- local_fixture_copy("biotek_neo2")

  expect_error(
    biotek_parse(
      data_xl = file.path(dir, "biotek-neo2_example_data.xlsx"),
      layout_csv = file.path(dir, "layout.csv"),
      timeseries = FALSE
    ),
    "only parse timeseries"
  )
})
