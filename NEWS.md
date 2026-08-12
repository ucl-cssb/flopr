# flopr 0.5.0

This release follows a full code review of the package. It bundles several
bug fixes that change computed output, two deliberate breaking behaviour
changes, and one non-breaking API rename with a deprecation shim - kept as
a single release so the breaking changes are clearly flagged in one place.

## Breaking changes

* `calibrate_flu()` now actually uses an exact conversion-factor match when
  one exists in `conversion_factors_csv` at the precisely measured gain,
  instead of always falling through to an interpolated model fit. A bug in
  a `tryCatch(finally = ...)` call meant the exact match was silently
  discarded on every call. **Calibrated fluorescence values may change**
  for anyone whose `conversion_factors_csv` has an exact gain match - this
  is expected and more accurate than the previous behaviour.
* `process_fcs()` now stops with a clear error when no bacteria are found
  in a sample, instead of silently continuing with an empty/undefined
  `flowFrame`. `process_fcs_dir()` now warns and skips the affected file
  (recording it in `data_summary.csv` with `channel =
  "SKIPPED_no_bacteria_found"`) rather than halting the whole batch or
  silently continuing with the negative control's normalisation disabled.
  The previous guard was dead code that never actually fired.
* `infinite_parse()` and `spark_parse()` now default to `timeseries = TRUE`,
  matching `cytation_parse()`/`biotek_parse()`. Previously they defaulted
  to `FALSE`. Any existing call that omits `timeseries` will now parse in
  timeseries mode - pass `timeseries = FALSE` explicitly to keep the old
  behaviour.

## Deprecated

* `cytation_parse()`, `infinite_parse()`, and `spark_parse()`'s first
  argument is renamed from `data_csv` to `data_file` (it always accepted
  `.xls`/`.xlsx` too, so the old name undersold what it takes). The old
  `data_csv` name still works via a deprecation warning and will be
  removed in a future major version. Positional calls are unaffected.
* `sparkParse()` (the camelCase alias for `spark_parse()`) is soft
  deprecated - it still works, but now emits a deprecation warning
  pointing at `spark_parse()`.

## Bug fixes

* `cytation_parse()` and `spark_parse()` no longer overwrite their raw
  `.csv`/`.xls`/`.xlsx` input file with their own parsed output.
* `infinite_parse()` no longer silently no-ops on its default
  (`timeseries = FALSE`) argument.
* `process_plate()` no longer crashes when calibrating more than one
  fluorophore with `to_MEFL = TRUE`.
* `get_calibration()` no longer crashes on a successful multi-channel
  `mef_peaks` fit.
* `generate_cfs()` no longer resolves `filter()` to `stats::filter()`
  instead of `dplyr::filter()` in certain environments.
* `spark_parse()` no longer risks an infinite loop on malformed input, and
  no longer corrupts data on a single-measurement-block file.

## Other changes

* CRAN-readiness: real `Title`/`Description` fields, all declared
  dependencies actually used and all used dependencies declared, package
  data/examples excluded from the built package, runnable `@examples` on
  every exported function, and a clean `R CMD check --as-cran`.
* Added a `testthat`-based regression test suite (previously none existed)
  and a GitHub Actions CI workflow running `R CMD check` on every push/PR.
* Various documentation accuracy fixes and internal maintainability
  clean-up (deduplicated parsing helpers, split oversized files, cleared
  deprecated `dplyr`/`ggplot2` call patterns).
