# flopr 0.6.0

Prompted by a real crash: `process_plate(..., to_MEFL = TRUE)` could fail
with an opaque `Error in \`$<-.data.frame\`(...): replacement has 0 rows`
whenever a calibration run's raw column naming didn't match the current
experiment's - `od_name`/`flu_names` previously had to do double duty as
both the raw column name *and* the calibration lookup key. This release
decouples those two roles and adds support for multiple simultaneous OD
readings (e.g. OD600 and OD700).

## Breaking changes

* All OD output columns now use the same `normalised_<name>`/
  `calibrated_<name>` pattern fluorescence columns already used, instead of
  the previous fixed `normalised_OD`/`calibrated_OD` names. **Any script
  reading those exact column names needs updating** - e.g.
  `od_name = "OD700"` now produces `normalised_OD700`/`calibrated_OD700`,
  not `normalised_OD`/`calibrated_OD`. The `_OD.pdf` diagnostic plot
  filename is similarly now `_<od_name>.pdf` per OD reading.

## New features

* `process_plate()`'s `od_name` can now be a vector, to normalise/calibrate
  more than one OD reading in a single call (e.g.
  `od_name = c("OD600", "OD700")`). The first element is always the
  "primary" OD - the one autofluorescence correction is modelled against -
  regardless of how many are given; the rest are normalised/calibrated
  alongside it without affecting fluorescence normalisation.
* New `od_calib_names`/`flu_calib_names` arguments to `process_plate()`
  (and `od_calib_name`/`flu_calib_name` on the internal `calibrate_od()`/
  `calibrate_flu()`) let the calibration lookup key be set independently of
  the experiment's own raw column name. Both default to the raw name, so
  nothing changes unless you use them - only needed when
  `conversion_factors_csv` was generated from a calibration run that used
  different column naming than this experiment.
* `flu_calib_name`/`flu_calib_names` can carry a trailing `" TOP"` (e.g.
  `"GFP TOP"`) to select top-read rather than bottom-read calibration data,
  for plate readers that support reading fluorescence from either side of
  the well. Previously `calibrate_flu()`'s gain-interpolation fit blended
  top-read and bottom-read calibration rows for the same fluorophore
  together into one curve (or crashed outright, since a `"... TOP"` measure
  string doesn't parse as a plain gain number) - they're now kept separate.

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
