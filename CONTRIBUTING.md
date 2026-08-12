# Contributing to flopR

## Documentation

This package's `man/*.Rd` files are generated from roxygen2 comments in `R/`
— don't edit files under `man/` directly. After changing any roxygen comment
(`#'` lines) or a function's arguments, regenerate the docs and `NAMESPACE`
before committing:

```r
devtools::document()
```

Commit the resulting `man/` and `NAMESPACE` changes alongside your `R/`
change. `man/` drifting out of sync with the actual function signatures is a
recurring source of confusing bugs for users reading `?function_name`.

## Tests

Run the test suite with:

```r
devtools::test()
```
