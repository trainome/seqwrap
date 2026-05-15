## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes for the reviewer

* There are no published references describing the methods in this
  package, so no `<doi:...>` / `<arXiv:...>` references are included
  in the DESCRIPTION. The package provides a general-purpose iteration
  framework for fitting user-supplied models target-by-target rather
  than implementing a specific published method.

* One of the two vignettes (`fitting-models-with-seqwrap`) is shipped
  pre-built via the `R.rsp::asis` engine because a fresh re-render
  fits many models in parallel and would substantially exceed the
  CRAN check time budget. The corresponding `.qmd` source is excluded
  from the build via `.Rbuildignore`. The second vignette
  (`fitting-lme4-nlme-models-with-seqwrap`) is built live by quarto
  during `R CMD check`.
