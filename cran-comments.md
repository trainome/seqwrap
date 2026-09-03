## R CMD check results

-- R CMD check results --------------------- seqwrap 0.8.1 ----
Duration: 3m 25.9s

0 errors v | 0 warnings v | 0 notes v

Checked with `R CMD check --as-cran` on Windows 11, R 4.6.0, including the
vignettes, the tests and the PDF manual.

## Notes for the reviewer

* This is an update to a package already on CRAN (0.7.0). It has no reverse
  dependencies. The submission covers two development releases, 0.8.0 and
  0.8.1, which are both described in `NEWS.md`.

* The release changes several defaults and the layout of two slots in the
  returned object, in order to keep memory use bounded when iterating over very
  large target sets. The changes and the corresponding migration steps are
  documented in `NEWS.md` and in a dedicated vignette,
  `migrating-from-seqwrap-0-7`.

* New functions `seqwrap_priors()` and `dispersion_evaluation()` build
  empirical Bayes priors for `glmmTMB` from a completed run. `graphics` and
  `grDevices` have been added to Imports for the `plot()` method on the
  priors object.

* The "checking Rd contents" NOTE reported for 0.7.0 on the Debian clang and
  gcc flavours (Rd files without a usage section: `print.seqwrapResults.Rd`,
  `seqwrapResults.Rd`, `swcontainer.Rd`) is fixed. The two S7 class
  constructors now document their usage, and the `print()` method no longer
  carries an arguments section.

* `DHARMa` has moved from Imports to Suggests. It is now used only by
  `residual_diagnostics()`, which is guarded by `requireNamespace()`.

* There are no published references describing the methods in this
  package, so no `<doi:...>` / `<arXiv:...>` references are included
  in the DESCRIPTION. The package provides a general-purpose iteration
  framework for fitting user-supplied models target-by-target rather
  than implementing a specific published method.

* One of the three vignettes (`fitting-models-with-seqwrap`) is shipped
  pre-built via the `R.rsp::asis` engine because a fresh re-render
  fits many models in parallel and would substantially exceed the
  CRAN check time budget. The corresponding `.qmd` source is excluded
  from the build via `.Rbuildignore`. The other vignettes
  (`fitting-lme4-nlme-models-with-seqwrap` and `migrating-from-seqwrap-0-7`)
  are built during `R CMD check`.

* Examples, tests and vignettes built during check use at most two cores.
