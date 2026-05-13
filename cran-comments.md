## Submission

This is a new submission of the `seqwrap` package.

`seqwrap` provides an item-by-item wrapper for fitting user-specified
regression models to high-dimensional data (e.g. RNA-seq read counts),
delegating the fit to any model-fitting function the user supplies
(`stats::lm`, `glmmTMB::glmmTMB`, `nlme::lme`, etc.).

## Test environments

* Local: Windows 11, R 4.5.2 — 0 errors, 0 warnings, 0 notes
  (other than the expected "new submission" note).
* win-builder: R-devel and R-release (planned before release).
* macOS-builder: R-release (planned before release).

## R CMD check results

0 errors | 0 warnings | 1 note

* The only note is the standard "New submission" note.

## Downstream dependencies

There are no reverse dependencies.

## Notes on package contents

* The package ships two real RNA-seq example data sets (`dungan_counts`,
  `pillon_counts`) sourced from GEO (GSE195707, GSE202295) for use in the
  vignettes and documentation; both are stored as xz-compressed `.rda`
  files and used to demonstrate the API.
* Vignettes use a setup chunk that checks for optional dependencies
  (`edgeR`, `gt`, `cowplot`, `ggtext`, etc., several of which are
  Bioconductor packages) and gracefully skips evaluation when they are
  not installed on the build machine.
* Parallel computation in vignettes is capped at 2 cores and respects
  `_R_CHECK_LIMIT_CORES_`.
