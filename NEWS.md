# seqwrap (development version)

# seqwrap 0.7.0

- Updated documentation for CRAN submission. 

# seqwrap 0.6.3

- CRAN-readiness pass:
  * Title case in DESCRIPTION ("Item-by-Item").
  * Bumped `Depends: R (>= 4.1)` to reflect use of the native pipe in vignettes.
  * Replaced a non-canonical CRAN URL in the main vignette.
  * Added `@return` documentation to `seqwrap_summarise()`, `simcounts()`,
    `seqwrapResults`, and `swcontainer`.
  * Marked internal helpers (`seqwrap_check`, `fit_fun`, `fit_fun_lme`,
    `data_helper`, `seqwrap_mtf`) with `@noRd` so they no longer generate
    user-facing Rd pages.
  * Replaced `seqwrap:::simcounts` with `seqwrap::simcounts` in tests.
  * Fixed handling of `cores = "max"` and `cores = NULL` in `seqwrap()`;
    invalid values now error early.
  * Vignettes gracefully skip evaluation when optional Suggests packages
    (including Bioconductor `edgeR`) are unavailable.
  * Minor typo and roxygen-tag fixes.
- Removed the Pillon data set (moved to trainome/seqwrappaper)


# seqwrap 0.6.2

- Bugs
  * Fixed a bug in printing seqwrap_results where the variable k was not declared in the function environment. 
  * Fixed SUGGEST dependencies for Vignette build.

# seqwrap 0.6.1

- Updated the README


# seqwrap 0.6.0

- Bugs
  * Internal fix, passing "library(broom.mixed)" to clusterEvalQ is not recommended. Changed to requireNamespace("seqwrap").
  * Internal fix, removed triple colon on exported functions in vignette. 
  * Removed getwd() in saving model output on disk (seqwrap_mtf.R).
  * Updated @importFrom stats in simcounts2.R


# seqwrap 0.5.0

## News

- Bugfixes

# seqwrap 0.4.0

## News

-   Added support for `nlme::lme` and `nlme::gls`, see the vignette for the use of `additional_vars` when working with `nlme::gls`
-   `targetdata` now supports a list of data frames. Target-specific data frames are made available by their column names.

# seqwrap 0.2.0

This update has focused on improving the workflow of using seqwrap. A new function (`seqwrap_compose`) allows for the user to collect all data needed to iterate over targets to fit models. `seqwrap` can still be specified using the same arguments, but can also use a `swcontainer` created with `seqwrap_compose`. In `seqwrap` we only need to specify e.g. the number of cores and return/save models.

## Breaking changes

-   The `fittin_fun` argument in seqwrap has been replaced by `modelfun`.

## New features

-   `seqwrap_compose` let's you collect all data elements and arguments needed to run iterative modelling with `seqwrap` without initializing it.

-   `seqwrap:::simcounts` was created as an internal function used in testing. It creates a simulated data set of counts based on variation across genes in a set of parameters.

-   A new set of classes and methods has been written using the S7 OOP system. This means that data is validated to prevent errors in setting up the models and data.

-   `targetdata` is now available as an argument in `seqwrap` and `seqwrap_compose` making it possible to supply target-wise values used in the arguments. E.g. setting the dispersion parameter to a fixed value.

-   Using a `swcontainer` object as the first argument in `seqwrap` followed by a named argument will lead to an update of the `swcontainer` object before any modelling.

-   `seqwrap_summarise` efficiently combine data frames from summary and evaluation functions.

-   and more...

## Known limitations

-   
