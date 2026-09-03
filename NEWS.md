# seqwrap 0.8.1

## New features
- New `seqwrap_priors()` builds empirical Bayes priors for `glmmTMB` from a
  completed run. Fixed-effect estimates are summarised across targets into
  `normal()` priors (centered on zero by default, or on the mean across
  targets), random-effect standard deviations into a `gamma()` prior, and the
  log dispersion is modelled as a loess trend against log mean count so that
  each target receives its own dispersion prior. The result is a list with one
  prior data frame per target, ready to be passed as `targetdata` to
  `seqwrap_compose()`, and prints a summary of the priors it holds. This
  replaces the manual construction shown in earlier versions of the vignette.
  A `plot()` method, also reachable through `seqwrap_priors(plot = TRUE)`,
  draws each prior over the estimates it was built from: histograms with the
  prior density for fixed effects and random-effect standard deviations, and
  the log dispersion against log mean count with the fitted trend. The
  per-target estimates behind the panels are kept in the `data` attribute of
  the object so that other displays can be built from them.
- New `dispersion_evaluation()` is an evaluation function for `glmmTMB` models
  that records the log dispersion, its standard error and the log mean count of
  each target, which is what `seqwrap_priors()` needs to build the dispersion
  prior.

- A list passed as `targetdata` can now be named by target identifier, in
  which case elements are matched to targets by name regardless of order, and
  the list may hold more targets than `data` (for example priors built from a
  full run applied to a subset). `seqwrap_priors()` returns such a list. An
  unnamed list is still matched by position, as before, but `seqwrap()` now
  warns that it does so, since a misordered list is otherwise silent. A named
  list that lacks a target is rejected with a message naming the missing
  targets.

## Improvements and bug fixes
- `seqwrap()` now reports, in its verbose header, how many times the progress
  bar will advance and in what step size. The bar moves once per batch of
  chunks (one chunk per core), so with the default chunk size it advances in
  about four steps; the message points to `chunk_size` for finer updates.
- Target-wise data (`targetdata`) was paired with the wrong targets when the
  target identifiers in `data` were not in sorted order. Targets are grouped by
  sorted identifier internally, while target-wise data was taken by row, so row
  `j` of `targetdata` ended up with the `j`-th target in sorted order rather
  than with row `j` of `data`. Target-wise data now follows the same
  reordering as the data, for data frames and lists, with and without chunking
  or `subset`.

# seqwrap 0.8.0

## Improvements and bug fixes
- `seqwrap()` now dispatches data to workers in chunks of targets rather than
  one target at a time. The new `chunk_size` argument controls
  the granularity and defaults to a value chosen from the number of targets and
  workers. Per-target data frames are built on the worker, so the full
  per-target expansion of the data no longer exists in the parent process.
  Results are unchanged by chunking.
- New `cache` argument to `seqwrap()` controlling how summaries and evaluations
  are accumulated:
  * `"none"` is the previous behavior, one data frame per target.
  * `"memory"` (default) binds them into a single data frame per slot while 
  still on the worker. This reduces memory usage.
  * `"disk"` additionally writes each chunk to `cache_path` (the session
    temporary directory by default), keeping memory flat regardless of the
    number of targets.
  `seqwrap_summarise()` reads all three representations transparently.
- `seqwrap()` and `seqwrap_compose()` does no longer attempt to support 
  a `DGEList`. An informative message is given if attempting to use a `DGEList`
  as input in `seqwrap()` or `seqwrap_compose()`. `edgeR` remains in Suggests, 
  since the vignettes use it for library size normalisation.
- `seqwrap()` and `seqwrap_compose()` with `rownames = TRUE` produced wrong 
  results. The setting was stored on the container but never passed to the 
  internal data preparation step, so the first sample column was used as the 
  target identifier.
- `exported` objects from `seqwrap()` and `seqwrap_compose()` were not available 
  to `summary_fun` and `eval_fun`. The list was sent to the workers under the 
  name `exported` rather than having its contents placed there. Elements of 
  `exported` must now be named, which is checked.
- Using `subset = integer(0)` in `seqwrap` failed when constructing the 
  results object. Missing or blank target identifiers are now rejected with 
  a clear message instead of failing inside a worker with "argument is of 
  length zero".
- the help pages for `seqwrapResults`, `swcontainer` and the `print()` method 
  no longer document arguments without a `\usage` section, which raised a NOTE 
  in the CRAN checks for 0.7.0.
- The check for whether the fitting function is `nlme::lme` is now made once per
  run rather than once per target which reduces the computational cost.
- `save_models = TRUE` wrote every model to a single file named `.rds`, 
  overwriting it once per target, because the target identifier was not
  available inside the per-target worker. Models are now written to one file per
  target, named after the target, with characters that are not valid in file
  names replaced.
- `seqwrap()` now warns when `return_models = TRUE` would retain more than ten
  models: "You are saving N model objects which could be very memory intensive.
  If this was intended, ignore this warning." Prototyping on a small subset,
  which is what `return_models` is meant for, stays quiet.
- New `targets` property on `seqwrapResults` lists every target identifier in
  fitting order. This is the index of the run, and remains available
  when models are not returned and summaries are cached to disk. Targets that
  completed without issues are found by 
  `setdiff(x@targets, seqwrap_errors(x)$target)`.
- `seqwrap_summarise()` now has `drop_warnings`, defaulting to FALSE so that
  targets which raised a warning are still reported. Set it to TRUE to remove
  targets that warned at any stage, or pass a character vector of stages
  (`"fit"`, `"summary"`, `"evaluation"`). The identifiers of removed targets are
  returned in a `dropped` element and the count is reported, so the choice is
  visible and reversible. Note that targets which warn are do so because of 
  different mechanisms. The `@errors` slot contain error and warning codes from 
  the specific fitting algorithm.
- `seqwrap_summarise(errors = TRUE)` now returns the errors it announced. The
  argument previously printed "Attempting to summarise errors and warnings" and
  returned nothing. The new `errors` element holds one row per condition raised,
  with columns `target`, `stage`, `type` and `message`, covering only the
  targets that raised something.
- When a warning was raised by the fitting algorithm `seqwrap_mtf()` (internal)
  discarded the fitted model. Conditions are now recorded with 
  `withCallingHandlers()`, which resumes execution, so only a hard errors leads
  to lost results. The same fix applies to warnings raised inside `summary_fun`
  and `eval_fun`.
- Summary and evaluation functions are no longer called when model fitting
  fails. Passing NULL to a summary function tended to produce a degenerate
  zero-column result rather than an honest NULL: `broom.mixed::tidy(NULL)`
  returns a 0 x 0 tibble, and such targets then disappeared silently from
  combined summaries. A failed fit now yields NULL in both slots, with the
  cause recorded in `errors`.
- Combining results no longer fails when per-target data frames have different
  columns. `rbind()` errors outright in that case, so a handful of targets
  returning a differently shaped summary could discard an entire run's results,
  or abort the run itself when caching was enabled. Missing columns are filled
  with NA instead, and zero-column results are dropped.
- `generic_evaluation()` no longer simulates residuals with DHARMa. It now
  reads the convergence diagnostics that the fitting algorithm itself reports,
  which is roughly 400 times cheaper per target: for 800,000 targets on a
  single core this is about 10 minutes rather than about 69 hours.
- `generic_evaluation()` is now an S3 generic, so the diagnostics are read in
  whatever form the engine supplies them. Methods are provided for `glmmTMB`,
  `lme4` (`merMod`, covering `lmer` and `glmer`), `nlme` (`lme` and `gls`),
  `mgcv` (`gam` and `bam`), `MASS::glm.nb` (`negbin`), `glm` and `lm`, with a
  default method for anything else. All methods return the same columns 
  (`engine`, `converged`, `code`, `message`, `singular`, `iterations`) 
  so results remain compatible between runs.
  * `converged` and `singular` are reported separately. A fit can satisfy the
    optimiser while still being degenerate: `glmmTMB` regularly returns
    convergence code 0 alongside a non-positive-definite Hessian, and `lme4`
    returns code 0 for singular random effect structures. Screening large sets
    of fits usually means requiring `converged & !singular`.

- DHARMa moved from Imports to Suggests; `mgcv` added to Suggests.
  
## New features
- New `cache` property on `seqwrapResults` describing the on-disk cache.
- `seqwrap_errors()` was added to make it easier to handle errors and 
  warnings. The function validates `stage` and `type` even when the run raised
  no conditions.
- New `residual_diagnostics()` provides the previous DHARMa based behaviour (see
  changes to `generic_evaluation()` above) as an opt-in evaluation function, 
  with a configurable number of simulations.
- New `seqwrap_cache_clear()` removes the chunk files written by
  `cache = "disk"`.


## Breaking changes

- `seqwrap` argument `return_models` now defaults to FALSE. Retaining every
  fitted model is the largest memory cost a run can incur, and it was easy to
  leave on unintentionally after prototyping. Pass `return_models = TRUE`
  explicitly when the fitted objects are needed, or use `save_models` with
  `model_path` to keep them on disk instead.
- `seqwrap` used with `cache` set to the default `"memory"` changes `@summaries`
  and `@evaluations` to hold a single data frame with a `target` column rather
  than one data frame per target. 
  * `seqwrap_summarise()` returns exactly the same combined data frames as
    before, so the documented workflow is unaffected. Only code reading
    `@summaries` or `@evaluations` directly needs changing:
    `x@summaries[["gene1"]]` becomes
    `x@summaries[x@summaries$target == "gene1", ]`.
  * When `cache = "none"` restores the previous per-target lists. This is also
    required when a summary or evaluation function returns something that
    cannot be row-bound, such as a list or a model object; combining now fails
    with an explicit message naming the offending target.

- The `errors` slot of `seqwrapResults` is now in long form, simplified data 
  frame holding one row per condition raised with the columns `target`, 
  `stage`, `type`, `class` and `message`. It previously held six list columns of
  condition objects, with an entry for every target whether or not anything
  went wrong. The new shape is both easier to filter and much smaller.
  Code of the form `x@errors$warnings_fit[[i]]` should become
  `seqwrap_errors(x, stage = "fit", type = "warning")`.
  * Condition objects are no longer retained. They cost roughly 1.6 kB each
    because they carry their originating call. The message and the condition
    class are kept instead.
  * The new `class` column makes it possible to separate specific warning types
    from generic ones without matching on message text.

- The columns in the `evaluations` slot differ when `eval_fun`
  is left as NULL. Code that read `uniformity`, `dispersion` or `outliers` from
  the default evaluation should pass `eval_fun = residual_diagnostics`.





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
  * Fixed a bug in printing seqwrap_results where the variable k was not 
  declared in the function environment. 
  * Fixed SUGGEST dependencies for Vignette build.

# seqwrap 0.6.1

- Updated the README


# seqwrap 0.6.0

- Bugs
  * Internal fix, passing "library(broom.mixed)" to clusterEvalQ is not
  recommended. Changed to requireNamespace("seqwrap").
  * Internal fix, removed triple colon on exported functions in vignette. 
  * Removed getwd() in saving model output on disk (seqwrap_mtf.R).
  * Updated @importFrom stats in simcounts2.R


# seqwrap 0.5.0



# seqwrap 0.4.0

## News

-   Added support for `nlme::lme` and `nlme::gls`, see the vignette for the 
use of `additional_vars` when working with `nlme::gls`
-   `targetdata` now supports a list of data frames. Target-specific data 
frames are made available by their column names.

# seqwrap 0.2.0

This update has focused on improving the workflow of using seqwrap. A new 
function (`seqwrap_compose`) allows for the user to collect all data needed 
to iterate over targets to fit models. `seqwrap` can still be specified using
the same arguments, but can also use a `swcontainer` created with 
`seqwrap_compose`. In `seqwrap` we only need to specify e.g. the number of 
cores and return/save models.

## Breaking changes

- The `fittin_fun` argument in seqwrap has been replaced by `modelfun`.

## New features

- `seqwrap_compose` let's you collect all data elements and arguments needed
to run iterative modelling with `seqwrap` without initializing it.

- `seqwrap:::simcounts` was created as an internal function used in testing. 
  It creates a simulated data set of counts based on variation across genes in
  a set of parameters.

- A new set of classes and methods has been written using the S7 OOP system. 
This means that data is validated to prevent errors in setting up the models 
and data.

- `targetdata` is now available as an argument in `seqwrap` and 
`seqwrap_compose` making it possible to supply target-wise values used in the 
arguments. E.g. setting the dispersion parameter to a fixed value.

- Using a `swcontainer` object as the first argument in `seqwrap` followed 
by a named argument will lead to an update of the `swcontainer` object before
any modelling.

- `seqwrap_summarise` efficiently combine data frames from summary and 
evaluation functions.


