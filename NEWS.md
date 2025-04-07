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
