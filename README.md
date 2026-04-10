# seqwrap - Item-By-item Iterative Model Fitting


The `seqwrap` R package allows you to model high-dimensional data using
an item-by-item strategy with a specific model engine. Using `seqwrap`,
you can fit, e.g., gene/transcript data from an RNA sequencing
experiment to a model specified with random effects or custom
distributions. This is possible because `seqwrap` efficiently iterates
over all items (e.g., genes) and fits the data to the same model
formulation using established R packages for regression modeling.

See the package Vignettes for examples.

## Installing

`seqwrap` can be installed from GitHub:

``` r
remotes::install_github("trainome/seqwrap")
```
