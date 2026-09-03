# Global settings for the seqwrap package
#
# onLoad
.onLoad <- function(...) {
  S7::methods_register()

  # The S7 print method for seqwrapResults leaves a `print` binding in this
  # namespace. R takes that as a locally defined generic and registers any
  # S3method(print, ...) directive from NAMESPACE in the namespace's own
  # methods table, which base::print never consults. Register against the
  # base generic explicitly instead.
  registerS3method(
    "print",
    "seqwrap_priors",
    sw_print_priors,
    envir = asNamespace("base")
  )
}
