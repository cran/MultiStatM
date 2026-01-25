# Deprecated and removed functions
# These stubs are kept for backward compatibility and provide
# informative error messages. They will be removed in a future release.


#' Deprecated function
#'
#' `SampleGC()` has been replaced by the more efficient
#' \code{GramCharlier()} function.
#'
#' @keywords deprecated internal
#' @export
SampleGC <- function(...) {
  stop(
    "`SampleGC()` has been replaced by the more efficient `GramCharlier()`.\n",
    "Please update your code.",
    call. = FALSE
  )
}


#' Deprecated function
#'
#' `SampleEdg()` has been replaced by the more efficient
#' \code{Edgeworth()} function.
#'
#' @keywords deprecated internal
#' @export
SampleEdg <- function(...) {
  stop(
    "`SampleEdg()` has been replaced by the more efficient `Edgeworth()`.\n",
    "Please update your code.",
    call. = FALSE
  )
}


#' Deprecated function
#'
#' `SampleHermiteN()` has been removed. The \code{HermiteN()} function now
#' supports both univariate and multivariate Hermite polynomials up to order N.
#'
#' @keywords deprecated internal
#' @export
SampleHermiteN <- function(...) {
  stop(
    "`SampleHermiteN()` has been removed.\n",
    "Use `HermiteN()` instead, which now supports univariate and multivariate cases.",
    call. = FALSE
  )
}


#' Deprecated function
#'
#' `UnivMomCum()` has been renamed to \code{MargMomCum()}.
#'
#' @keywords deprecated internal
#' @export
UnivMomCum <- function(...) {
  stop(
    "`UnivMomCum()` has been renamed to `MargMomCum()`.\n",
    "Please update your code.",
    call. = FALSE
  )
}
