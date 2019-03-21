#' Quantile for the double-truncated normal distribution, x ~ N(mean,
#' sd^2) * 1(abs(x) > t)
#'
#' @param x vector of values at which to evaluate density
#' @param mean vector of means
#' @param sd vector of standard deviations
#' @param t truncation point
#' @param log.p logical; if TRUE, probability given as log(p)
#' @return The density of the double-truncated normal distribution at x
#' @examples
#' qtnorm(3, 0, 1, 1)
#' qtnorm(0.5, 0, 1, 1)
#'
#' @export
qtnorm <- function(p, mean, sd, t, log.p = FALSE) {
  require(dplyr)

  if (log.p) p <- exp(p)

  # Tail areas of corresponding (nontruncated) normal distribution
  C <- pnorm((-t - mean) / sd) + 1 - pnorm((t - mean) / sd)

  case_when(p <= pnorm((-t - mean) / sd) / C ~ mean + sd * qnorm(C * p),
            TRUE ~ mean + sd * qnorm(C * p +
                                       pnorm((  t - mean) / sd) -
                                       pnorm(( -t - mean) / sd)))

}


