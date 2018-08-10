
#' Density for the doulbe-truncated normal distribution, x ~ N(mean,
#' sd^2) * 1(abs(x) > t)
#'
#' @param x vector of values at which to evaluate density
#' @param mean vector of means
#' @param sd vector of standard deviations
#' @param t truncation point
#' @param log logical; if TRUE, probability given as log(p)
#' @return The density of the double-truncated normal distribution at x
#' @examples
#' dtnorm(3, 0, 1, 1)
#' dtnorm(0.5, 0, 1, 1)
#'
#' @export

dtnorm <- function(x, mean, sd, t, log.p = FALSE) {
  require(dplyr)

  # Tail areas of corresponding (nontruncated) normal distribution
  C <- pnorm((-t - mean) / sd) + 1 - pnorm((t - mean) / sd)

  density <- case_when(abs(x) > t ~ dnorm(x, mean, sd, log) / C,
                       TRUE ~ 0)

  if (log.p) {
    return(log(density))
  } else {
    return(density)
  }

}
