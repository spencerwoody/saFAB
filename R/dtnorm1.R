
#'  Density for the doulbe-truncated normal distribution, x ~ N(mean,
#' sd^2) * 1(abs(x) > t)
#'
#'
#' 
#' @title dtnorm
#' @param x vector of values at which to evaluate density 
#' @param mean vector of means
#' @param sd vector of standard deviations
#' @param t truncation value
#' @param log.p logical; if TRUE, probability given as log(p)
#' @return The density of the double-truncated normal distribution at x 
#' @author Spencer Woody
#'
#' @export
dtnorm1 <- function(x, mean, sd, t, log.p = FALSE) {
  require(dplyr)

  # Tail areas of corresponding (nontruncated) normal distribution
  C <- 1 - pnorm((t - mean) / sd)

  density <- case_when(x > t ~ dnorm(x, mean, sd, log.p),
                       TRUE ~ 0)

  if (log.p) {
    return(density - log(C))
  } else {
    return(density / C)
  }

}
