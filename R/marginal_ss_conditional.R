#' Density of marginal distribution of y ~ N(theta, sigma^2),
#' where theta has spike-slab prior, under joint selection
#'
#' @param y value at which to evaluate
#' @param sigma standard deviation of sampling distribution
#' @param t truncation point
#' @param p proportion for spike part
#' @param tau standard deviation of normal
#'
marginal_ss_conditional <- function(y, sigma, t, p, tau) {

  int_part <- integrate(int_fun_gauss,
                        lower = -Inf,
                        upper = +Inf,
                        y = y, sigma = sigma, tau = tau, t = t)$value

  density <- (p * int_part +
                (1 - p) * dnorm(y, 0, sigma) / Pr_S_conditional(0, sigma, t)) *
    ifelse(abs(y) > t, 1, 0)

  return(density)

}
