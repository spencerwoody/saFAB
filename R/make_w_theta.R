#' Create optimal spending function for the specified marginal distribution
#'
#' @param marginal_fun function for marginal density of y
#' @param sigma standard deviation for sampling model
#' @param t truncation point
#' @param ... other arguments to marginal function
#' @param alpha confidence level
#' @param theta_min
#' @param theta_max
#' @param num_theta length of theta grid
#' @param epsilon starting point for root finding
#'
make_w_theta <- function(marginal_fun, sigma, t, ...,
                         alpha = 0.05,
                         theta_min = -7, theta_max = 7,
                         num_theta = 5000, epsilon = 1e-10, verbose = F) {

  require(dplyr)

  theta_vec <- seq(theta_min, theta_max, length.out = num_theta)
  w_vec <- rep(NA, num_theta)

  if (verbose) {
    prog <- progress_estimated(num_theta)
  }

  for (i in 1:num_theta) {
    w_vec[i] <- uniroot(Hprime_w_safab,
                        lower = epsilon, upper = 1 - epsilon,
                        theta = theta_vec[i],
                        sigma = sigma, t = t, alpha = alpha,
                        marginal_fun = marginal_fun, ...)$root

    if (verbose) {
      prog$tick()$print()
    }
  }

  return(
    data.frame(
      theta = theta_vec,
      w = w_vec
    )
  )


}
