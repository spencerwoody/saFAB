
## TODO: edit Hprime function to logit-scale w?

##' Create optimal spending function for the specified marginal distribution
##' 
##' @param marginal_fun function for marginal density of y
##' @param sigma standard deviation for sampling model
##' @param t truncation point
##' @param ... additional arguments to marginal function
##' @param alpha confidence level
##' @param theta_min minimal theta to create spending function
##' @param theta_max maximum theta to create spending function
##' @param num_theta length of theta grid
##' @param epsilon starting point for root finding
##' @export
make_w_theta1 <- function(marginal_fun, sigma, t, ...,
                         alpha = 0.05,
                         theta_min = -7, theta_max = 7,
                         num_theta = 5000, epsilon = 1e-10,
                         verbose = FALSE) {

  require(dplyr)

  ## Define grid of theta
  theta_vec <- seq(theta_min, theta_max, length.out = num_theta)
  w_vec <- rep(NA, num_theta)

  if (verbose) prog <- progress_estimated(num_theta)

  for (i in 1:num_theta) {

    ## The w for this theta is find by solving for root of derivative
    ## for loss function
    w_vecI <- try(
      uniroot(Hprime_w_safab1,
              lower = epsilon, upper = 1 - epsilon,
              theta = theta_vec[i],
              sigma = sigma, t = t, alpha = alpha,
              marginal_fun = marginal_fun, ...)$root,
      silent = TRUE
    )

    ## If there's an error, it means that the root finder failed;
    ## probably because the real root is close to 
    if (inherits(w_vecI, "try-error")) {
      if (theta_vec[i] < t) {
        w_vecI <- epsilon
      } else {
        w_vecI <- 1 - epsilon
      }
    }
    
    w_vec[i] <- w_vecI

    if (verbose) prog$tick()$print()
    
  }

  if (verbose) cat("\n")

  return(
    data.frame(
      theta = theta_vec,
      w = w_vec
    )
  )


}
