#' Derivative of objective function to be minimized, as described by...
#'
#' @param w weight corresponding to theta; i.e., w(theta)
#' @param theta the mean of sampling distribution y
#' @param sigma standard deviation of sampling distribution for y
#' @param t the truncation, i.e. y ~ N(theta, sigma) * 1(abs(y) > t)
#' @param alpha confidence level
#' @param marginal_fun the function for the marginal distribution
#' @param ... optional arguments to the marginal function
#'
#' @export
Hprime_w_safab <- function(w, theta, sigma, t, alpha, marginal_fun, ...) {


  F_inv_1 <- qtnorm(alpha * w + 1 - alpha, theta, sigma, t)
  F_inv_2 <- qtnorm(alpha * w            , theta, sigma, t)

  marginal_fun(F_inv_1, sigma = sigma, t = t, ...)  /
    dtnorm(F_inv_1, theta, sigma, t) -
    marginal_fun(F_inv_2, sigma = sigma, t = t, ...) /
    dtnorm(F_inv_2, theta, sigma, t)

}
