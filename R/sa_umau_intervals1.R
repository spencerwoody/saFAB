
##' Selection-adjusted UMAU intervals (one-sided selection region)
##'
##' Selection-adjusted UMAU intervals
##' @title Selection-adjusted UMAU intervals
##' @param y Value for which to make a confidence interval
##' @param thetaRange Range of theta values for which to construct
##'   acceptance region
##' @param sd standard deviation of normal
##' @param t truncation, i.e., S = {i: |y_i| > t}
##' @param alpha Confidence level
##' @return
##' @author Spencer Woody
##'
##' @export
sa_umau_intervals1 <- function(y, thetaRange = c(-15, 15), sd, t,
                              alpha = 0.10) {

  rootFun <- function(theta, y, p, sd, t) {
    y - qtnorm1(p, theta, sd, t)
  }

  lower <- uniroot(f = rootFun, interval = thetaRange,
                   y = y, p = 1 - alpha / 2, sd = sd, t = t)$root  

  upper <- uniroot(f = rootFun, interval = thetaRange,
                   y = y, p = alpha / 2, sd = sd, t = t)$root

  c(lower, upper)

}
