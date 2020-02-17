#' Use predictive recursion to obtain a
#'
#' @param y sample
#' @param mu0 foobar
#' @param sig0 foobar
#' @param nulltype foobar
#' @return a list, where my_fit is the fitted marginal or y
#'
#' @export
pr_fitter <- function(y, mu0=NULL, sig0=NULL, nulltype = 'theoretical') {
  require(FDRreg)
  require(splines)

  # sampling model of y
  if(nulltype == 'theoretical') {
    if(missing(mu0)) mu0 = 0
    if(missing(sig0)) sig0 = 1
  } else if(nulltype == 'empirical') {
    null_fit <- efron(y, nmids = 200, df = 15,
                     nulltype = 'empirical')
    mu0 <- null_fit$mu0
    sig0 <- null_fit$sig0
  } else {
    error('nulltype must be either "theoretical" or "empirical"')
  }

  # fit predictive recursion
  pr_fit <- prfdr(y, mu0, sig0)

  # interpolate on the log scale to get the prior
  logprior_fit <- splinefun(pr_fit$x_grid,
                            log(pr_fit$pitheta_grid),
                            method='natural')
  prior_fit <- function(theta) exp(logprior_fit(theta))

  # interpolate on the log scale to get the marginal density
  logmy_fit <- splinefun(pr_fit$x_grid, log(pr_fit$fmix_grid), method='natural')
  my_fit <- function(y) exp(logmy_fit(y))

  return(list(
    prior_fit = prior_fit,
    my_fit = my_fit,
    pr_fit = pr_fit,
    sig0 = sig0
  ))

}
