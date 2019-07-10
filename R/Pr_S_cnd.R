#' Probability of selection event, conditonal on theta
#'
#' @param theta the given value of mean of sampling distribution for y
#' @param sigma standard deviation for y
#' @param t truncation point
#' 
#' 
#'
Pr_S_cnd <- function(theta, sigma, t) {
  pnorm((-t - theta) / sigma) + 1 - pnorm((t - theta) / sigma)
}
