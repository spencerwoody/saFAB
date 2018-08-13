
# FAB functions ---------------------------------------------------------

# Derivative of Bayes risk (expected Lebesque measure) for spike-slab prior,
# w.r.t. w, for calculating w(theta), the optimal spending functions
Hprime_fab <- function(w, theta, p, sigma, tau, alpha) {
  
  p * sigma * alpha / sqrt(sigma^2 + tau^2) * 
    exp(- {theta - sigma * qnorm(alpha * (1 - w))}^2 / (2 * sigma^2 + tau^2) + 
    {qnorm(alpha * (1 - w))}^2 / 2) +
    (1 - p) * alpha * 
    exp(- {theta - sigma * qnorm(alpha * (1 - w))}^2 / (2 * sigma^2) + 
    {qnorm(alpha * (1 - w))}^2 / 2) - 
    p * sigma * alpha / sqrt(sigma^2 + tau^2) * 
    exp(- {theta - sigma * qnorm(1 - alpha * w)}^2 / (2 * sigma^2 + tau^2) + 
    {qnorm(1 - alpha * w)}^2 / 2) -
    (1 - p) * alpha * 
    exp(- {theta - sigma * qnorm(1 - alpha * w)}^2 / (2 * sigma^2) + 
    {qnorm(1 - alpha * w)}^2 / 2)
  
}

##' Make spending function for FAB setting, with spike-slab prior
##'
##' 
##' @title make_w_fab
##' @param sigma standard deviation of sampling distribution
##' @param p prior probabilitiy of signal
##' @param tau standard deviation of Gaussian slab
##' @param alpha confidence level
##' @param theta_min minimal value of theta to create spending function
##' @param theta_max maximal value of theta to create spending function
##' @param num_theta number of thetas over which to create grid 
##' @param epsilon tolerance level
##' @param verbose logical; if TRUE, print progress bar
##' @return a vector of w's, corresponding to the spending function
##' @author Spencer Woody
##'
##' @export
make_w_fab <- function(sigma, p, tau, alpha = 0.1, theta_min = -7,
                       theta_max = 7, num_theta = 5000,
                       epsilon = 1e-5, verbose = F) {
  require(dplyr)

  theta_vec <- seq(theta_min, theta_max, length.out = num_theta)
  
  w_theta_fab <- rep(NA, num_theta)
  
  prog <- progress_estimated(num_theta)
  
  for (i in 1:length(theta_vec)) {
    
    w_theta_fab[i] <- uniroot(f = Hprime_fab, 
                              interval = c(epsilon, 1 - epsilon), 
                              theta = theta_vec[i],
                              p = p,
                              sigma = sigma, tau = tau, alpha = alpha)$root
    if (verbose) {
      prog$tick()$print()
    }
  }
  
  return(data.frame(theta = theta_vec,
                    w = w_theta_fab))
  
}

##' Lower bound of FAB interval
##'
##' 
##' @title theta_l_fab
##' @param y single observation
##' @param theta_vec vector of thetas
##' @param w_theta vector of w's for thetas (i.e., spending function)
##' @param sigma standard deviation of sampling distribution
##' @param alpha confidence level
##' @return lower bound of FAB interval
##' @author Spencer Woody
##'
##' @export
theta_l_fab <- function(y, theta_vec, w_theta, sigma, alpha = 0.1) {
  theta_l_fab <- y + sigma * qnorm(alpha * (1 - w_theta)) - theta_vec
  
  lower_quantile_fab <- y + 
    sigma * qnorm(alpha * (1 - w_theta[which.min(abs(theta_l_fab))]))
  
  return(lower_quantile_fab)
}

##' Upper bound of FAB interval
##'
##' 
##' @title theta_u_fab
##' @param y single observation
##' @param theta_vec vector of thetas
##' @param w_theta vector of w's for thetas (i.e., spending function)
##' @param sigma standard deviation of sampling distribution
##' @param alpha confidence level
##' @return upper bound of FAB interval
##' @author Spencer Woody
##'
##' @export
theta_u_fab <- function(y, theta_vec, w_theta, sigma, alpha = 0.1) {
  
  theta_u_fab <- y + sigma * qnorm(1 - alpha * w_theta) - theta_vec
  
  upper_quantile_fab <- y + 
    sigma * qnorm(1 - alpha * w_theta[which.min(abs(theta_u_fab))])
  
  return(upper_quantile_fab)
}

