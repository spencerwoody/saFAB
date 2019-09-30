// Copyight James G. Scott 2014
// Published under GPL3 licence

#include <Rcpp.h>


// [[Rcpp::export]]
List PredictiveRecursionFDR(NumericVector z, IntegerVector sweeporder,
    NumericVector grid_x, NumericVector theta_guess,
    double mu0 = 0.0, double sig0 = 1.0, double nullprob=0.95, double decay = -0.67) {
  // z: z statistics
  // sweeporder: an ordering of the points in z, usually 5-10 stacked permutations of 1 ... n
  // grid_x: a grid of points at which the alternative density will be approximated
  // theta_guess: an initial guess for the sub-density under the alternative hypothesis
  // nullprob: an initial guess for the fraction of null cases
  // decay: the stochastic-approximation decay parameter, should be in (-1, -2/3)

  // Set-up
  int n = sweeporder.size();
  int k, gridsize = grid_x.size();
  NumericVector theta_subdens(clone(theta_guess));
  double pi0 = nullprob;
  NumericVector joint1(gridsize);
  NumericVector ftheta1(gridsize);
  double m0, m1, mmix, cc;

  // Begin sweep through the data
  for(int i=0; i<n; i++) {
    k = sweeporder[i];
    if(i % 200 == 0) Rcpp::checkUserInterrupt();  
    cc = pow(3.0+(double)i, decay);
    joint1 = dnorm(grid_x, z[k]-mu0, sig0) * theta_subdens;
    m0 = pi0*R::dnorm(z[k] - mu0, 0.0, sig0, 0);
    m1 = trapezoid(grid_x, joint1);
    mmix = m0 + m1;
    pi0 = (1.0-cc)*pi0 + cc*(m0/mmix);
    ftheta1 = joint1/mmix;
    theta_subdens = (1.0-cc)*theta_subdens + cc*ftheta1;
  }

  // Now calculate marginal distribution along the grid points
  NumericVector y_mix(gridsize);
  NumericVector y_signal(gridsize);
  for(int j=0; j<gridsize; j++) {
    joint1 = dnorm(grid_x, grid_x[j] - mu0, sig0)*theta_subdens;
    m0 = pi0*R::dnorm(grid_x[j], mu0, sig0, 0);
    m1 = trapezoid(grid_x, joint1);
    y_mix[j] = m0 + m1;
    y_signal[j] = m1/(1.0-pi0);
  }

  return Rcpp::List::create(Rcpp::Named("grid_x")=grid_x,
          Rcpp::Named("theta_subdens")=theta_subdens,
          Rcpp::Named("pi0")=pi0,
          Rcpp::Named("y_mix")=y_mix,
          Rcpp::Named("y_signal")=y_signal
          );
}

