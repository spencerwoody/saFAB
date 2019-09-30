
// Copyight James G. Scott 2014
// Published under GPL3 licence

#include <Rcpp.h>


// [[Rcpp::export]]
List eval_pr_dens(NumericVector z, double mu0, NumericVector sig0, NumericVector grid_x, NumericVector grid_theta) {
  // z: z statistics
  // mu0: mean of Gaussian
  // sig0: vector of standard errors, se(z[i])
  // grid_x: a grid of points at which the alternative density pi(theta) is evaluated
  // grid_theta: the (unnormalized) mixing density pi(theta) at each point in grid_x
  // this function will evaluate the predictive density of each z point

  // Set-up
  int n = z.size();
  int gridsize = grid_x.size();
  NumericVector fsignal_z(n);
  NumericVector joint1(gridsize);
  double norm_constant = trapezoid(grid_x, grid_theta);

  // Begin sweep through the data, each time integrating by trap rule
  for(int i=0; i<n; i++) {
    if(i % 200 == 0) Rcpp::checkUserInterrupt();  
    joint1 = dnorm(grid_x, z[i] - mu0, sig0[i]) * grid_theta;
    fsignal_z[i] = trapezoid(grid_x, joint1)/norm_constant;
  }
  return Rcpp::List::create(Rcpp::Named("fsignal_z")=fsignal_z
          );
}

