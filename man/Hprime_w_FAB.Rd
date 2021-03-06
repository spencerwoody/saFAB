% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hprime_w_FAB.R
\name{Hprime_w_FAB}
\alias{Hprime_w_FAB}
\title{Derivative of objective function to be minimized, as described by...}
\usage{
Hprime_w_FAB(w, theta, sigma, alpha, marginal_fun, ...)
}
\arguments{
\item{w}{weight corresponding to theta; i.e., w(theta)}

\item{theta}{the mean of sampling distribution y}

\item{sigma}{standard deviation of sampling distribution for y}

\item{alpha}{confidence level}

\item{marginal_fun}{the function for the marginal distribution}

\item{...}{optional arguments to the marginal function}

\item{t}{the truncation, i.e. y ~ N(theta, sigma) * 1(abs(y) > t)}
}
\description{
Derivative of objective function to be minimized, as described by...
}
