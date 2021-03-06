% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inRange.R
\name{inRange}
\alias{inRange}
\title{inRange}
\usage{
inRange(x, intervals)
}
\arguments{
\item{x}{the value to check if it's in the intervals}

\item{intervals}{a vector of even length; the odd elements
correspond to the lower end of the interval, and the even
elements are the upper ends of the intervals}
}
\value{
logical; TRUE if x is in any of intervals
}
\description{
Check whether x is in the range of any one of several disjoint
intervals
}
\author{
Spencer Woody
}
