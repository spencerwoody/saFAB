% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rss_nonzero.R
\name{rss_nonzero}
\alias{rss_nonzero}
\title{rss}
\usage{
rss_nonzero(n, p, mu, tau)
}
\arguments{
\item{n}{number of samples to generate}

\item{p}{proportion of zeros}

\item{tau}{standard deviation of normal (slab) prior}
}
\value{
n samples from the spike-and-slab prior
}
\description{
Generate n samples from the spike-slab normal prior, a mixture of
a normal slab and point mass spike at zero
}
\details{
theta ~ p * 1(theta = 0) + (1 - p) * N(theta; 0, tau^2)
}
\author{
Spencer Woody
}
