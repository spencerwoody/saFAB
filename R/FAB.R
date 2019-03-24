##' Create FAB intervals, given the spending function
##'
##' 
##' @title Create FAB intervals, given the spending function 
##' @param theta 
##' @param w 
##' @param sigma 
##' @param alpha 
##' @param yMin 
##' @param yMax 
##' @param yNum 
##' @param verbose 
##' @return 
##' @author Spencer Woody
##'
##' @export
FAB <- function(theta, w, sigma, alpha = 0.10, yMin = -5, yMax = 5,
                yNum = 5000, verbose = TRUE) {
  require(dplyr)
  require(rootSolve)
  require(splines)
  require(ggplot2)
  require(latex2exp)

  numTheta <- length(theta)

  # Create vectors for upper and lower acceptance regions
  Al <- rep(NA, numTheta)
  Au <- rep(NA, numTheta)

  if (verbose) cat("\nCreating acceptance regions...\n")

  # Fill create acceptance regions
  prog <- progress_estimated(numTheta)
  for (i in 1:numTheta) {

    # Tail regions
    wl <- alpha * w[i]
    wu <- 1 - alpha + wl

    # saFAB acceptance regions
    Al[i] <- qnorm(wl, theta[i], sigma)
    Au[i] <- qnorm(wu, theta[i], sigma)

    # Print progress update
    if (verbose) prog$tick()$print()

  }

  if (verbose) cat("\n")

  # Make data.frame for acceptance regions
  Adf <- data.frame(
    theta = rep(theta, 2),
    A = c(Al, Au),
    ul = c(rep("lower", numTheta), rep("upper", numTheta))
  )

  # Create spline functions for acceptance regions
  AlFun <- splinefun(theta, Al)
  AuFun <- splinefun(theta, Au)

  # Create confidence intervals along a sequence of y values
  yGrid <- seq(yMin, yMax, length.out = yNum)

  Cdf <- data.frame(
    y = yGrid,
    intervals = I(vector("list", length(yGrid))),
    numIntervals = NA,
    intervalLength = NA
  )

  if (verbose) cat("\nCreating confidence intervals...\n")
  prog <- progress_estimated(yNum)
  for (i in 1:length(yGrid)) {

    # Temporary auxillary functions
    Lfun <- function(x) {
      AlFun(x) - yGrid[i]
    }

    Ufun <- function(x) {
      AuFun(x) - yGrid[i]
    }


    # Roots of auxillary functions
    rootsL <- uniroot.all(Lfun, range(theta))
    rootsU <- uniroot.all(Ufun, range(theta))

    # Intervals come from roots of auxillary functions
    Cdf$intervals[[i]] <- c(rootsU, rootsL)

    # Number of disjoint intervals
    Cdf$numIntervals[i] <- length(Cdf$intervals[[i]]) / 2

    # Length of intervals
    Cdf$intervalLength[i] <- totalLength(Cdf$intervals[[i]])

    # Print progress bar
    if (verbose) prog$tick()$print()
  }

  if (verbose) cat("\n")

  # Output the confidence interval dataframe for plotting

  # Confidence interval dataframe for plotting
  CdfPlotting <- data.frame(
    y = rep(NA, sum(Cdf$numIntervals)),
    lower = NA,
    upper = NA
  )

  rowCount <- 1

  for (j in 1:nrow(Cdf)) {
    numIntervalsJ <- Cdf$numIntervals[j]

    Rows <- (rowCount):(rowCount + numIntervalsJ - 1)

    CdfPlotting$y[Rows] <- rep(Cdf$y[j], numIntervalsJ)

    CdfPlotting$lower[Rows] <- Cdf$intervals[[j]][seq(1, numIntervalsJ * 2, by = 2)]
    CdfPlotting$upper[Rows] <- Cdf$intervals[[j]][seq(2, numIntervalsJ * 2, by = 2)]

    rowCount <- rowCount + numIntervalsJ

  }

  CdfPlot <- ggplot(CdfPlotting) +
    geom_linerange(aes(y, ymin = lower, ymax = upper), alpha = 0.25) +
    labs(x = "y", y = TeX("$\\theta \\in C(y)$"))

  # Output
  return(list(
    Adf = Adf,
    Cdf = Cdf,
    CdfPlotting = CdfPlotting,
    CdfPlot = CdfPlot
  ))

}
