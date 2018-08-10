#' Create selection-adjusted FAB intervals, given the spending function
#'
#' @param theta vector of thetas from spending function
#' @param w vector for spending function
#' @param t truncation point
#' @param alpha confidence level
#' @param yMin minimal value of y for the confidence interval
#' @param yMax maxmimal value of y for the confidence interval
#' @param yNum length of vector of y's at which to make confidence interval
#' @param verbose logical; if TRUE (default), print progress bars
#'
#' @export

saFAB <- function(theta, w, t, sigma, alpha = 0.10,
                  yMin = -5, yMax = 5, yNum = 5000, verbose = TRUE) {
  require(dplyr)
  require(rootSolve)
  require(splines)

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
    Al[i] <- qtnorm(wl, theta[i], sigma, t)
    Au[i] <- qtnorm(wu, theta[i], sigma, t)

    # Print progress update
    if (verbose) prog$tick()$print()

  }

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

  # Output
  return(list(
    Adf = Adf,
    Cdf = Cdf,
    CdfPlotting =
  ))

}
