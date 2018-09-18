
##' Check whether x is in the range of any one of several disjoint
##' intervals
##'
##' 
##' @title inRange
##' @param x the value to check if it's in the intervals
##' @param intervals a vector of even length; the odd elements
##'   correspond to the lower end of the interval, and the even
##'   elements are the upper ends of the intervals
##' @return logical; TRUE if x is in any of intervals
##' @author Spencer Woody
##'
##' @export
inRange <- function(x, intervals) {
  any(x >= intervals[seq(1, length(intervals), by = 2)] &
        x <= intervals[seq(2, length(intervals), by = 2)])
}
