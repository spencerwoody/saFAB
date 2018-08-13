##' Sum length of disjoint intervals
##'
##' 
##' @title totalLength
##' @param intervals 
##' @return Sum of lengths of disjoint intervals 
##' @author Spencer Woody
##'
##' @export
totalLength <- function(intervals) {
  sum(intervals[seq(2, length(intervals), by = 2)] -
        intervals[seq(1, length(intervals), by = 2)])
}
