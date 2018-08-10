#' Calculate total length of disjoint intervals
#'
#' @param intervals vector containing start and end points of intervals
#'

totalLength <- function(intervals) {
  sum(intervals[seq(2, length(intervals), by = 2)] -
        intervals[seq(1, length(intervals), by = 2)])
}
