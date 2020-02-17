
##' @export
ptnorm1 <- function(q, mean = 0, sd = 1, t = 0, lower.tail = TRUE,
                    log.p = FALSE) {

  C <- 1 - pnorm((t - mean) / sigma)

  (pnorm((y - mean) / sd) - pnorm((t - mean) / sd)) / C * ifelse(y > t, 1, 0)

}
