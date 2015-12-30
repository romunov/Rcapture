#' Merge results from profileCI function
#
#' To summarize results more efficients, use this function to merge them into one object.
#'
#' @param x List. A list of \code{closedp} objects to be merged into a single \code{matrix}.
#' 
#' @author Roman Lustrik \email{roman.lustrik@@biolitika.si}
#' @export

rbindCIs <- function(x) {
  list.of.res <- lapply(x, "[[", "results")
  num.samp.units <- lapply(x, "[[", "n")
  list.of.res <- mapply(list.of.res, num.samp.units, FUN = function(x1, x2) {
    cbind(x1, sample.size = x2)
  }, SIMPLIFY = FALSE)
  do.call("rbind", list.of.res)
}
