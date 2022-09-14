


#' Graphical Summary of the estimated W matrix
#'
#' @param x \code{estimateW} object.
#' @param cols Main colors to use for the plot
#' @param breaks Breaks for the colours
#' @param main Legend title
#' @param \dots further arguments are passed on to the invoked
#'
#' @import plot.matrix
#' @export
plot.estimateW = function(x,
                          cols = c("white","lightgrey","black"),
                          breaks=c(0,0.5,0.75,1),
                          main="Posterior incl. prob. of W",
                          ...) {
  w_pip <- apply(x$postw>0,c(1,2),mean)
  graphics::par(mar=c(3.0, 3.0, 3.0, 3.0))
  plot(w_pip,
       col=cols,main = main,
       breaks=breaks,border = NA,...)
  graphics::abline(nrow(x$postw)+1,-1)
}

#' Graphical summary of a generated spatial weight matrix
#'
#' @param x \code{sim_sdmw} object
#' @param cols Main colors to use for the plot
#' @param breaks Breaks for the colours
#' @param main Legend title
#' @param \dots further arguments are passed on to the invoked
#'
#' @import plot.matrix
#' @export
plot.sim_sdmw = function(x,
                          cols = c("white","lightgrey","black"),
                          breaks=c(0,0.5,0.75,1),
                          main="Spatial weight matrix",
                          ...) {
  graphics::par(mar=c(3.0, 3.0, 3.0, 3.0))
  plot(x$W,
       col=cols,main = main,
       breaks=breaks,border = NA,...)
  graphics::abline(nrow(x$W)+1,-1)
}
