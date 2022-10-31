


#' Graphical Summary of the estimated W matrix
#'
#' @param x \code{estimateW} object.
#' @param cols Main colors to use for the plot
#' @param breaks Breaks for the colors
#' @param \dots further arguments are passed on to be invoked
#'
#' @import plot.matrix
#' @export
plot.estimateW = function(x,
                          cols = c("white","lightgrey","black"),
                          breaks=c(0,0.5,0.75,1),
                          ...) {
  w_pip <- apply(x$postw>0,c(1,2),mean)
  graphics::par(mar=c(3.0, 3.0, 3.0, 3.0))
  plot(w_pip,
       col=cols,main="Posterior incl. prob. of W",
       breaks=breaks,border = NA,...)
  graphics::abline(nrow(x$postw)+1,-1)
}

#' Graphical summary of a generated spatial weight matrix
#'
#' @param x \code{sim_dgp} object
#' @param \dots further arguments are passed on to the invoked
#'
#' @import plot.matrix
#' @export
plot.sim_dgp = function(x, ...) {
  graphics::par(mar=c(3.0, 3.0, 3.0, 3.0))
  W = as.matrix(x$W); W[W>0] = 1;
  plot(W,
       col=c("white","black"),main = "Spatial weight matrix",
       breaks=c(0,0.5,1),
       border = NA,...)
  graphics::abline(nrow(W)+1,-1)
}

#' @exportS3Method
plot.exoW <- function(x, ...) {
  plot(c(x$postr),type="l", main = "Rho posterior draws")
}

#' @exportS3Method
plot.normalgamma <- function(x, ...) {
  plot(c(x$posts),type="l", main = "Sigma posterior draws")
}
