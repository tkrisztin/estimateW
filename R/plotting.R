


#' Graphical Summary of the estimated W matrix
#'
#' @param x \code{estimateW} object.
#' @param \dots further arguments are passed on to the invoked
#'
#' @import plot.matrix
#' @export
plot.estimateW = function(x,...) {
  w_pip <- apply(x$postw>0,c(1,2),mean)
  #colnames(w_pip) <- rownames(w_pip) <- shp$ISO3
  #w_pip2 <- w_pip[order(shp$LON,decreasing = FALSE),order(shp$LON,decreasing = FALSE)]
  #apply(x$postb,MARGIN = 1,FUN=quantile, c(0.05,0.5,0.95))

  upper_graph_pip=0.75
  graphics::par(mar=c(3.0, 3.0, 3.0, 3.0))
  plot(w_pip,
       col=c("white","lightgrey","black"),
       breaks=c(0,0.5,upper_graph_pip,1),
       xlab="",ylab="",main="",border=NA,font=3,cex.axis=0.75,las=2)
  graphics::abline(nrow(x$X)+1,-1)
}
