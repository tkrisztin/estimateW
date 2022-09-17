#' @exportS3Method
print.estimateW <- function(x, probs = c(.05,.5,.95),...) {
  if (is.null(x[["postr"]])) {
    model_type = "SLX"} else if (is.null(x[["Z"]])) {
      model_type = "SAR"} else {model_type = "SDM"}
  cat("Model type: ", model_type, "\n", sep = "")
  cat("Draws (total): ", x$param$nretain, " (", x$param$niter, ")\n", sep = "")
  coefs = rbind(x$postb,x$posts)
  if (!is.null(x[["postr"]])) {coefs = rbind(coefs,x$postr)}
  cat("Coefficients:\n")
  print.default(format(t(apply(coefs,c(1),stats::quantile,probs)), digits = 3L),
                print.gap = 2L, quote = FALSE)
  cat("\n")
  if (model_type != "SLX") {
    cat("Direct effect:\n")
    print.default(format(t(apply(x$post.direct,c(1),stats::quantile,probs)), digits = 3L),
                  print.gap = 2L, quote = FALSE)
    cat("Indirect effect:\n")
    print.default(format(t(apply(x$post.indirect,c(1),stats::quantile,probs)), digits = 3L),
                  print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  # invisible(x)
}


#' @exportS3Method
summary.estimateW <- function(object, ...) {
  coefs = rbind(object$postb,object$posts)
  if (!is.null(object[["postr"]])) {coefs = rbind(coefs,object$postr)}
  summary( t(coefs) )
}

#' @exportS3Method
summary.sim_dgp <- function(object, ...) {
  dat1 = data.frame(Y = object$Y,object$X)
  if (!is.null(object[["Z"]])) dat1 = cbind(dat1,Z = object$Z)
  summary( dat1 )
}


#' @exportS3Method
print.sim_dgp <- function(x, probs = c(.05,.5,.95),...) {
  coefs = x$para$posts
  if (!is.null(x$para[["beta1"]])) coefs = c(coefs,x$para$beta1)
  if (!is.null(x$para[["beta2"]])) coefs = c(coefs,x$para$beta2)
  if (!is.null(x$para[["beta3"]])) coefs = c(coefs,x$para$beta3)
  if (!is.null(x$para[["posts"]])) coefs = c(coefs,x$para$posts)
  cat("Simulated data generating process:\n")
  print.default(format(coefs, digits = 3L),
                print.gap = 2L, quote = FALSE)
  cat("\n")
  # invisible(x)
}
