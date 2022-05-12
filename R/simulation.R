


#' Simulating a SAR data generating process
#'
#' @param n Number of spatial observations
#' @param tt Number of time observations
#' @param k Number of explanatory variables
#' @param rho The true rho parameter
#' @param beta The \eqn{k x 1} vector of true betas
#' @param sigma2 The true \eqn{\sigma^2} parameter
#' @param n_neighbor Number of neighbours for the \eqn{W} matrix
#' @param SYMMETRIC Should the generated \eqn{W} matrix be symmetric? (default: FALSE)
#'
#' @return A list with the generated \eqn{X}, \eqn{Y} and \eqn{W} and a list of parameters.
#' @export sarwsim
#'
#' @examples
#' dgp_dat = sarwsim(n =20, tt = 10, k=3, rho = .5, beta = c(1,-1,1), sigma2 = .5)
#'
sarwsim = function(n, tt, k, rho, beta, sigma2 = .5, n_neighbor = 7, SYMMETRIC = FALSE) {
  # Randomly generated spatial patterns
  xy <- cbind(stats::runif(n),stats::runif(n))
  W<-getWknn(xy,n_neighbor)

  # one-forward, one-behind W
  # W = matrix(0,n,n)
  # W[-1,-n] = W[-1,-n] + diag(n-1)
  # W[-n,-1] = W[-n,-1] + diag(n-1)
  # W = W + W%*%W + W%*%W%*%W
  # W[W>0] = 1; diag(W) = 0
  # W = W/rowSums(W)

  if (SYMMETRIC) {
    W = t(W) + W
    W = as.matrix(W/rowSums(W))
    W[W>0] = 1
  }
  diag(W) = 0

  A <- as.matrix(diag(n) - rho*W)
  Ainv <- solve(A)
  X <- cbind(1,matrix(stats::rnorm( n*(k-1) * tt),n*tt,k-1))
  Y = kronecker(Matrix::.sparseDiagonal(tt),Ainv) %*% (
    X %*% beta + stats::rnorm(n*tt, mean = 0,sd = sqrt(sigma2)))

  ret <- list(Y = as.matrix(Y),
              X = as.matrix(X),
              W = W,
              para = list(rho = rho,
                          beta = beta,
                          sigma2 = sigma2,
                          rho = rho,
                          xy = xy,
                          n_neighbor = n_neighbor))
  return(ret)
}

#' Constructs a nearest neighbour spatial weight matrix using lat/long coordinates and geodesic distance.
#'
#' @param xy An \eqn{n x 2} matrix of x and y coordinates.
#' @param k Number of nearest neighbours
#'
#' @return A sparse \eqn{n x n} spatial weight matrix
getWknn<-function(xy,k) {
  #construct Wmatrix
  n<- NROW(xy)
  nn<-FNN::knn.index(xy,k=k)
  W<-Matrix::Matrix(0,n,n)
  cseq = c(1:n)
  for (i in 1:k) {
    tmp1<-cbind( cseq,nn[,i],1/k)
    tmp2<-cbind( cseq,cseq,0)
    tmp3<-rbind(tmp1,tmp2)
    if (i==1) {
      W<-Matrix::sparseMatrix(i=tmp3[,1],j=tmp3[,2],x=tmp3[,3])
    } else {
      W <- W + Matrix::sparseMatrix(i=tmp3[,1],j=tmp3[,2],x=tmp3[,3])
    }
  }
  diag(W) = 0
  return(W)
}
