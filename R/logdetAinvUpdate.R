


#' Efficient update of the log-determinant and the matrix inverse
#'
#' While updating the elements of the spatial weight matrix in SAR and SDM type spatial models in a
#' MCMC sampler, the log-determinant term has to be regularly updated, too.
#' When the binary elements of the adjacency matrix are treated unknown, the Matrix Determinant Lemma
#' and the Sherman-Morrison formula are used for computationally efficient updates.
#'
#'
#' Let \eqn{A = (I_n - \rho W)} be an invertible \eqn{n} by \eqn{n} matrix. \eqn{v} is an \eqn{n} by \eqn{1}
#' column vector of real numbers and \eqn{u} is a binary vector containing a single one and zeros otherwise.
#' Then the Matrix Determinant Lemma states that:
#'
#' \deqn{A + uv' = (1 + v'A^{-1}u)det(A)}.
#'
#' This provides an update to the determinant, but the inverse of \eqn{A} has to be updated as well.
#' The Sherman-Morrison formula proves useful:
#'
#' \deqn{(A + uv')^{-1} = A^{-1} \frac{A^{-1}uv'A^{-1}}{1 + v'A^{-1}u}}.
#'
#' Using these two formulas, an efficient update of the spatial projection matrix determinant can be achieved.
#'
#' @param ch_ind vector of non-negative integers, between 1 and \eqn{n}. Denotes which rows of \eqn{A}
#'  should be updated.
#' @param diff a numeric \code{length(ch_ind)} by \code{n} matrix. This value will be added to the corresponding rows of \eqn{A}.
#' @param AI numeric \eqn{n} by \eqn{n} matrix that is the inverse of \eqn{A = (I_n - \rho W)}. This inverse will
#' be updated using the Sherman-Morrison formula.
#' @param logdet single number that is the log-determinant of the matrix \eqn{A}. This log-determinant
#' will be updated through the Matrix Determinant Lemma.
#'
#' @return A list containing the updated \eqn{n} by \eqn{n} matrix \eqn{A^{-1}}, as well as the
#' updated log determinant of \eqn{A}
#'
#' @export
#'
#' @references
#' Sherman, J., and Morrison, W. J. (1950) Adjustment of an inverse matrix corresponding to a
#' change in one element of a given matrix. \emph{The Annals of Mathematical Statistics}, \bold{21(1)},
#' 124-127.
#'
#' Harville, D. A. (1998) Matrix algebra from a statistician's perspective. Taylor & Francis.
logdetAinvUpdate <- function(ch_ind, diff, AI, logdet) {
  for (i in 1:length(ch_ind)) {
    ii <- ch_ind[i]
    diffAI <- t(diff[i, ]) %*% AI
    if (diffAI[ii] <= -1) {
      return(list(AI = AI, logdet = NA))
    } else {
      logdet <- log(1 + diffAI[ii]) + logdet
    }
    AI <- AI - AI[, ii] %*% diffAI / as.double(1 + diffAI[ii])
  }
  return(list(AI = AI, logdet = logdet))
}
