


#' Efficient update of the log-determinant and the matrix inverse using the
#' Matrix Determinant Lemma and the Sherran-Morrison formula
#'
#' @param ch_ind The index of rows to update
#' @param diff A \eqn{ch_ind x n} matrix. This value will be added to the corresponding rows of \eqn{A}.
#' @param AI A \eqn{n x n} inverse of the matrix \eqn{A}, where \eqn{A = I - \rho W}.
#' @param logdet The log-determinant of \eqn{A}.
#'
#' @return A list containing \eqn{A^-1} and the log determinant of \eqn{A}
update_ldetAI <- function(ch_ind, diff, AI, logdet) {
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
