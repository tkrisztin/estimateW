% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbinompdf.R
\name{bbinompdf}
\alias{bbinompdf}
\title{The beta binomial pdf for sparsity priors}
\usage{
bbinompdf(x, nsize, a, b, min_k = 0, max_k = nsize)
}
\arguments{
\item{x}{Number of neighbors (scalar)}

\item{nsize}{Maximal number of elements}

\item{a}{Scalar prior parameter a}

\item{b}{Scalar prior parameter b}

\item{min_k}{Minimum number of elements (defaults to 0)}

\item{max_k}{Maximum number of elements (defaults to nsize + 1)}
}
\value{
Density of neighbors
}
\description{
The beta binomial pdf for sparsity priors
}
