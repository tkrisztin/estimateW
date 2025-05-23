% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6_samplers.R
\docType{class}
\name{beta_sampler}
\alias{beta_sampler}
\title{An R6 class for sampling slope parameters}
\format{
An \code{\link[R6]{R6Class}} generator object
}
\description{
This class samples slope parameters with a Gaussian prior from the conditional posterior.
Use the \link{beta_priors} class for setup.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{beta_prior}}{The current \code{\link{beta_priors}}}

\item{\code{curr_beta}}{The current value of \eqn{\beta}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-beta_sampler-new}{\code{beta_sampler$new()}}
\item \href{#method-beta_sampler-sample}{\code{beta_sampler$sample()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-beta_sampler-new"></a>}}
\if{latex}{\out{\hypertarget{method-beta_sampler-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{beta_sampler$new(beta_prior)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{beta_prior}}{The list returned by \code{\link{beta_priors}}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-beta_sampler-sample"></a>}}
\if{latex}{\out{\hypertarget{method-beta_sampler-sample}{}}}
\subsection{Method \code{sample()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{beta_sampler$sample(Y, X, curr_sigma)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{Y}}{The \eqn{N} by \eqn{1} matrix of responses}

\item{\code{X}}{The \eqn{N} by \eqn{k} design matrix}

\item{\code{curr_sigma}}{The variance parameter \eqn{\sigma^2}}
}
\if{html}{\out{</div>}}
}
}
}
