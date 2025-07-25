% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_est_gmm.R
\name{par_list2vec}
\alias{par_list2vec}
\alias{par_vec2list}
\title{Converting parameters between list or vector form}
\usage{
par_list2vec(lpar, d, type = "exa")

par_vec2list(vpar, d, type = "exa")
}
\arguments{
\item{lpar}{List of parameters.}

\item{d}{Dimensionality of the system of time-series.}

\item{type}{One among \code{exa}, \code{asy}, \code{cau} to determine which version of the model the user wants to estimate
among exact, asymptotic, and causal.}

\item{vpar}{Vector of parameters.}
}
\value{
Numeric vector or List of numeric objects.
}
\description{
Functions that convert the way in which parameters are stored, from vector to list
or vice-verse.
}
