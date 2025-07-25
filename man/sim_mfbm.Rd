% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{sim_mfbm}
\alias{sim_mfbm}
\title{Simulation of the mfBm}
\usage{
sim_mfbm(
  n = 50,
  H = c(0.2, 0.2),
  sig = c(1, 1),
  rho = matrix(c(1, 0.5, 0.5, 1), nc = 2),
  eta = matrix(c(0, 0.5, -0.5, 0), nc = 2),
  delta = 1,
  m = 1,
  gn = FALSE,
  plot = FALSE,
  print = TRUE,
  choix = NULL,
  forceEta = FALSE
)
}
\arguments{
\item{n}{Length of the time series.}

\item{H}{Vector of Hurst coefficients.}

\item{sig}{Vector of scale coefficients.}

\item{rho}{Symmetric matrix of instantaneous correlation coefficients.}

\item{eta}{Antisymmetric matrix of asymmetry coefficients.}

\item{delta}{Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years).}

\item{m}{Number of simulations.}

\item{gn}{Logical. If TRUE, the increments process is returned (multivariate fractional Gaussian noise).}
}
\value{
list of m numeric matrices of size  d x n
}
\description{
This function simulates a multivariate fractional Brownian motion of dimension \eqn{d}
using circulant embedding. The input parameters must satisfy the condition
for the existence of the covariance matrix, which is verified by the function.
}
\references{
Amblard, Pierre-Olivier, et al. "Basic properties of the multivariate fractional Brownian motion." arXiv preprint arXiv:1007.0828 (2010).
}
\author{
Prof. J.F. Coeurjolly
}
