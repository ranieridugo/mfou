% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_mfbm.R
\name{cov_ijts_mfbm}
\alias{cov_ijts_mfbm}
\title{Cross-covariance function of the mfBm}
\usage{
cov_ijts_mfbm(ti, tj, Hi, Hj, rhoij, etaij, sigmai = 1, sigmaj = 1)
}
\arguments{
\item{ti}{Time of component i;}

\item{tj}{Time of component j;}

\item{Hi}{Hurst coefficient of component i.}

\item{Hj}{Hurst coefficient of component j.}

\item{rhoij}{Correlation coefficient (covariance at \eqn{t=1}).}

\item{etaij}{Asymmetry coefficient.}

\item{sigmai}{Scale coefficient of component i (square root of variance at \eqn{t=1}).}

\item{sigmaj}{Scale coefficient of component j (square root of variance at \eqn{t=1}).}
}
\value{
Numeric value.
}
\description{
Cross-covariance function of the multivariate fractional Brownian motion,
based on model parameters:
\deqn{\text{Cov}\left(B_{t_i}^{H_i},B_{t_j}^{H_j}\right).}
}
\references{
Amblard, Pierre-Olivier, et al. "Basic properties of the multivariate fractional Brownian motion." arXiv preprint arXiv:1007.0828 (2010).
}
