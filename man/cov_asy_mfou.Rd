% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_mfou.R
\name{cov_asy_iik_mfou}
\alias{cov_asy_iik_mfou}
\alias{cov_asy_ijk_mfou}
\alias{cov_asy_jik_mfou}
\title{Asymptotic cross-covariance function of the mfOU process}
\usage{
cov_asy_iik_mfou(nui, Hi, vari, k)

cov_asy_ijk_mfou(nui, nuj, Hi, Hj, rhoij, etaij, cov0ij, k)

cov_asy_jik_mfou(nui, nuj, Hi, Hj, rhoij, etaij, cov0ij, k)
}
\arguments{
\item{nui}{Diffusion coefficient of component i, \eqn{\nu_i}.}

\item{Hi}{Hurst coefficient of component i, \eqn{H_i}.}

\item{vari}{Variance of component i.}

\item{k}{Lag, \eqn{k > 0}.}

\item{nuj}{Diffusion coefficient of component j, \eqn{\nu_j}.}

\item{Hj}{Hurst coefficient of component j, \eqn{H_j}.}

\item{rhoij}{Correlation coefficient of the underlying mfBm, \eqn{\rho_{ij}}.}

\item{etaij}{Asymmetry coefficient of the underlying mfBm, \eqn{\eta_{ij}}.}

\item{cov0ij}{Covariance between component i and j.}
}
\value{
Numeric
}
\description{
Calculation of the asymptotic \eqn{\left(\alpha\to 0\right)} cross-covariance of
the multivariate fractional Ornstein-Uhlenbeck (mfOU) process,
based on model parameters:
\deqn{\text{Cov}\left(Y_{t+k}^i,Y_t^j\right)=\text{Cov}\left(Y_{t}^i,Y_t^j\right)+
 \frac{\nu_i\nu_j}{2}\left(\rho_{i,j}+\eta_{i,j}\right)k^{H_i+H_j}+o(1),}
\eqn{k\in\mathbb{R}}.
\itemize{
\item \code{cov_asy_iik_mfou} calculates the autocovariance at lag k;
\item \code{cov_asy_ijk_mfou} calculates the covariance between \eqn{Y_{t+k}^i} and \eqn{Y_t^j};
\item \code{cov_asy_jik_mfou} calculates the covariance between \eqn{Y_t^i} and \eqn{Y_{t+k}^j}.
}
}
