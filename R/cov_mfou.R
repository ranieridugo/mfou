#' Sample cross-covariance estimator
#'
#' @description This function calculates the sample estimate of the cross-covariance from
#' observations:
#' - `cov_iik` calculates the sample autocovariance at lag k for a time series \eqn{y^i},
#' \deqn{\frac{1}{n-k}\sum_{t=1}^{n-k} y^i_{t+k} y^i_{t}-\mu_i^2;}
#' - `cov_ijk` calculates the sample cross-covariance between two time series of
#' the same length, \eqn{y^i} at time t + k and \eqn{y^j} at time t,
#' \deqn{\frac{1}{n-k}\sum_{t=1}^{n-k} y^i_{t+k} y^j_{t}-\mu_i\mu_j,}
#' \eqn{k\in\mathbb{Z}}.
#'
#' @param yi Numeric vector representing a time series.
#' @param yj Numeric vector of the same length of yi representing a time series.
#' @param k Lag (integer) at which to calculate the cross-covariance.
#' @param mui Mean of yi. If NULL, the sample average is used.
#' @param muj Mean of yj. If NULL, the sample average is used.
#'
#' @returns Numeric value.
#' @export
#' @rdname ccf_sample
#'
#' @examples
#' x = sin(1 : 100) + rnorm(1 : 100)
#' y = cos(1 : 100) + rnorm(1 : 100)
#' cov_iik(y, 0)
#' cov_ijk(x, y, 0)
#' cov_ijk(x, y, 1)
cov_iik <- function(yi, k, mui = NULL) {
  cov_ijk(yi, yi, k, mui = mui, muj = mui)
}

#' @export
#' @rdname ccf_sample
cov_ijk <- function(yi, yj, k, mui = NULL, muj = NULL) {
  if (length(yi) != length(yj)) stop("yi and yj must have same length.")
  y = cbind(yi, yj)
  n = nrow(y)
  if (is.null(mui) | is.null(muj)) {
    mui = mean(yi, na.rm = TRUE)
    muj = mean(yj, na.rm = TRUE)
  }
  if (k > 0) {
    mean(y[(1 + k) : n, 1] * y[1 : (n - k), 2] - mui * muj, na.rm = TRUE)
  } else if (k < 0) {
    k = - k
    mean(y[(1 + k) : n, 2] * y[1 : (n - k), 1] - mui * muj, na.rm = TRUE)
  } else if (k == 0) {
    mean(y[, 1] * y[, 2] - mui * muj, na.rm = TRUE)
  }
}

#' Variance function of the fOU process
#'
#' @description Variance of the fractional Ornstein-Uhlenbeck (fOU) process, based on model parameters,
#' \eqn{\text{Var}\left(Y_{t}\right)}.
#'
#' @param a Speed of mean reversion coefficient.
#' @param nu Diffusion coefficient.
#' @param H Hurst coefficient.
#'
#' @returns Numeric value.
#' @export
var_fou <- function(a, nu, H) {
  nu ^ 2 / (2 * a ^ (2 * H)) * gamma(1 + 2 * H)
}

#' Cross-covariance function of the mfOU process
#'
#' @description This function calculates the covariance function of the multivariate fractional
#' Ornstein-Uhlenbeck (mfOU) process, based on model parameters:
#' - `cov_iik_mfou` calculates the autocovariance at lag k of component i,
#' \deqn{\text{Cov}\left(Y_{t+k}^i,Y_t^i\right);}
#' - `cov_ij0_mfou` calculates the covariance between component i and j,
#' \deqn{\text{Cov}\left(Y_t^i,Y_t^j\right);}
#' - `cov_ijk_mfou` calculates the covariance between component i at time \eqn{t+k} and
#'  component j at time \eqn{t}, \deqn{\text{Cov}\left(Y_t^i,Y_t^j\right);}
#' - `cov_jik_mfou` calculates the covariance between component i at time \eqn{t} and
#'  component j at time \eqn{t+k}, \deqn{\text{Cov}\left(Y_t^i,Y_t^j\right);}
#'  \eqn{k\in\mathbb{R}}.
#'
#' @param ai Speed of mean reversion of component i.
#' @param aj Speed of mean reversion of component j.
#' @param nui Diffusion coefficient of component i.
#' @param nuj Diffusion coefficient of component j.
#' @param Hi Hurst coefficient of component i.
#' @param Hj Hurst coefficient of component j.
#' @param rhoij Correlation coefficient of the underlying mfBm.
#' @param etaij Asymmetry coefficient of the underlying mfBm.
#' @param k Lag, \eqn{k >= 0}.
#'
#' @returns Numeric value.
#' @export
#' @rdname cov_mfou
cov_iik_mfou <- function(ai, nui, Hi, k) {
  var = var_fou(ai, nui, Hi)
  if(k < 0) k = - k
  if (k == 0) {
    value = var
  } else if (k > 0) {
    integrand = function(y) exp(- abs(y)) * abs(k * ai + y) ^ (2 * Hi)
    value = nui ^ 2 / (2 * ai ^ (2 * Hi)) *
      (0.5 * integrate(f = integrand, lower = - Inf, upper = Inf)$value -
         abs(ai * k) ^ (2 * Hi))
  }
  return(unname(value))
}

#' @export
#' @rdname cov_mfou
cov_ij0_mfou <- function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij) {
  H = Hi + Hj
  if(length(H) > 1){
    value = 1/(2 * (ai + aj)) *
      gamma(H + 1) * nui * nuj * ((ai ^ (1 - H) + aj ^ (1 - H)) * rhoij +
                                    (aj ^ (1 - H) - ai ^ (1 - H)) * etaij)
  } else {
    if(H != 1){
      value = 1/(2 * (ai + aj)) *
        gamma(H + 1) * nui * nuj * ((ai ^ (1 - H) + aj ^ (1 - H)) * rhoij +
                                      (aj ^ (1 - H) - ai ^ (1 - H)) * etaij)
    } else if (H == 1) {
      value = ((nui * nuj)/(ai + aj)) * (rhoij + etaij * 0.5 *
                                           (log(aj) - log(ai)))
    }
  }
  return(unname(value))
}

#' @export
#' @rdname cov_mfou
cov_ijk_mfou <- function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij, k) {
  if (k < 0) stop("'k' must be eqnual or greater than zero")
  H = Hi + Hj
  cov = cov_ij0_mfou(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij)
  if(k == 0) {
    value = cov
  } else if (k > 0) {
    if(H != 1) {
      value = exp(- ai * k) * cov + nui * nuj * exp(- ai * k) *
        H * (H - 1) * (rhoij + etaij) * 0.5 * I_ij(ai, aj, Hi, Hj, k)
    } else if (H == 1){
      value = exp(- ai * k) * cov - nui * nuj * exp(- ai * k) *
        etaij * 0.5 * I_ij(ai, aj, Hi, Hj, k)
    }
  }
  return(unname(value))
}

#' @export
#' @rdname cov_mfou
cov_jik_mfou <- function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij, k) {
    if (k < 0) stop("'k' must be equal or greater than zero")
    H = Hi + Hj
    cov = cov_ij0_mfou(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij)
    if(k == 0) {
      value = cov
    } else if (k > 0) {
      if(H != 1) {
        value =
          exp(- aj * k) * cov + nui * nuj * exp(- aj * k) *
          H * (H - 1) * (rhoij - etaij) * 0.5 * I_ji(ai, aj, Hi, Hj, k)
      } else if (H == 1){
        value =
          exp(- aj * k) * cov + nui * nuj * exp(- aj * k) *
          etaij * 0.5 * I_ji(ai, aj, Hi, Hj, k)
      }
    }
    return(unname(value))
}

#' Autocovariance function of the fOU
#'
#' @description
#' These functions calculate the autocovariance function of the fractional
#' Ornstein-Uhlenbeck process for a number of lags, based on on model parameters.
#' - `acf_fou` calculates the exact autocovariances;
#' - `acf_asy_fou`calculates the asymptotic autocovariances arising in the regime of
#' small speed of mean reversion.
#'
#'
#' @param a Speed of mean reversion coefficient.
#' @param nu Diffusion coefficient.
#' @param H Hurst coefficient.
#' @param lag_max Integer. The autocovariance is calculated for all integers between
#' \eqn{0} and `lag_max`.
#' @param delta Time step in the discrete representation of the process (e.g., 1/252 for daily observations if time is in years).
#'
#' @returns Numeric vector.
#' @export
acf_fou =
  function(a, nu, H, lag_max = 50, delta = 1/252){
    acf =
      sapply(X = c(0 : lag_max) * delta, FUN = cov_iik_mfou,
             ai = a, nui = nu, Hi = H)|> unname()
    names(acf) <- c(0 : lag_max)
    return(acf)
  }

acf_asy_fou =
  function(nu, H, var, lag_max = 50, delta = 1/252){
    acf =
      sapply(X = c(0 : lag_max) * delta, FUN = cov_asy_iik_mfou,
             nui = nu, Hi = H, vari = var)|> unname()
    names(acf) <- c(0 : lag_max)
    return(acf)
  }

#' Cross-covariance function of the mfOU
#'
#' @description
#' These functions calculate the cross-covariance function of the multivariate
#' fractional Ornstein-Uhlenbeck process for a number of lags, based on on model parameters.
#' - `ccf_mfou` calculates the exact cross-covariances;
#' - `ccf_asy_mfou`calculates the asymptotic cross-covariances arising in the regime of
#' small speed of mean reversion coefficients.
#'
#' @param ai Speed of mean reversion of component i.
#' @param aj Speed of mean reversion of component j.
#' @param nui Diffusion coefficient of component i.
#' @param nuj Diffusion coefficient of component j.
#' @param Hi Hurst coefficient of component i.
#' @param Hj Hurst coefficient of component j.
#' @param rhoij Correlation coefficient of the underlying mfBm.
#' @param etaij Asymmetry coefficient of the underlying mfBm.
#' @param lag_max Integer. The cross-covariance is calculated for all integers between
#' \eqn{0} and `lag_max`.
#' @param delta Time step in the discrete representation of the process (e.g., 1/252 for daily observations if time is in years).
#'
#' @returns Numeric value.
#' @export
ccf_mfou <-
  function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij, lag_max = 50, delta = 1/252) {
    ccf =
      c(
        sapply(X = c(lag_max : 1) * delta,
               FUN = cov_jik_mfou,
               ai = ai, aj = aj,
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij),
        cov_ij0_mfou(ai = ai, aj = aj,
                   nui =  nui, nuj = nuj,
                   Hi = Hi, Hj = Hj,
                   rhoij = rhoij, etaij = etaij),
        sapply(X = c(1 : lag_max) * delta,
               FUN = cov_ijk_mfou,
               ai = ai, aj = aj,
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij)
      )|> unname()
    names(ccf) <- c(- c(lag_max : 1), 0, c(1 : lag_max))
    return(ccf)
  }

ccf_asy_mfou <-
  function(nui, nuj, Hi, Hj, rhoij, etaij, cov0ij, lag_max = 50, delta = 1/252) {
    ccf =
      c(
        sapply(X = c(lag_max : 1) * delta,
               FUN = cov_asy_jik_mfou,
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij, cov0ij = cov0ij),
        cov0ij,
        sapply(X = c(1 : lag_max) * delta,
               FUN = cov_asy_ijk_mfou,
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij, cov0ij = cov0ij)
      )|> unname()
    names(ccf) <- c(- c(lag_max : 1), 0, c(1 : lag_max))
    return(ccf)
  }

#' Integral appearing in the cross-covariance function
#'
#' @param ai Speed of mean reversion of component i, \eqn{\alpha_i}.
#' @param aj Speed of mean reversion of component j, \eqn{\alpha_j}.
#' @param Hi Hurst coefficient of component i, \eqn{H_i}.
#' @param Hj Hurst coefficient of component j, \eqn{H_j}.
#' @param k Lag, \eqn{k > 0}.
#'
#' @returns Numeric
#' @export
#' @rdname i_ccf
I_ij <- function(ai, aj, Hi, Hj, k) {
  H = Hi + Hj
  integrand1 = function(u){u^(H - 2) * exp(- aj * u)}
  integrand2 = function(u){u^(H - 2) * (exp(ai * u) - exp(- aj * u))}
  value = 1 / (ai + aj) * (
    (exp((ai + aj) * k) - 1) *
      integrate(f = integrand1, lower = k, upper = Inf)[1]$value +
      integrate(f = integrand2, lower = 0 , upper = k)[1]$value
  )
  return(value)
}

#' @export
#' @rdname i_ccf
I_ji <- function (ai, aj, Hi, Hj, k) {
  H = Hi + Hj
  integrand1 = function(u){u^(H - 2) * exp(- ai * u)}
  integrand2 = function(u){u^(H - 2) * (exp(aj * u) - exp(- ai * u))}
  value = 1/(ai + aj) * (
    (exp((ai + aj) * k) - 1) *
      integrate(f = integrand1, lower = k, upper = Inf)[1]$value +
      integrate(f = integrand2, lower = 0 , upper = k)[1]$value
  )
  return(value)
}

#' Asymptotic cross-covariance function of the mfOU process
#'
#' @description Calculation of the asymptotic \eqn{\left(\alpha\to 0\right)} cross-covariance of
#'  the multivariate fractional Ornstein-Uhlenbeck (mfOU) process,
#'  based on model parameters:
#'  \deqn{\text{Cov}\left(Y_{t+k}^i,Y_t^j\right)=\text{Cov}\left(Y_{t}^i,Y_t^j\right)+
#'  \frac{\nu_i\nu_j}{2}\left(\rho_{i,j}+\eta_{i,j}\right)k^{H_i+H_j}+o(1),}
#'  \eqn{k\in\mathbb{R}}.
#'  - `cov_asy_iik_mfou` calculates the autocovariance at lag k;
#'  - `cov_asy_ijk_mfou` calculates the covariance between \eqn{Y_{t+k}^i} and \eqn{Y_t^j};
#'  - `cov_asy_jik_mfou` calculates the covariance between \eqn{Y_t^i} and \eqn{Y_{t+k}^j}.
#'
#' @param nui Diffusion coefficient of component i, \eqn{\nu_i}.
#' @param nuj Diffusion coefficient of component j, \eqn{\nu_j}.
#' @param Hi Hurst coefficient of component i, \eqn{H_i}.
#' @param Hj Hurst coefficient of component j, \eqn{H_j}.
#' @param rhoij Correlation coefficient of the underlying mfBm, \eqn{\rho_{ij}}.
#' @param etaij Asymmetry coefficient of the underlying mfBm, \eqn{\eta_{ij}}.
#' @param k Lag, \eqn{k > 0}.
#' @param vari Variance of component i.
#' @param cov0ij Covariance between component i and j.
#'
#' @returns Numeric
#' @export
#' @rdname cov_asy_mfou
cov_asy_iik_mfou <- function(nui, Hi, vari, k) {
  value =
    vari - 0.5 * nui ^ 2 * k ^ (2 * Hi)
  return(value)
}

#' @export
#' @rdname cov_asy_mfou
cov_asy_ijk_mfou <- function(nui, nuj, Hi, Hj, rhoij, etaij, cov0ij, k) {
  value =
    cov0ij - 0.5 * nui * nuj * (rhoij + etaij) * k ^ (Hi + Hj)
  return(value)
}

#' @export
#' @rdname cov_asy_mfou
cov_asy_jik_mfou <- function(nui, nuj, Hi, Hj, rhoij, etaij, cov0ij, k) {
  value =
    cov0ij - 0.5 * nui * nuj * (rhoij - etaij) * k ^ (Hi + Hj)
  return(value)
}
