#' Forecasting the mfBm using Gaussian conditioning
#'
#' @description
#' Using a standard result for Gaussian random vectors, these functions allow
#' us to compute point forecasts and associated standard errors in the framework
#' of the multivariate fractional Brownian motion of dimension \eqn{d}.
#'
#' Given that the vector that we want to predict, \eqn{X:=B_{(t+h)\Delta}\in\mathbb{R}^d}, and the
#' history of information, \eqn{Y:=\left(B_{t\Delta},\dots, B_{\Delta}\right)\in\mathbb{R}^{d\times t}},
#' are jointly Gaussian, we use that:
#'
#' \deqn{\mu_{X|Y}=\mu_X+\Sigma_{XY}\Sigma_{YY}^{-1}\left(Y-\mu_Y\right),}
#' \deqn{\Sigma_{XX|Y}=\Sigma_{XX}-\Sigma_{XY}\Sigma_{YY}^{-1}\Sigma_{XY}^T,}
#' where we used the notation \eqn{\mu_Z=\mathbb{E}(Z)} and
#' \eqn{\Sigma_{ZU} = \mathbb{E}((Z-\mu_Z)(U-\mu_U)^T)}.
#'
#' - `cov_jhdt_mfbm` calculates the matrix \eqn{\Sigma_{XY}} described above;
#' - `cov_dt_mfbm` calculates the matrix \eqn{\Sigma_{YY}} described above;
#' - `cov_0_mfbm` calculates the matrix \eqn{\Sigma_{XX}} described above;
#' - `fcast_point_mfbm` calculates point forecasts, i.e. the vector \eqn{\mu_{X|Y}} described above;
#' - `fcast_se_mfbm` calculates standard errors of the forecasts, i.e. the matrix \eqn{\Sigma_{XX|Y}} described above;
#'
#' The functions `fcast_point_mfbm` and `fcast_se_mfbm` are not efficient in
#' rolling window applications with fixed parameters. In which case,
#' efficiency gains could often be obtained by avoiding redundant matrix operations.
#'
#' @param n Integer. Length of the time series.
#' @param h Integer. Forecast horizon (e.g. 1, 5, 10 steps ahead).
#' @param H Vector of Hurst coefficients.
#' @param sig Vector of scale coefficients.
#' @param rho Symmetric matrix of instantaneous correlation coefficients.
#' @param eta Antisymmetric matrix of asymmetry coefficients.
#' @param mu Vector of means. Default is NULL, in which case it is computed as sample average.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years).
#' @param my Numeric matrix containing the time series in its columns.
#' @param j Component of the system. Integer between 1 and \eqn{d} or column name.
#' @param plot Logical. Default is FALSE. If TRUE, a panel of plots will exhibit
#' the weights of the linear combination.
#' @returns Numeric matrices.
#' @export
#' @rdname fcast_mfbm
#' @references
#' Bibinger, Markus, Jun Yu, and Chen Zhang. "Modeling and Forecasting Realized Volatility with Multivariate Fractional Brownian Motion." arXiv preprint arXiv:2504.15985 (2025).
#'
#' https://en.wikipedia.org/wiki/Multivariate_normal_distribution
cov_jhdt_mfbm <- function(n, h, H, sig, rho, eta, j, delta) {
  d = length(H)
  if(!all(length(H), nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  sxy = numeric()
  for (m in 1 : n) {
    for (l in 1 : d) {
      tmp = cov_ijts_mfbm(ti = (n + h) * delta,
                          tj = m * delta,
                          Hi = H[j],
                          Hj = H[l],
                          rhoij = rho[j ,l],
                          etaij = eta[j, l],
                          sigmai = sig[j],
                          sigmaj = sig[l])
      sxy = c(sxy, tmp)
    }
  }
  return(sxy)
}

#' @export
#' @rdname fcast_mfbm
cov_dt_mfbm <- function(n, H, sig, rho, eta, delta) {
  d = length(H)
  if(!all(length(H), nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
      ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  sigma_yy = matrix(NA, nrow = d * n, ncol = d * n)
  for (i in (1 : n)) {
    for (l in (1 : d)) {
      foo = numeric()
      for (s in (i : n)) {
        for (j in (ifelse(s == i, l, 1) : d)) {
          tmp = cov_ijts_mfbm(ti = i * delta,
                              tj = s * delta,
                              Hi = H[l],
                              Hj = H[j],
                              rhoij = rho[l ,j],
                              etaij = eta[l, j],
                              sigmai = sig[l],
                              sigmaj = sig[j])
          foo = c(foo, tmp)
        }
      }
      n_foo = (d - l + 1 + d * (n - i))
      sigma_yy[(i - 1) * d + l, (d * n - n_foo + 1) : (d * n)] = foo # CHECK
    }
  }
  sigma_yy[lower.tri(sigma_yy)] <- t(sigma_yy)[lower.tri(sigma_yy)]
  return(sigma_yy)
}

#' @export
#' @rdname fcast_mfbm
cov_0_mfbm <- function(n, H, sig, rho, eta, delta) {
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  cov = matrix(NA, nc = d, nr = d)
  for (i in 1 : d){
    for (j in i : d) {
      cov[i, j] = cov_ijts_mfbm(n * delta, n * delta, H[i], H[j],
                                rho[i, j], eta[i, j],
                                sigmai = sig[i], sigmaj = sig[j])
    }
  }
  cov[lower.tri(cov)] = t(cov)[lower.tri(cov)]
  return(cov)
}

#' @export
#' @rdname fcast_mfbm
fcast_point_mfbm <- function(my, h, H, sig, rho, eta, delta,
                             mu = NULL, plot = FALSE) {
  if(sum(is.na(my)) > 0) stop("Error: your time series contain NAs.")
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions doesn't match.")
  fcast = rep(NA, d)
  n = nrow(my)
  if(is.null(mu)) mu = apply(my, 2, mean)
  syy = cov_dt_mfbm(n, H, sig, rho, eta, delta)
  syy_inv = solve(syy)
  if(plot == TRUE) par(mfrow = c(d, d), mar = c(2, 2, 0.2, 0.2), oma = c(0, 1.8, 0, 0))
  for(j in 1 : d) {
    sxy = cov_jhdt_mfbm(n, h, H, sig, rho, eta, j, delta)
    w = sxy %*% syy_inv
    if(plot == TRUE) {
      for (i in 1 : d){
        foo = w[seq(i, length(w), by = d)]
        plot(foo, type = "l")
        abline(h = 0, col = "red", lwd = 0.1)
        if(i == 1) mtext(colnames(my)[j], side = 2, line = 2.5, cex = 1)
      }
    }
    fcast[j] =
      w %*%  c(t(my) - mu) + mu[j]
  }
  if(plot == TRUE) par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
  names(fcast) <- names(H)
  return(fcast)
}

#' @export
#' @rdname fcast_mfbm
fcast_se_mfbm <- function(n, h, H, sig, rho, eta, delta) {
  sxx = cov_0_mfbm(n, H, sig, rho, eta, delta)
  sxy = t(sapply(1 : d, cov_jhdt_mfbm, n = n, h = h, H = H, sig = sig,
               rho = rho, eta = eta, delta = delta))
  syy = cov_dt_mfbm(n, H, sig, rho, eta, delta)
  syy_inv = solve(syy)
  return(sxx - sxy %*% syy_inv %*% t(sxy))
}

#' Forecasting the mfOU using Gaussian conditioning
#'
#' @description
#' Using a standard result for Gaussian random vectors, these functions allow
#' us to compute point forecasts and associated standard errors in the framework
#' of the multivariate fractional Ornstein-Uhlenbeck process of dimension \eqn{d}.
#'
#' Given that the vector that we want to predict, \eqn{X:=Y_{(n+h)\Delta}\in\mathbb{R}^d}, and the
#' history of information, \eqn{Y:=\left(Y_{n\Delta},\dots, Y_{\Delta}\right)\in\mathbb{R}^{d\times n}},
#' are jointly Gaussian, we use that:
#'
#' \deqn{\mu_{X|Y}=\mu_X+\Sigma_{XY}\Sigma_{YY}^{-1}\left(Y-\mu_Y\right),}
#' \deqn{\Sigma_{XX|Y}=\Sigma_{XX}-\Sigma_{XY}\Sigma_{YY}^{-1}\Sigma_{XY}^T,}
#' where we used the notation \eqn{\mu_Z=\mathbb{E}(Z)} and
#' \eqn{\Sigma_{ZU} = \mathbb{E}((Z-\mu_Z)(U-\mu_U)^T)}.
#'
#' - `cov_jhdt_mfou` calculates the matrix \eqn{\Sigma_{XY}} described above;
#' - `cov_dt_mfou` calculates the matrix \eqn{\Sigma_{YY}} described above;
#' - `cov_0_mfou` calculates the matrix \eqn{\Sigma_{XX}} described above;
#' - `fcast_point_mfou` calculates point forecasts, i.e. the vector \eqn{\mu_{X|Y}} described above;
#' - `fcast_se_mfou` calculates standard errors of the forecasts, i.e. the matrix \eqn{\Sigma_{XX|Y}} described above;
#'
#' The functions `fcast_point_mfou` and `fcast_se_mfou` are not efficient in
#' rolling window applications with fixed parameters. In which case,
#' efficiency gains could often be obtained by avoiding redundant matrix operations.
#'
#' @param n Integer. Length of the time series.
#' @param h Integer. Forecast horizon (e.g. 1, 5, 10 steps ahead).
#' @param a Vector of speed of mean reversion coefficients.
#' @param nu Vector of diffusion coefficients.
#' @param H Vector of Hurst coefficients.
#' @param rho Symmetric matrix of instantaneous correlation coefficients.
#' @param eta Antisymmetric matrix of asymmetry coefficients.
#' @param mu Vector of long term means. Default is NULL, in which case it is computed as sample average.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years).
#' @param my Numeric matrix containing the time series in its columns.
#' @param j Component of the system. Integer between 1 and \eqn{d} or column name.
#' @param plot Logical. Default is FALSE. If TRUE, a panel of plots will exhibit
#' the weights of the linear combination.
#'
#' @returns Numeric matrices.
#' @export
#' @rdname fcast_mfou
#' @references
#' https://en.wikipedia.org/wiki/Multivariate_normal_distribution
#' @examples
#' # generating the data
#' n = 1000
#' delta = 1/252;
#' H = c(0.4, 0.2, 0.6)
#' rho = matrix(c(1, .2, .7, 0.2, 1, 0.6, 0.7, 0.6, 1), nc = 3)
#' eta = matrix(c(0, .15, 0, - 0.15, 0, - 0.02, 0, 0.02, 0), nc = 3)
#' a = c(1, 2, 0.05); nu = c(8, 2, 3); mu = c(0, 0, 0)
#' mfou = sim_mfou_exa(n, H, rho, eta, a, nu, mu, delta)
#'
#' # rolling window forecast using n observations at a time
#' d = ncol(mfou)
#' n = 512; h = 1
#' fcast = cbind(mfou, matrix(NA, nrow = nrow(mfou), ncol = d))
#' colnames(fcast) <- c(paste0("x_", 1 : d), paste0("f_", 1 : d))
#' syy = cov_dt_mfou(n = n, a = a, nu = nu, H = H, rho = rho, eta = eta, delta = delta)
#' sxy = t(sapply(X = 1 : d, FUN = cov_jhdt_mfou,
#'                n = n, h = h, a = a, nu = nu, H = H, rho = rho, eta = eta, delta = delta))
#' syy_inv = solve(syy)
#' w = sxy %*% syy_inv
#' for(f in (n + 1) : nrow(fcast)) {
#'   tmp = t(fcast[(f - n) : (f - 1), 1 : d])
#'   fcast[f, (d + 1) : (2 * d)] =
#'     w %*% c(tmp)
#' }
#'
#' # plot
#' fcast |>
#'   as.data.frame() |>
#'   dplyr::mutate(t = dplyr::row_number()) |>
#'   tidyr::pivot_longer(!t,
#'     values_to = "val",
#'     names_to = c("s", "d"),
#'     names_pattern = "([^_]+)_([^_]+)") |>
#'   dplyr::filter(val != 0) |>
#'   dplyr::mutate(
#'     s = factor(s, levels = c("x", "f"))  # ensures 'f' is drawn *after* 'x') |>
#'   ggplot2::ggplot(ggplot2::aes(x = t, y = val, col = s)) +
#'   ggplot2::geom_line() +
#'   ggplot2::facet_wrap(~ d, nrow = 3) +
#'   ggplot2::theme_minimal()
cov_dt_mfou <- function(n, a, nu, H, rho, eta, delta) {
  d = length(H)
  if(!all(length(H), nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  sy = matrix(NA, nr = n * d, nc = n * d)
  sk <- function(k){
    cov = matrix(NA, nc = d, nr = d)
    for (i in 1 : d){
      for (j in 1 : i) {
        cov[i, j] = cov_jik_mfou(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i,j], eta[i, j], k * delta) # opposite because of the minus sign in the notes
        if(i != j) cov[j, i] = cov_ijk_mfou(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i,j], eta[i, j], k * delta)
      }
    }
    return(cov)
  }
  sy[1 : d, ] = do.call("cbind", sapply(0 : (n - 1), sk, simplify = FALSE))
  for(k in 2 : n) {
    sy[((k - 1) * d + 1) : (k * d), ((k - 1) * d + 1) : (n * d)] <-
      sy[1 : d, 1 : (n * d - (k - 1) * d)]
  }
  sy[lower.tri(sy)] <- t(sy)[lower.tri(sy)]
  return(sy)
}

#' @export
#' @rdname fcast_mfou
cov_dt_fou <- function(n, a, nu, H, delta) {
  sy = matrix(NA, n, n)
  sy[1, ] = sapply((0 : (n - 1)) * delta, cov_iik_mfou, ai = a, nui = nu, Hi = H)
  for(k in 2 : n) {
    sy[k, k : n] <-
      sy[1, 1 : (n - (k - 1))]
  }
  sy[lower.tri(sy)] <- t(sy)[lower.tri(sy)]
  return(sy)
}

#' @export
#' @rdname fcast_mfou
cov_jhdt_mfou <- function(n, h, a, nu, H, rho, eta, delta, j) {
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  sxy = numeric()
  gk <- function(k) {
    tmp = numeric(d)
    for(i in 1 : d) {
      tmp[i] = cov_ijk_mfou(a[j], a[i], nu[j], nu[i], H[j], H[i], rho[j, i], eta[j, i], k * delta)
    }
    return(tmp)
  }
  sxy = unlist(sapply((n + h - 1) : h, gk, simplify = FALSE))
  return(sxy)
}


#' @export
#' @rdname fcast_mfou
cov_0_mfou <- function(a, nu, H, rho, eta) {
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  cov = matrix(NA, nc = d, nr = d)
  for (i in 1 : d){
    for (j in i : d) {
      cov[i, j] = cov_ij0_mfou(a[i], a[j], nu[i], nu[j], H[i], H[j],
                                rho[i, j], eta[i, j])
    }
  }
  cov[lower.tri(cov)] = t(cov)[lower.tri(cov)]
  return(cov)
}


#' @export
#' @rdname fcast_mfou
fcast_point_mfou <- function(my, h, a, nu, H, rho, eta, delta,
                             mu = NULL, plot = FALSE) {
  if(sum(is.na(my)) > 0) stop("Error: your time series contain NAs.")
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions doesn't match.")
  fcast = rep(NA, d)
  n = nrow(my)
  if(is.null(mu)) mu = apply(my, 2, mean)
  syy = cov_dt_mfou(n, a, nu, H, rho, eta, delta)
  syy_inv = solve(syy)
  if(plot == TRUE) par(mfrow = c(d, d), mar = c(2, 2, 0.2, 0.2), oma = c(0, 1.8, 0, 0))
  for(j in 1 : d) {
    sxy = cov_jhdt_mfou(n, h, a, nu, H, rho, eta, delta, j)
    w = sxy %*% syy_inv
    if(plot == TRUE) {
      for (i in 1 : d){
        foo = w[seq(i, length(w), by = d)]
        plot(foo, type = "l")
        abline(h = 0, col = "red", lwd = 0.1)
        if(i == 1) mtext(colnames(my)[j], side = 2, line = 2.5, cex = 1)
      }
    }
    fcast[j] =
      w %*%  c(t(my) - mu) + mu[j]
  }
  if(plot == TRUE) par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
  names(fcast) <- names(H)
  return(fcast)
}

#' @export
#' @rdname fcast_mfou
fcast_se_mfou <- function(n, h, a, nu, H, rho, eta, delta) {
  sxx = cov_0_mfou(a, nu, H, rho, eta)
  sxy = t(sapply(1 : d, cov_jhdt_mfou, n = n, h = h, a = a, nu = nu, H = H,
                 rho = rho, eta = eta, delta = delta))
  syy = cov_dt_mfou(n, a, nu, H, rho, eta, delta)
  syy_inv = solve(syy)
  return(sxx - sxy %*% syy_inv %*% t(sxy))

}

#' Forecasting the slow mfOU using Gaussian conditioning
#'
#' @description
#' Using a standard result for Gaussian random vectors, these functions allow
#' us to compute point forecasts and associated standard errors in the framework
#' of the slow multivariate fractional Ornstein-Uhlenbeck process of dimension \eqn{d}.
#' This process is constructed by leveraging the small mean reversion coefficient
#' asymptotic result of the cross-covariance function .
#'
#' Given that the vector that we want to predict, \eqn{X:=Y_{(n+h)\Delta}\in\mathbb{R}^d}, and the
#' history of information, \eqn{Y:=\left(Y_{n\Delta},\dots, Y_{\Delta}\right)\in\mathbb{R}^{d\times n}},
#' are jointly Gaussian, we use that:
#'
#' \deqn{\mu_{X|Y}=\mu_X+\Sigma_{XY}\Sigma_{YY}^{-1}\left(Y-\mu_Y\right),}
#' \deqn{\Sigma_{XX|Y}=\Sigma_{XX}-\Sigma_{XY}\Sigma_{YY}^{-1}\Sigma_{XY}^T,}
#' where we used the notation \eqn{\mu_Z=\mathbb{E}(Z)} and
#' \eqn{\Sigma_{ZU} = \mathbb{E}((Z-\mu_Z)(U-\mu_U)^T)}.
#'
#' - `cov_jhdt_mfou_asy` calculates the matrix \eqn{\Sigma_{XY}} described above;
#' - `cov_dt_mfou_asy` calculates the matrix \eqn{\Sigma_{YY}} described above;
#' - `cov_0_mfou_asy` calculates the matrix \eqn{\Sigma_{XX}} described above;
#' - `fcast_point_mfou_asy` calculates point forecasts, i.e. the vector \eqn{\mu_{X|Y}} described above;
#' - `fcast_se_mfou_asy` calculates standard errors of the forecasts, i.e. the matrix \eqn{\Sigma_{XX|Y}} described above;
#'
#' The functions `fcast_point_mfou_asy` and `fcast_se_mfou_asy` are not efficient in
#' rolling window applications with fixed parameters. In which case,
#' efficiency gains could often be obtained by avoiding redundant matrix operations.
#'
#' @param n Integer. Length of the time series.
#' @param h Integer. Forecast horizon (e.g. 1, 5, 10 steps ahead).
#' @param nu Vector of diffusion coefficients.
#' @param H Vector of Hurst coefficients.
#' @param rho Symmetric matrix of instantaneous correlation coefficients.
#' @param eta Antisymmetric matrix of asymmetry coefficients.
#' @param cov0 Symmetric matrix of lag 0 covariances.
#' @param mu Vector of long term means. Default is NULL, in which case it is computed as sample average.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years).
#' @param my Numeric matrix containing the time series in its columns.
#' @param j Component of the system. Integer between 1 and \eqn{d} or column name.
#' @param plot Logical. Default is FALSE. If TRUE, a panel of plots will exhibit
#' the weights of the linear combination.
#'
#' @returns Numeric matrices.
#' @export
#' @rdname fcast_mfou_asy
cov_dt_mfou_asy <- function(n, nu, H, rho, eta, cov0, delta) {
  d = length(H)
  if(!all(length(H), nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  sy = matrix(NA, nr = n * d, nc = n * d)
  sk <- function(k){
    cov = matrix(NA, nc = d, nr = d)
    for (i in 1 : d){
      for (j in 1 : i) {
        cov[i, j] = cov_asy_ijk_mfou(nu[i], nu[j], H[i], H[j], rho[i,j], eta[i, j], cov0[i, j], k * delta)
        if(i != j) cov[j, i] = cov_asy_jik_mfou(nu[i], nu[j], H[i], H[j], rho[i,j], eta[i, j], cov0[i, j], k * delta)
      }
    }
    return(cov)
  }
  sy[1 : d, ] = do.call("cbind", sapply(0 : (n - 1), sk, simplify = FALSE))
  for(k in 2 : n) {
    sy[((k - 1) * d + 1) : (k * d), ((k - 1) * d + 1) : (n * d)] <-
      sy[1 : d, 1 : (n * d - (k - 1) * d)]
  }
  sy[lower.tri(sy)] <- t(sy)[lower.tri(sy)]
  return(sy)
}

#' @export
#' @rdname fcast_mfou_asy
cov_jhdt_mfou_asy <- function(n, h, nu, H, rho, eta, cov0, delta, j) {
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions don't match.")
  sxy = numeric()
  gk <- function(k) {
    tmp = numeric(d)
    for(i in 1 : d) {
      tmp[i] = cov_asy_jik_mfou(nu[j], nu[i], H[j], H[i], rho[j, i], eta[j, i], cov0[i, j], k * delta)
      # calculations suggest ij because it's the first component moving forward. However, it works with ji :/
    }
    return(tmp)
  }
  sxy = unlist(sapply((n + h - 1) : h, gk, simplify = FALSE))
  return(sxy)
}

#' @export
#' @rdname fcast_mfou_asy
fcast_point_mfou_asy <- function(my, h, nu, H, rho, eta, cov0, delta,
                             mu = NULL, plot = FALSE) {
  if(sum(is.na(my)) > 0) stop("Error: your time series contain NAs.")
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions doesn't match.")
  fcast = rep(NA, d)
  n = nrow(my)
  if(is.null(mu)) mu = apply(my, 2, mean)
  syy = cov_dt_mfou_asy(n, nu, H, rho, eta, cov0, delta)
  syy_inv = solve(syy)
  if(plot == TRUE) par(mfrow = c(d, d), mar = c(2, 2, 0.2, 0.2), oma = c(0, 1.8, 0, 0))
  for(j in 1 : d) {
    sxy = cov_jhdt_mfou_asy(n, h, nu, H, rho, eta, cov0, delta, j)
    w = sxy %*% syy_inv
    if(plot == TRUE) {
      for (i in 1 : d){
        foo = w[seq(i, length(w), by = d)]
        plot(foo, type = "l")
        abline(h = 0, col = "red", lwd = 0.1)
        if(i == 1) mtext(colnames(my)[j], side = 2, line = 2.5, cex = 1)
      }
    }
    fcast[j] =
      w %*%  c(t(my) - mu) + mu[j]
  }
  if(plot == TRUE) par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
  names(fcast) <- names(H)
  return(fcast)
}


#' @export
#' @rdname fcast_mfou_asy
fcast_se_mfou <- function(n, h, nu, H, rho, eta, cov0, delta) {
  sxy = t(sapply(1 : d, cov_jhdt_mfou_asy, n = n, h = h, nu = nu, H = H,
                 rho = rho, eta = eta, cov0 = cov0, delta = delta))
  syy = cov_dt_mfou_asy(n, nu, H, rho, eta, cov0, delta)
  syy_inv = solve(syy)
  return(cov0 - sxy %*% syy_inv %*% t(sxy))

}
