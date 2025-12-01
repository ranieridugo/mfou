#' Estimator of H
#'
#' @description
#' This function calculates the change of frequency estimator (COF) of the Hurst
#'  coefficient of a fractional Ornstein-Uhlenbeck process. The result is restricted between
#' \eqn{10^{-3}} and \eqn{1 - 10^{-3}} in case of finite sample degeneracies.
#'
#' @param y Numeric vector representing a time series.
#' @param verbose Logical. If TRUE, warnings are printed if \eqn{H<0} or \eqn{H>1}.
#'
#' @returns Numeric value.
#' @export
#'
#' @references Wang, Xiaohu, Weilin Xiao, and Jun Yu. "Modeling and forecasting realized volatility with the fractional Ornstein????????Uhlenbeck process." Journal of Econometrics 232.2 (2023): 389-415.
est_H_fou <- function(y, verbose = FALSE){
  x = na.omit(y)
  x = x[is.finite(x)]
  n = length(x)
  num = sum((x[5 : n] - 2 * x[3 : (n - 2)] + x[1 : (n - 4)]) ^ 2)
  den = sum((x[3 : n] - 2 * x[2 : (n - 1)] + x[1 : (n - 2)]) ^ 2)
  H = 0.5 * log2(num/den)
  if(H > 0 & H < 1){
    return(H)
  } else if (H < 0) {
    if(verbose == TRUE){
      print("Warning: H < 0")
    }
    return(1e-3)
  } else if (H > 1) {
    if(verbose == TRUE){
      print("Warning: H > 1")
    }
    return(1 - 1e-3)
  }
}

#' Estimator of \eqn{\nu}
#'
#' @description
#' This function calculates the method of moment estimator of the diffusion
#' coefficient of a fractional Ornstein-Uhlenbeck process.
#' @param y Numeric vector representing a time series.
#' @param H Hurst coefficient or its estimate.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years), \eqn{\delta\in\mathbb{R}}.
#'
#' @returns Numeric value.
#' @export
#'
#' @references Wang, Xiaohu, Weilin Xiao, and Jun Yu. "Modeling and forecasting realized volatility with the fractional Ornstein????????Uhlenbeck process." Journal of Econometrics 232.2 (2023): 389-415.
est_nu_fou <- function(y, H, delta = 1/252){
  x = na.omit(y)
  x = x[is.finite(x)]
  n = length(x)
  num = sum((x[3 : n] - 2 * x[2 : (n - 1)] + x[1 : (n - 2)]) ^ 2)
  den = (n * (4 - 2 ^ (2 * H)) * delta ^ (2 * H))
  nu = sqrt(num / den)
  return(nu)
}

#' Estimator of \eqn{\mu}
#'
#' @description
#' This function calculates the method of moment estimator of the long-term mean
#' of a fractional Ornstein-Uhlenbeck process. It corresponds to a sample average.
#'
#'
#' @param y Numeric vector representing a time series.
#'
#' @returns Numeric value.
#' @export
#'
#' @references Wang, Xiaohu, Weilin Xiao, and Jun Yu. "Modeling and forecasting realized volatility with the fractional Ornstein????????Uhlenbeck process." Journal of Econometrics 232.2 (2023): 389-415.
est_mu_fou <- function(y){
  x = y[is.finite(y)]
  return(mean(x, na.rm = TRUE))
}

#' Estimator of \eqn{\alpha}
#'
#' @description
#' This function calculates the method of moment estimator of the speed of mean
#' reversion coefficient in a fractional Ornstein-Uhlenbeck process.
#'
#'
#' @param y Numeric vector representing a time series.
#' @param H Hurst coefficient or its estimate.
#' @param nu Diffusion coefficient or its estimate.
#'
#' @returns Numeric value.
#' @export
#'
#' @references Wang, Xiaohu, Weilin Xiao, and Jun Yu. "Modeling and forecasting realized volatility with the fractional Ornstein????????Uhlenbeck process." Journal of Econometrics 232.2 (2023): 389-415.
est_a_fou <- function(y, H, nu){
  x = na.omit(y)
  x = x[is.finite(x)]
  n = length(x)
  num = n * sum(x ^ 2) - sum(x) ^ 2
  den = n ^ 2 * nu ^ 2 * H * gamma(2 * H)
  alpha = (num / den) ^ (- 1 / (2 * H))
  return(alpha)
}

#' Two-step estimator of fOU
#'
#' @description
#' This function calculates the two-step estimator of the parameters of a
#' fractional Ornstein-Uhlenbeck (fOU) process. It combines the Change of
#' Frequency estimator of H with a method of moment estimator for the remaining
#' parameters \eqn{\nu}, \eqn{\alpha}, and \eqn{\mu}.
#'
#'
#' @param y Numeric vector representing a time series.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years), \eqn{\delta\in\mathbb{R}}.
#' @param H Optional parameter corresponding to the Hurst coefficient. If provided,
#' the procedure is one-step.
#'
#' @returns Numeric vector.
#' @export
#'
#' @references Wang, Xiaohu, Weilin Xiao, and Jun Yu. "Modeling and forecasting realized volatility with the fractional Ornstein????????Uhlenbeck process." Journal of Econometrics 232.2 (2023): 389-415.
est_mm_fou <-
  function(y, delta = 1/252, H = NULL){
    n = length(na.omit(y))
    y[is.finite(y)]
    if(!is.null(H)){
      H.p = H
    } else (
      H.p = est_H_fou(y = y)
    )
    nu.p = est_nu_fou(y = y, H = H.p, delta = delta)
    alpha.p = est_a_fou(y = y, H = H.p, nu = nu.p)
    mu.p = est_mu_fou(y = y)
    out = c("H" = H.p, "nu" = nu.p, "alpha" = alpha.p, "mu" = mu.p)
    return(out)
  }

#' Method of moment estimator for \eqn{\rho} and \eqn{\eta} of the mfOU process
#'
#' @description
#' These functions calculates the method of moment estimator of the parameters
#'  \eqn{\rho_{i,j}} and \eqn{\eta_{i,j}} characterizing the cross-covariance
#'  function of the multivariate fractional Ornstein-Uhlenbeck process. It is
#'  obtained by inverting the cross-covariances at lags \eqn{k} and \eqn{-k}.
#' - `est_rhoij_etaij_mfou` uses exact conditions;
#' - `est_asy_rhoij_etaij_mfou` uses conditions from the approximate cross-covariance
#'  function in the regime of small speed of mean reversion.
#'
#'
#' @param yi Numeric vector representing a time series.
#' @param yj Numeric vector of the same length of yi representing a time series.
#' @param ai Speed of mean reversion of component i.
#' @param aj Speed of mean reversion of component j.
#' @param nui Diffusion coefficient of component i.
#' @param nuj Diffusion coefficient of component j.
#' @param Hi Hurst coefficient of component i.
#' @param Hj Hurst coefficient of component j.
#' @param k lag for the covariance function that is inverted, \eqn{k>0}.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years), \eqn{\delta\in\mathbb{R}}.
#'
#' @returns Numeric vector.
#' @rdname est_ij
#' @export
#'
#' @references Dugo, Ranieri, Giacomo Giorgio, and Paolo Pigato. "The multivariate fractional Ornstein-Uhlenbeck process." arXiv preprint arXiv:2408.03051 (2024).
# est_rhoij_etaij_mfou <-
#   function(yi, yj, ai, aj, nui, nuj, Hi, Hj,
#            k = 1, delta = 1/252) {
#     if (length(yi) != length(yj)) stop("yi and yj must have same length.")
#     H = Hi + Hj
#     my = na.omit(cbind(yi, yj))
#     my = my[is.finite(rowSums(my)), ]
#     n = nrow(my)
#     Iijk = I_ij(ai, aj, Hi, Hj, delta * k)
#     Ijik = I_ji(ai, aj, Hi, Hj, delta * k)
#     g120 = cov(my[, 1], my[, 2])
#     g12k = cov(my[(k + 1) : n, 1], my[1 : (n - k), 2])
#     g21k = cov(my[1 : (n - k), 1], my[(k + 1) : n, 2])
#     foo = 1 / (H * (H - 1) * Iijk * Ijik * nui * nuj)
#     rhoij = - foo *
#       (g120 * Iijk - exp(aj * delta * k) * g21k * Iijk + g120 * Ijik -
#          exp(ai * k * delta) * g12k * Ijik)
#     etaij = foo *
#       (g120 * Iijk - exp(aj * k * delta) * g21k * Iijk -
#          g120 * Ijik + exp(ai * k * delta) * g12k * Ijik)
#     out = c("rhoij" = unname(rhoij), "etaij" = unname(etaij))
#     return(out)
#   }
est_rhoij_etaij_mfou <-
  function(yi, yj, ai, aj, nui, nuj, Hi, Hj,
           k = 1, delta = 1/252) {
    if (length(yi) != length(yj)) stop("yi and yj must have same length.")
    H = Hi + Hj
    my = na.omit(cbind(yi, yj))
    my = my[is.finite(rowSums(my)), ]
    n = nrow(my)
    Iijk = I_ij(ai, aj, Hi, Hj, delta * k)
    Ijik = I_ji(ai, aj, Hi, Hj, delta * k)

    g120 = cov(my[, 1], my[, 2])
    g12k = cov(my[(k + 1) : n, 1], my[1 : (n - k), 2])
    g21k = cov(my[1 : (n - k), 1], my[(k + 1) : n, 2])

    a1 = - (Iijk + Ijik) / (nui * nuj * H * (H - 1) * Iijk * Ijik)
    a2 = 1 / (nui * nuj * H * (H - 1) * exp(- ai * k * delta) * Iijk)
    a3 = 1 / (nui * nuj * H * (H - 1) * exp(- aj * k * delta) * Ijik)
    b1 = (Iijk - Ijik) / (nui * nuj * H * (H - 1) * Iijk * Ijik)

    rhoij = a1 * g120 + a2 * g12k + a3 * g21k
    etaij = b1 * g120 + a2 * g12k - a3 * g21k

    out = c("rhoij" = unname(rhoij), "etaij" = unname(etaij))
    return(out)
  }

#' @rdname est_ij
#' @export
est_asy_rhoij_etaij_mfou <-
  function(yi, yj, nui, nuj, Hi, Hj,
           k = 1, delta = 1/252) {
    if (length(yi) != length(yj)) stop("yi and yj must have same length.")
    H = Hi + Hj
    my = na.omit(cbind(yi, yj))
    my = my[rowSums(my) != Inf, ]
    my = my[rowSums(my) != -Inf, ]
    n = nrow(my)
    rhoij = (2 * cov(my[, 1], my[, 2]) -
               cov(my[(1 + k) : n, 1], my[1 : (n - k), 2]) -
               cov(my[(1 + k) : n, 2], my[1 : (n - k), 1])) /
      (nui * nuj * (k * delta) ^ H)
    etaij = (cov(my[(1 + k) : n, 2], my[1 : (n - k), 1]) -
               cov(my[(1 + k) : n, 1], my[1 : (n - k), 2]) ) /
      (nui * nuj * (k * delta) ^ H)
    out = c("rhoij" = unname(rhoij), "etaij" = unname(etaij))
    return(out)
  }

#' Method of moment estimator for the matrices \eqn{\rho} and \eqn{\eta} of the mfOU process
#'
#' @description
#' These functions calculate the method of moment estimator of the parameter
#' matrices \eqn{\rho} and \eqn{\eta} characterizing the cross-covariance function
#' of the multivariate fractional Ornstein-Uhlenbeck process. These are
#' obtained by inverting the cross-covariances at lag \eqn{k} and \eqn{-k}.
#' - `est_rho_eta_mfou` uses exact conditions;
#' - `est_asy_rho_eta_mfou` uses conditions from the approximate cross-covariance
#'  function in the regime of small speed of mean reversion.
#'
#' @param my Numeric matrix containing the time series in its columns.
#' @param a Numeric vector containing the speed of mean reversion coefficients.
#' @param nu Numeric vector containing the diffusion coefficients.
#' @param H Numeric vector containing the Hurst coefficients.
#' @param k lag for the covariance function that is inverted, \eqn{k>0}.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years), \eqn{\delta\in\mathbb{R}}.
#'
#' @returns List of two numeric matrices.
#' @rdname est_cov
#' @export
#'
#' @references Dugo, Ranieri, Giacomo Giorgio, and Paolo Pigato. "The multivariate fractional Ornstein-Uhlenbeck process." arXiv preprint arXiv:2408.03051 (2024).
est_rho_eta_mfou <- function(my, a, nu, H, k = 1, delta = 1/252){
  if(is.null(colnames(my))) colnames(my) <- paste0("x", 1 : ncol(my))
  mrho <- meta <-
    matrix(NA, nrow = ncol(my), ncol = ncol(my))
  colnames(mrho) <- colnames(meta) <- colnames(my)
  rownames(mrho) <- rownames(meta) <- colnames(my)
  for(i in 2 : ncol(my)){
    for(j in 1 : (i - 1)){
      tmp = est_rhoij_etaij_mfou(my[, i], my[, j], a[i], a[j], nu[i], nu[j], H[i], H[j],
                                 k = k, delta = delta)
      mrho[i, j] = tmp[1]
      meta[i, j] = tmp[2]
    }
  }
  diag(mrho) <- 1
  diag(meta) <- 0
  mrho[upper.tri(mrho, diag = FALSE)] = t(mrho)[upper.tri(mrho, diag = FALSE)]
  meta[upper.tri(meta, diag = FALSE)] = - t(meta)[upper.tri(meta, diag = FALSE)]
  return(list(
    "rho" = mrho,
    "eta" = meta))
}

#' @rdname est_cov
#' @export
est_asy_rho_eta_mfou <- function(my, nu, H, k = 1, delta = 1/252){
  mrho <- meta <-
    matrix(NA, nrow = ncol(my), ncol = ncol(my))
  colnames(mrho) <- colnames(meta) <- colnames(my)
  rownames(mrho) <- rownames(meta) <- colnames(my)
  for(i in 2 : ncol(my)){
    for(j in 1 : (i - 1)){
      tmp = est_asy_rhoij_etaij_mfou(my[, i], my[, j], nu[i], nu[j], H[i], H[j],
                                     k = k, delta = delta)
      mrho[i, j] = tmp[1]
      meta[i, j] = tmp[2]
    }
  }
  diag(mrho) <- 1
  diag(meta) <- 0
  mrho[upper.tri(mrho, diag = FALSE)] = t(mrho)[upper.tri(mrho, diag = FALSE)]
  meta[upper.tri(meta, diag = FALSE)] = - t(meta)[upper.tri(meta, diag = FALSE)]
  return(list(
    "rho" = mrho,
    "eta" = meta))
}

#' Two-step estimator of the mfOU process
#'
#' @description
#' This function calculates the method of moment estimator of all the parameters
#' of a multivariate fractional Ornstein-Uhlenbeck process of dimension \eqn{d}.
#' The parameters govering the univariate marginal distributions,
#' \eqn{H_i,\ \nu_i,\ \alpha_i,\ \text{and}\ \mu_i,\ i=1,\dots,d}, are obtained in a
#' first step, whereas those determining the joint dynamics, \eqn{\rho\ \text{and}\ \eta},
#' are obtained pairwise using the estimates from the previous step as input.
#'
#' @param my Numeric matrix containing the time series of the multivariate system as columns.
#' @param k lag for the covariance function that is inverted, \eqn{k>0}.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years), \eqn{\delta\in\mathbb{R}}.
#' @param type One between `"exa"` and `"asy"` to determine whether the exact cross-covariance condition
#' or the asymptotic one is to be used.
#'
#' @returns List of numeric matrices
#' @export
#' @references
#' Dugo, Ranieri, Giacomo Giorgio, and Paolo Pigato. "The multivariate fractional Ornstein-Uhlenbeck process." arXiv preprint arXiv:2408.03051 (2024).
#'
#' Wang, Xiaohu, Weilin Xiao, and Jun Yu. "Modeling and forecasting realized volatility with the fractional Ornstein????????Uhlenbeck process." Journal of Econometrics 232.2 (2023): 389-415.
est_mm_mfou <-  function(my, delta = 1/252, k = 1, type = "exa"){
  if(is.null(colnames(my))) colnames(my) <- paste0("y", 1 : ncol(my))
  # Univariate estimates
  mUniv = apply(my, 2, est_mm_fou, delta = delta)
  # Pairwise estimates of rhoij and etaij
  if (type == "exa") {
    append(est_rho_eta_mfou(my, mUniv[3, ], mUniv[2, ], mUniv[1, ],
                            k = k, delta = delta),
           list("univ" = mUniv), after = 0)
  } else if (type == "asy") {
    append(
      append(
        est_asy_rho_eta_mfou(my, mUniv[2, ], mUniv[1, ],
                             k = k, delta = delta),
        list("univ" = mUniv), after = 0),
      list("cov" = cov(my, use = "complete.obs")))
  }
}

#' Generalized Method of Moments estimation of the mfOU process
#'
#' @description
#' This function estimates the parameters of a multivare fractional Ornstein-Uhlenbeck process
#' of dimension \eqn{d}, say \eqn{\theta}, by leveraging on the Generalized Method of Moments (GMM) methodology.
#' The GMM consists in minimizing a square distance between sample estimates of cross-covariances,
#' \eqn{\hat\gamma}, and theoretical ones, \eqn{\gamma(\theta)}. The vectors \eqn{\hat\gamma}
#' and \eqn{\gamma(\theta)} include cross-covariances evaluated for all combinations
#' \eqn{i,j=1,\dots,d} and all lags \eqn{k\in\mathcal{L}}.
#' The loss function being minimized is
#' \deqn{\left(\hat\gamma-\gamma(\theta)\right)^T W_n \left(\hat\gamma-\gamma(\theta)\right),}
#' for a positive semi-definite matrix \eqn{W_n}.
#' The methodology considers lags \eqn{\mathcal{L}=(0,1,2,3,5,20,50)} as in the paper.
#'
#' @param my Numeric matrix containing the time series of the multivariate system in its columns.
#' @param type One among `exa`, `asy`, `cau` to determine which version of the model the user wants to estimate
#' among exact, asymptotic, and causal.
#' @param step2 Logical. If TRUE, 2-step GMM is implemented. If FALSE, only the first step with the weighting
#' matrix given by the identity matrix weights is performed.
#' @param verbose Logical. If TRUE, optimization details are printed for each iteration.
#'
#' @returns List of numeric objects containing estimated parameters.
#' @export
#' @references
#' Dugo, Ranieri, Giacomo Giorgio, and Paolo Pigato. "Multivariate Rough Volatility." arXiv preprint arXiv:2412.14353 (2024).
#' @examples
#' # Simulated data
#' mfou <-
#'   sim_mfou(n = 10000, delta = 1/252, m = 1,
#'            H = c(0.1, 0.2, 0.3),
#'            a = c(1, 0.5, 1),
#'            nu = c(1, 0.5, 1),
#'            mu = c(0, 0, 0),
#'            rho = matrix(c(1, 0.6, 0.7, 0.6, 1, 0.8, 0.7, 0.8, 1),
#'                         ncol = 3),
#'            eta = matrix(c(0, -0.1, -0.2, 0.1, 0, -0.3, 0.2, 0.3, 0),
#'                      ncol = 3))[2000 : 10000, ]
#' est_gmm_mfou(mfou)
#'
#' # Log-realized volatilities
#' my = t(t(lrv[, c(2, 3)]) - apply(lrv[, c(2, 3)], 2, mean))
#' est_gmm_mfou(my)
est_gmm_mfou <- function(my,
                         lag_max = 5, lag_add = c(20, 50), delta = 1/252,
                         type = "exa", step2 = TRUE,
                         verbose = TRUE) {
  d = ncol(my); dd = d * (d - 1) / 2
  if(type == "exa") {
    vlb = c(rep(0.001, d), rep(0.01, d), rep(0.01, d), rep(- 0.999, dd), rep(- 2, dd))
    vub = c(rep(5, d), rep(10, d), rep(0.95, d), rep(0.999, dd), rep(2, dd))
    tet0 = est_mm_mfou(my, delta, type = "exa")
  } else if(type == "asy"){
    vlb = c(rep(0.01, d), rep(0.01, d), rep(0.01, d), rep(- 3, dd), rep(- 0.99, dd), rep(- 2, dd))
    vub = c(rep(3, d), rep(10, d), rep(0.95, d), rep(3, dd), rep(0.99, dd), rep(2, dd))
    tet0 = est_mm_mfou(my, delta, type = "asy")
  } else if(type == "cau"){
    vlb = c(rep(0.001, d), rep(0.01, d), rep(0.01, d), rep(- 0.999, dd))
    vub = c(rep(5, d), rep(10, d), rep(0.95, d), rep(0.999, dd))
    tet0 = est_mm_mfou(my, delta, type = "exa")
  }
  vtet0 = par_list2vec(tet0, d, type)
  itrace = ifelse(verbose == TRUE, 6, 0)
  # step 1
  nl = length(c(0 : 5, c(20, 50)))
  nw = nl * d + nl * dd + (nl - 1) * dd
  res1 =
    optim(
      par = vtet0,
      fn = (function(par){
        loss = errors_gmm(par, my, type, ts = FALSE,
                          lag_max, lag_add, delta)
        W = diag(nw)
        return(loss %*% W %*% loss)
      }),
      gr = (function(par){
        loss = errors_gmm(par, my, type, ts = FALSE,
                          lag_max, lag_add, delta)
        W = diag(nw)
        jac = jacobian_gmm(par, my, type,
                           lag_max, lag_add, delta)
        return(2 * t(jac) %*% W %*% loss)
      }),
      method = "L-BFGS-B",
      lower = vlb, upper = vub,
      control = list(maxit = 3000, trace = itrace)
    )
  if(step2 == TRUE){
    # weighting matrix
    S = cov_mom_nw(res1$par, my, type, verbose = verbose)
    W = solve(diag(diag(S)))
    # step2
    res2 =
      optim(
        par = vtet0,
        fn = (function(par){
          loss = errors_gmm(par, my, type, ts = FALSE,
                            lag_max, lag_add, delta)
          return(loss %*% W %*% loss)
        }),
        gr = (function(par){
          loss = errors_gmm(par, my, type, ts = FALSE,
                            lag_max, lag_add, delta)
          jac = jacobian_gmm(par, my, type,
                             lag_max, lag_add, delta)
          return(2 * t(jac) %*% W %*% loss)
        }),
        method = "L-BFGS-B",
        lower = vlb,
        upper = vub,
        control = list(maxit = 3000, trace = itrace)
      )
    out = par_vec2list(res2$par, d, type)
  } else if(step2 == FALSE) {
    out = par_vec2list(res1$par, d, type)
  }
  return(out)
}
