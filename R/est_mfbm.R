#' Method of moment estimation of the mfBm
#'
#' @description
#' Method of moment estimator of the multivariate fractional Brownian motion.
#' - `est_H_fbm` estimates the Hurst coefficient of a fractional Brownian motion;
#' - `est_sig_fbm`  estimates the scale coefficient of a fractional Brownian
#'  motion (also its variance at \eqn{t=1});
#' - `est_rhoij_mfbm` estimates the (contemporaneous) correlation coefficient between
#' two fractional Brownian motions;
#' - `est_etaij_mfbm` estimates the asymmetry coefficient between
#' two fractional Brownian motions;
#' - `est_rho_eta_mfbm` combines the two functions above to estimates the matrices
#' \eqn{\rho} and \eqn{\eta} for a multivariate fractional Brownian motion of
#' arbitrary dimension.
#'
#'
#' @param x Numeric vector representing a time series.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years).
#' @param xi Numeric vector representing a time series.
#' @param xj Numeric vector of the same length of xj representing a time series.
#' @param mx Numeric matrix containing the time series in its columns.
#'
#' @returns Numeric value, except for `est_rho_eta_mfbm` that returns a list of
#' numeric matrices.
#' @export
#' @rdname est_mfbm
#'
#' @references
#' Bibinger, Markus, Jun Yu, and Chen Zhang. "Modeling and Forecasting Realized Volatility with Multivariate Fractional Brownian Motion." arXiv preprint arXiv:2504.15985 (2025).
est_H_fbm <- function(x) {
  n = length(x)
  num = sum((x[3 : n] - x[1 : (n - 2)]) ^ 2, na.rm = TRUE)
  den = sum((x[2 : n] - x[1 : (n - 1)]) ^ 2, na.rm = TRUE)
  return(log(num/den) / (2 * log(2)))
}

#' @export
#' @rdname est_mfbm
est_sig_fbm <- function(x, H = NULL, delta = 1) {
  if(is.null(H)) H = est_H_fbm(x)
  n = length(x)
  foo = sum((x[2 : n] - x[1 : (n - 1)]) ^ 2, na.rm = TRUE)
  return(sqrt(foo / (n * delta ^ (2 * H))))
}

#' @export
#' @rdname est_mfbm
est_rhoij_mfbm <- function(xi, xj) {
  if(length(xi) != length(xj)) print("xi and xj must have same length.")
  n = length(xi)
  num = sum((xi[2 : n] - xi[1 : (n - 1)]) * (xj[2 : n] - xj[1 : (n - 1)]), na.rm = TRUE)
  den1 = sum((xi[2 : n] - xi[1 : (n - 1)]) ^ 2, na.rm = TRUE)
  den2 = sum((xj[2 : n] - xj[1 : (n - 1)]) ^ 2, na.rm = TRUE)
  return(num / sqrt(den1 * den2))
}

#' @export
#' @rdname est_mfbm
est_etaij_mfbm <- function(xi, xj) {
  if(length(xi) != length(xj)) print("xi and xj must have same length.")
  n = length(xi)
  dxi = xi[2 : n] - xi[1 : (n - 1)]
  dxj = xj[2 : n] - xj[1 : (n - 1)]
  ddxi = xi[3 : n] - xi[1 : (n - 2)]
  ddxj = xj[3 : n] - xj[1 : (n - 2)]
  num = sum(dxj[2 : (n - 1)] * dxi[1 : (n - 2)] - dxi[2 : (n - 1)] * dxj[1 : (n - 2)], na.rm = TRUE)
  den = sqrt(sum(ddxi[1 : (n - 2)] ^ 2, na.rm = TRUE) * sum(ddxj[1 : (n - 2)] ^ 2, na.rm = TRUE)) -
    2 * sqrt(sum(dxi[1 : (n - 1)] ^ 2, na.rm = TRUE) * sum(dxj[1 : (n - 1)] ^ 2, na.rm = TRUE))
  return(num / den)
}

#' @export
#' @rdname est_mfbm
est_rho_eta_mfbm <- function(mx){
  if(is.null(colnames(mx))) colnames(mx) <- paste0("x", 1 : ncol(mx))
  mrho <- meta <-
    matrix(NA, nrow = ncol(mx), ncol = ncol(mx))
  colnames(mrho) <- colnames(meta) <- colnames(mx)
  rownames(mrho) <- rownames(meta) <- colnames(mx)
  for(i in 2 : ncol(mx)){
    for(j in 1 : (i - 1)){
      mrho[i, j] = est_rhoij_mfbm(mx[, i], mx[, j])
      meta[i, j] = est_etaij_mfbm(mx[, i], mx[, j])
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
