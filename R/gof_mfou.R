#' Goodness of covariance fit for fOU
#'
#' @description
#' These functions evaluate (graphically) the goodness of fit of estimated
#' parameters of a fractional Ornstein-Uhlenbeck process in terms of the distance
#' between theoretical (model) and empirical (data) autocovariances.
#' - `gof_fou` considers exact autocovariances;
#' - `gof_asy_fou` considers asymptotic autocovariances arising in the regime
#' of small speed of mean reversion.
#'
#' @param y Numeric vector representing the time series.
#' @param est Vector of estimated parameters. In the order: Hurst coefficient,
#' diffusion coefficient, and speed of mean reversion (the latter is unnecessary
#' for the asymptotic version) or simply output from `est_mm_fou`
#' @param lag_max Positive integer. All the lags between 0 (included) and `lag_max` are considered.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years), \eqn{\delta\in\mathbb{R}}.
#'
#' @returns Plot and numeric matrix containing the theoretical and empirical cross-covariances.
#' @export
#' @rdname gof_fou
#'
#' @examples
#' y = lrv[, 2] - mean(lrv[, 2])
#' est = est_mm_fou(y)
#' gof_fou(y, est) # exact version
#' gof_asy_fou(y, est) # asymptotic version
gof_fou <- function(y, est, lag_max = 50, delta = 1/252){
  par(mfrow = c(1, 1), mar = c(3, 3.55, 0.2, 0.2), mgp = c(2.1, 0.8, 0), cex.lab = 1.2)
  acf_theo = acf_fou(est[3], est[2], est[1], lag_max, delta)
  acf_emp = sapply(X = 0 : lag_max, FUN = cov_ijk, yi = y, yj = y, mui = 0, muj = 0)
  plot(x = 0 : lag_max, y = acf_emp, type = "h",
       ylim = c(0, max(acf_theo, acf_emp)),
       ylab = "gamma_{i,i}(k)",
       xlab = "k")
  lines(x = 0 : lag_max, y = acf_theo, col = "red")
  mtext(paste0("Autocovariance"),
        side = 3, line = 0.5, adj = 0, cex = 1.2, font = 1)
  legend("topright", legend = c("empirical", "theoretical"),
         col = c("black", "red"), lwd = c(1, 2), pt.cex = 2, bty = "n")
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
  return(
    cbind(
      "theo" = acf_theo,
      "emp" = acf_emp
    )
  )
}

#' Goodness of covariance fit for mfOU
#'
#' @description
#' These functions evaluate (graphically) the goodness of fit of estimated
#' parameters of a multivariate fractional Ornstein-Uhlenbeck process in
#' terms of the distance between theoretical (model) and empirical (data)
#' cross-covariances in a pairwise fashion.
#' - `gof_mfou` considers exact cross-covariances;
#' - `gof_asy_mfou` considers asymptotic cross-covariances arising in the regime
#' of small speed of mean reversion.
#'
#' @param my Numeric matrix containing the time series of the multivariate system as columns.
#' @param est List of estimated parameters, possibly from `est_gmm_mfou` called on `my.
#' @param si Index identifying the component i of the system (or its column name in my)
#' @param sj Index identifying the component j of the system (or its column name in my)
#' @param lag_max Positive integer. All the lags between 0 (included) and `lag_max` are considered.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years), \eqn{\delta\in\mathbb{R}}.
#'
#' @returns Numeric matrix containing the theoretical and empirical cross-covariances and plots.
#' @export
#' @rdname gof_mfou
#'
#' @examples
#' # exact version
#' my = t(t(lrv[, c(2, 3)]) - apply(lrv[, c(2, 3)], 2, mean))
#' est = est_gmm_mfou(my)
#' gof_mfou(my, est)
#'
#' # asymptotic version
#' my = t(t(lrv[, c(2, 3)]) - apply(lrv[, c(2, 3)], 2, mean))
#' est = est_gmm_mfou(my, type = "asy")
#' gof_asy_mfou(my, est)
gof_mfou <- function(my, est, si = 1, sj = 2, lag_max = 50, delta = 1/252){
  par(mfrow = c(1, 1), mar = c(3, 3.55, 2, 0.2), mgp = c(2.1, 0.8, 0), cex.lab = 1.2)
  ccf_theo = ccf_mfou(est$univ[3, si], est$univ[3, sj],
                      est$univ[2, si], est$univ[2, sj],
                      est$univ[1, si], est$univ[1, sj],
                      est$rho[si, sj], est$eta[si, sj],
                      lag_max, delta)
  ccf_emp = sapply(X = - lag_max : lag_max, FUN = cov_ijk,
                   yi = my[, si], yj = my[, sj], mui = 0, muj = 0)
  names(ccf_emp) <- - lag_max : lag_max
  plot(x = - lag_max : lag_max, y = ccf_emp, type = "h",
       ylim = c(0, max(ccf_emp, ccf_theo)),
       ylab = "gamma_{i,j}(k)",
       xlab = "k", xaxt = "n", yaxt="n")
  lines(x = - lag_max : lag_max, y = ccf_theo, col = "red")
  mtext(paste0("Cross-covariance between i = ", si, " and j = ", sj),
        side = 3, line = 0.5, adj = 0, cex = 1.2, font = 1)
  legend("topright", legend = c("empirical", "theoretical"),
         col = c("black", "red"), lwd = c(1, 2), pt.cex = 2, bty = "n")
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
  return(
    cbind(
      "theo" = ccf_theo,
      "emp" = ccf_emp
      )
  )
}

#' Check for slow mean reversion
#'
#' @description
#' These functions verify (graphically) the regime of slow mean reversion, which consists
#' in small speed of mean reversion coefficients with respect to the time scale considered.
#' This regime is characterized by the linearity of the autocovariances and
#' cross-covariances when considered against a suitable power of the lags,
#' determined by the Hurst coefficients.
#' - `lin_acf_fou` evaluates the condition on the autocovariance of a univariate component;
#' - `lin_ccf_mfou` evaluates the condition on the cross-covariances between two components;
#'
#' The conditions should be considered for all components and pairs of components
#' before modeling the entire system in the asymptotic regime.
#'
#' @param y Numeric vector representing a time series.
#' @param H Hurst coefficient. Default is NULL, in which case H is obtained with the COF estimator.
#' @param yi Numeric vector representing a time series.
#' @param yj Numeric vector of the same length of yi representing a time series.
#' @param Hi Hurst coefficient of component i. Default is NULL, in which case Hi is obtained with the COF estimator.
#' @param Hj Hurst coefficient of component j. Default is NULL, in which case Hj is obtained with the COF estimator.
#' @param lag_max Positive integer. All the lags between 0 (included) and `lag_max` are considered.
#'
#' @returns Numeric matrices containing theoretical and empirical covariances.
#' @export
#' @rdname lin_cov
#'
#' @examples
#' # univariate case
#' y = lrv[, 2] - mean(y = lrv[, 2])
#' lin_acf_fou(y)
#'
#' # multivariate case
#' yi = lrv[, 1]
#' yj = lrv[, 2]
#' lin_ccf_mfou(yi, yj)
lin_acf_fou <- function(y, H = NULL, lag_max = 50){
  par(mfrow = c(1, 1), mar = c(3, 3.55, 0.2, 0.2), mgp = c(2.1, 0.8, 0), cex.lab = 1.2)
  if(is.null(H)) H = est_H_fou(y)
  x = (0 : lag_max) ^ (2 * H)
  emp = sapply(X = 0 : lag_max, FUN = cov_ijk, yi = y, yj = y, mui = 0, muj = 0)
  fit = lm(emp ~ x)$fitted
  plot(x = x, y = emp, type = "h",
       ylim = c(0, max(emp)),
       ylab = "gamma_{i,i}(k)",
       xlab = "k^{2H_i}")
  lines(x = x, y = fit, col = "red")
  legend("topright", legend = c("empirical", "linear fit"),
         col = c("black", "red"), lwd = c(1, 2), pt.cex = 2, bty = "n")
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
  }


#' @export
#' @rdname lin_cov
lin_ccf_mfou <- function(yi, yj, Hi = NULL, Hj = NULL, lag_max = 50){
  if(is.null(Hi)) Hi = est_H_fou(yi)
  if(is.null(Hj)) Hj = est_H_fou(yj)
  par(mfrow = c(1, 1), mar = c(3, 3.55, 0.2, 0.2), mgp = c(2.1, 0.8, 0), cex.lab = 1.2)
  x = (0 : lag_max) ^ (Hi + Hj)
  empi = sapply(X = 0 : lag_max, FUN = cov_ijk, yi = yi, yj = yj, mui = 0, muj = 0)
  empj = sapply(X = 0 : lag_max, FUN = cov_ijk, yi = yj, yj = yi, mui = 0, muj = 0)
  fiti = lm(empi ~ x)$fitted
  fitj = lm(empj ~ x)$fitted
  plot(x = c(- rev(x), x[- 1]), y = c(rev(empj), empi[- 1]),
       type = "h",
       ylab = "gamma_{i,j}(k)$",
       xlab = "sign{(k)}|k|^{H_i+H_j}")
  lines(x = c(- rev(x), x[- 1]),
        y = c(rep(NA, lag_max), fiti), col = "red")
  lines(x = c(- rev(x), x[- 1]),
        y = c(rev(fitj), rep(NA, lag_max)), col = "red")
  legend("topright", legend = c("empirical", "linear fit"),
         col = c("black", "red"), lwd = c(1, 2), pt.cex = 2, bty = "n")
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
}

#' @export
#' @rdname gof_fou
gof_asy_fou <- function(y, est, lag_max = 50, delta = 1/252){
  par(mfrow = c(1, 1), mar = c(3, 3.55, 0.2, 0.2), mgp = c(2.1, 0.8, 0), cex.lab = 1.2)
  acf_theo = acf_asy_fou(est[2], est[1], var(y, na.rm = TRUE), lag_max, delta)
  acf_emp = sapply(X = 0 : lag_max, FUN = cov_ijk, yi = y, yj = y, mui = 0, muj = 0)
  plot(x = 0 : lag_max, y = acf_emp, type = "h",
       ylim = c(0, max(acf_theo, acf_emp)),
       ylab = "gamma_{i,i}(k)",
       xlab = "k")
  lines(x = 0 : lag_max, y = acf_theo, col = "red")
  mtext(paste0("Autocovariance"),
        side = 3, line = 0.5, adj = 0, cex = 1.2, font = 1)
  legend("topright", legend = c("empirical", "theoretical asy"),
         col = c("black", "red"), lwd = c(1, 2), pt.cex = 2, bty = "n")
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
  return(
    cbind(
      "theo" = acf_theo,
      "emp" = acf_emp
    )
  )
}

#' @export
#' @rdname gof_mfou
gof_asy_mfou <- function(my, est, si = 1, sj =2, lag_max = 50, delta = 1/252){
  par(mfrow = c(1, 1), mar = c(3, 3.55, 2, 0.2), mgp = c(2.1, 0.8, 0), cex.lab = 1.2)
  ccf_theo = ccf_asy_mfou(
    est$univ[2, si], est$univ[2, sj],
    est$univ[1, si], est$univ[1, sj],
    est$rho[si, sj], est$eta[si, sj],
    est$cov[si, sj], lag_max, delta)
  ccf_emp = sapply(X = - lag_max : lag_max, FUN = cov_ijk, yi = my[, si], yj = my[, sj],
                   mui = 0, muj = 0)
  names(ccf_emp) <- - lag_max : lag_max
  plot(x = - lag_max : lag_max, y = ccf_emp, type = "h",
       ylim = c(0, max(ccf_emp, ccf_theo)),
       ylab = "gamma_{i,j}(k)",
       xlab = "k")
  lines(x = - lag_max : lag_max, y = ccf_theo, col = "red")
  mtext(paste0("Cross-covariance between ", si, " and ", sj),
        side = 3, line = 0.5, adj = 0, cex = 1.2, font = 1)
  legend("topright", legend = c("empirical", "theoretical asy"),
         col = c("black", "red"), lwd = c(1, 2), pt.cex = 2, bty = "n")
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
  return(
    cbind(
      "theo" = ccf_theo,
      "emp" = ccf_emp
    )
  )
}
