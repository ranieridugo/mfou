#' @export
errors_exa_gmm <- function(theta, my, ts = FALSE,
                           lag_max = 5, lag_add = c(20, 50), delta = 1/252) {
  theta = unname(theta)
  n = nrow(my)
  d = ncol(my)
  dd = d * (d - 1) / 2
  a <- theta[1 : d]
  nu <- theta[(d + 1) : (2 * d)]
  H <- theta[(2 * d + 1) : (3 * d)]
  rho <- eta <-
    matrix(nr = d, nc = d)
  rho[upper.tri(rho)] <- theta[(3 * d + 1) : (3 * d + dd)]
  eta[upper.tri(eta)] <- theta[(3 * d + dd + 1) : (3 * d + 2 * dd)]
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  eta[lower.tri(eta)] <- - t(eta)[lower.tri(eta)]
  diag(rho) <- 1
  diag(eta) <- 0

  cov_ijk_mom =
    function(i, j, l) {
      mom =
        my[(1 + l) : n, i] * my[1 : (n - l), j] -
        cov_ijk_mfou(a[i], a[j], nu[i], nu[j], H[i], H[j],
                     rho[i, j], eta[i, j], l * delta)
      fill = rep(NA, l)
      return(c(fill, mom))
    }

  vlagsij = c(0 : lag_max, lag_add)
  vlagsji = c(1 : lag_max, lag_add)
  lmoms = list()
  for(i in 1 : d){
    for(j in 1 : d){
      if(i <= j) {
        vlags = vlagsij
      } else if(i > j){
        vlags = vlagsji
      }
      lmoms[[(i - 1) * d + j]] <-
        do.call("cbind",
                lapply(X = vlags, FUN = cov_ijk_mom, i = i, j = j)
        )
      colnames(lmoms[[(i - 1) * d + j]]) <- paste(i, j, vlags, sep = "_")
    }
  }

  mmoms = do.call("cbind", lmoms)
  rm(lmoms)
  if(ts == TRUE) {
    return(mmoms[(max(0 : lag_max, lag_add) + 1)  : n, ]) #
  } else if (ts == FALSE) {
    return(apply(mmoms, 2, mean, na.rm = TRUE)) #
  }
}
#' @export
jacobian_exa_gmm <- function(theta, my,
                             lag_max = 5, lag_add = c(20, 50), delta = 1/252){
  theta = unname(theta)
  n = nrow(my)
  d = ncol(my)
  dd = d * (d - 1) / 2
  a <- theta[1 : d]
  nu <- theta[(d + 1) : (2 * d)]
  H <- theta[(2 * d + 1) : (3 * d)]
  rho <- eta <-
    matrix(nr = d, nc = d)
  rho[upper.tri(rho)] <- theta[(3 * d + 1) : (3 * d + dd)]
  eta[upper.tri(eta)] <- theta[(3 * d + dd + 1) : (3 * d + 2 * dd)]
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  eta[lower.tri(eta)] <- - t(eta)[lower.tri(eta)]
  diag(rho) <- 1
  diag(eta) <- 0

  vlags = c(0 : lag_max, lag_add)
  p = length(vlags)
  vpop = function(pos, arg){
    vec = rep(0, length(theta))
    vec[pos] = arg
    return(vec)
  }
  jac <- matrix(nr = 0, nc = length(theta))
  for(i in 1 : d){
    for(j in 1 : d){
      for(l in 1 : p){
        lag = vlags[l] * delta
        if(i == j & lag == 0){ # variance
          tmp = vpop(
            c(i, d + i, 2 * d + i),
            c(dVda(a[i], nu[i], H[i]),
              dVdv(a[i], nu[i], H[i]),
              dVdH(a[i], nu[i], H[i])))
        }else if(i == j & lag > 0){ # auto-covariance
          tmp = vpop(
            c(i, d + i, 2 * d + i),
            c(dGisda(a[i], nu[i], H[i], lag),
              dGisdv(a[i], nu[i], H[i], lag),
              dGisdH(a[i], nu[i], H[i], lag)))
        }else if(i != j & lag == 0){ # covariance (only defined for i < j)
          if(i < j) {
            k = (j - 1) * (j - 2) / 2 + i
            tmp = vpop(
              c(i, j, d + i, d + j, 2 * d + i, 2 * d + j, 3 * d + k, 3 * d + dd + k),
              c(dG0da1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0da2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dv1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dv2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dH(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dH(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0drho(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0deta(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j])))
          }else if(i > j){
            tmp = NULL
          }
        }else if(i != j & lag > 0){ # cross-covariance
          if(i < j) {
            sdeta = + 1
            k = (j - 1) * (j - 2) / 2 + i
          } else if (i > j){
            sdeta = - 1
            k = (i - 1) * (i - 2) / 2 + j
          }
          tmp = vpop(
            c(i, j, d + i, d + j, 2 * d + i, 2 * d + j, 3 * d + k, 3 * d + dd + k),
            c(
              dGs12da1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12da2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12dv1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12dv2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12dH(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12dH(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12drho(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              sdeta * dGs12deta(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag))
          )
        }
        jac = rbind(jac, name = tmp, deparse.level = 0)
        if(!is.null(tmp)) rownames(jac)[nrow(jac)] <- paste0(i, "_", j, "_", vlags[l])
      }
    }
  }
  return(- jac)
}
#' @export
errors_cau_gmm <- function(theta, my, ts = FALSE,
                           lag_max = 5, lag_add = c(20, 50), delta = 1/252) {
  theta = unname(theta)
  n = nrow(my)
  d = ncol(my)
  dd = d * (d - 1) / 2
  a <- theta[1 : d]
  nu <- theta[(d + 1) : (2 * d)]
  H <- theta[(2 * d + 1) : (3 * d)]
  rho <- eta <-
    matrix(nr = d, nc = d)
  rho[upper.tri(rho)] <- theta[(3 * d + 1) : (3 * d + dd)]
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  diag(rho) <- 1
  for(i in 1 : d){
    for(j in 1 : d){
      eta[i, j] = etaij_caus(rho[i, j], H[i], H[j])
    }
  }
  cov_ijk_mom =
    function(i, j, l) {
      mom =
        my[(1 + l) : n, i] * my[1 : (n - l), j] -
        cov_ijk_mfou(a[i], a[j], nu[i], nu[j], H[i], H[j],
                     rho[i, j], eta[i, j], l * delta)
      fill = rep(NA, l)
      return(c(fill, mom))
    }

  vlagsij = c(0 : lag_max, lag_add)
  vlagsji = c(1 : lag_max, lag_add)
  lmoms = list()
  for(i in 1 : d){
    for(j in 1 : d){
      if(i <= j) {
        vlags = vlagsij
      } else if(i > j){
        vlags = vlagsji
      }
      lmoms[[(i - 1) * d + j]] <-
        do.call("cbind",
                lapply(X = vlags, FUN = cov_ijk_mom, i = i, j = j)
        )
      colnames(lmoms[[(i - 1) * d + j]]) <- paste(i, j, vlags, sep = "_")
    }
  }

  mmoms = do.call("cbind", lmoms)
  rm(lmoms)
  if(ts == TRUE) {
    return(mmoms[(max(0 : lag_max, lag_add) + 1)  : n, ]) #
  } else if (ts == FALSE) {
    return(apply(mmoms, 2, mean, na.rm = TRUE)) #
  }
}
#' @export
jacobian_cau_gmm <- function(theta, my,
                             lag_max = 5, lag_add = c(20, 50), delta = 1/252){
  theta = unname(theta)
  n = nrow(my)
  d = ncol(my)
  dd = d * (d - 1) / 2
  a <- theta[1 : d]
  nu <- theta[(d + 1) : (2 * d)]
  H <- theta[(2 * d + 1) : (3 * d)]
  rho <- eta <- matrix(nr = d, nc = d)
  rho[upper.tri(rho)] <- theta[(3 * d + 1) : (3 * d + dd)]
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  diag(rho) <- 1
  for(i in 1 : d){
    for(j in 1 : d){
      eta[i, j] = etaij_caus(rho[i, j], H[i], H[j])
    }
  }

  vlags = c(0 : lag_max, lag_add)
  p = length(vlags)
  vpop = function(pos, arg){
    vec = rep(0, length(theta))
    vec[pos] = arg
    return(vec)
  }
  jac <- matrix(nr = 0, nc = length(theta))
  for(i in 1 : d){
    for(j in 1 : d){
      for(l in 1 : p){
        lag = vlags[l] * delta
        if(i == j & lag == 0){ # variance
          tmp = vpop(
            c(i, d + i, 2 * d + i),
            c(dVda(a[i], nu[i], H[i]),
              dVdv(a[i], nu[i], H[i]),
              dVdH(a[i], nu[i], H[i])))
        }else if(i == j & lag > 0){ # auto-covariance
          tmp = vpop(
            c(i, d + i, 2 * d + i),
            c(dGisda(a[i], nu[i], H[i], lag),
              dGisdv(a[i], nu[i], H[i], lag),
              dGisdH(a[i], nu[i], H[i], lag)))
        }else if(i != j & lag == 0){ # covariance (only defined for i < j)
          if(i < j) {
            k = (j - 1) * (j - 2) / 2 + i
            tmp = vpop(
              c(i, j, d + i, d + j, 2 * d + i, 2 * d + j, 3 * d + k),
              c(dG0da1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0da2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dv1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dv2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dH(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0dH(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j]),
                dG0cdrho(a[i], a[j], nu[i], nu[j], H[i], H[j])))
          }else if(i > j){
            tmp = NULL
          }
        }else if(i != j & lag > 0){ # cross-covariance
          if(i < j) {
            sdeta = + 1
            k = (j - 1) * (j - 2) / 2 + i
          } else if (i > j){
            sdeta = - 1
            k = (i - 1) * (i - 2) / 2 + j
          }
          tmp = vpop(
            c(i, j, d + i, d + j, 2 * d + i, 2 * d + j, 3 * d + k),
            c(
              dGs12da1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12da2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12dv1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12dv2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], eta[i, j], lag),
              dGs12cdH1(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], lag),
              dGs12cdH2(a[i], a[j], nu[i], nu[j], H[i], H[j], rho[i, j], lag),
              dGs12cdrho(a[i], a[j], nu[i], nu[j], H[i], H[j], lag))
          )
        }
        jac = rbind(jac, name = tmp, deparse.level = 0)
        if(!is.null(tmp)) rownames(jac)[nrow(jac)] <- paste0(i, "_", j, "_", vlags[l])
      }
    }
  }
  return(- jac)
}
#' @export
errors_asy_gmm <- function(theta, my, ts = FALSE,
                           lag_max = 5, lag_add = c(20, 50), delta = 1/252) {
  theta = unname(theta)
  n = nrow(my)
  d = ncol(my)
  dd = d * (d - 1) / 2
  nu <- theta[(d + 1) : (2 * d)]
  H <- theta[(2 * d + 1) : (3 * d)]
  rho <- eta <- cov <-
    matrix(nr = d, nc = d)
  cov[upper.tri(cov)] <- theta[(3 * d + 1) : (3 * d + dd)]
  diag(cov) <- theta[1 : d]
  rho[upper.tri(rho)] <- theta[(3 * d + dd + 1) : (3 * d + 2 * dd)]
  eta[upper.tri(eta)] <- theta[(3 * d + 2 * dd + 1) : (3 * d + 3 * dd)]
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  eta[lower.tri(eta)] <- - t(eta)[lower.tri(eta)]
  diag(rho) <- 1
  diag(eta) <- 0

  cov_ijk_mom =
    function(i, j, l) {
      mom =
        my[(1 + l) : n, i] * my[1 : (n - l), j] -
        cov_asy_ijk_mfou(nu[i], nu[j], H[i], H[j],
                         rho[i, j], eta[i, j], cov[i, j], l * delta)
      fill = rep(NA, l)
      return(c(fill, mom))
    }

  vlagsij = c(0 : lag_max, lag_add)
  vlagsji = c(1 : lag_max, lag_add)
  lmoms = list()
  for(i in 1 : d){
    for(j in 1 : d){
      if(i <= j) {
        vlags = vlagsij
      } else if(i > j){
        vlags = vlagsji
      }
      lmoms[[(i - 1) * d + j]] <-
        do.call("cbind",
                lapply(X = vlags, FUN = cov_ijk_mom, i = i, j = j)
        )
      colnames(lmoms[[(i - 1) * d + j]]) <- paste(i, j, vlags, sep = "_")
    }
  }
  mmoms = do.call("cbind", lmoms)
  rm(lmoms)
  if(ts == TRUE) {
    return(mmoms[(max(vlags) + 1)  : n, ]) #
  } else if (ts == FALSE) {
    return(apply(mmoms, 2, mean, na.rm = TRUE)) #
  }
}
#' @export
jacobian_asy_gmm <- function(theta, my,
                             lag_max = 5, lag_add = c(20, 50), delta = 1/252){
  theta = unname(theta)
  n = nrow(my)
  d = ncol(my)
  dd = d * (d - 1) / 2
  nu <- theta[(d + 1) : (2 * d)]
  H <- theta[(2 * d + 1) : (3 * d)]
  rho <- eta <- cov <-
    matrix(nr = d, nc = d)
  cov[upper.tri(cov)] <- theta[(3 * d + 1) : (3 * d + dd)]
  diag(cov) <- theta[1 : d]
  rho[upper.tri(rho)] <- theta[(3 * d + dd + 1) : (3 * d + 2 * dd)]
  eta[upper.tri(eta)] <- theta[(3 * d + 2 * dd + 1) : (3 * d + 3 * dd)]
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  rho[lower.tri(rho)] <- t(rho)[lower.tri(rho)]
  eta[lower.tri(eta)] <- - t(eta)[lower.tri(eta)]
  diag(rho) <- 1
  diag(eta) <- 0

  vlags = c(0 : lag_max, lag_add)
  p = length(vlags)
  vpop = function(pos, arg){
    vec = rep(0, length(theta))
    vec[pos] = arg
    return(vec)
  }
  jac <- matrix(nr = 0, nc = length(theta))
  for(i in 1 : d){
    for(j in 1 : d){
      for(l in 1 : p){
        lag = vlags[l] * delta
        if(i == j & lag == 0){ # variance
          tmp = vpop(i, 1)
        }else if(i == j & lag > 0){ # auto-covariance
          tmp = vpop(
            c(i, d + i, 2 * d + i),
            c(1,
              - nu[i] * lag ^ (2 * H[i]),
              - nu[i] ^ 2 * lag ^ (2 * H[i]) * log(lag)))
        }else if(i != j & lag == 0){ # covariance (only defined for i < j)
          if(i < j) {
            k = (j - 1) * (j - 2) / 2 + i
            tmp = vpop(3 * d + k, 1)
          }else if(i > j){
            tmp = NULL
          }
        }else if(i != j & lag > 0){ # cross-covariance
          if(i < j) {
            sdeta = + 1
            k = (j - 1) * (j - 2) / 2 + i
          } else if (i > j){
            sdeta = - 1
            k = (i - 1) * (i - 2) / 2 + j
          }
          tmp = vpop(
            c(d + i, d + j, 2 * d + i, 2 * d + j, 3 * d + k, 3 * d + dd + k, 3 * d + 2 * dd + k),
            c(
              - (rho[i, j] + eta[i, j]) / 2 * nu[j] * lag ^ (H[i] + H[j]),
              - (rho[i, j] + eta[i, j]) / 2 * nu[i] * lag ^ (H[i] + H[j]),
              - (rho[i, j] + eta[i, j]) / 2 * nu[i] * nu[j] * lag ^ (H[i] + H[j]) * log(lag),
              - (rho[i, j] + eta[i, j]) / 2 * nu[i] * nu[j] * lag ^ (H[i] + H[j]) * log(lag),
              1,
              - 0.5 * nu[i] * nu[j] * lag ^ (H[i] + H[j]),
              - sdeta * 0.5 * nu[i] * nu[j] * lag ^ (H[i] + H[j])))
        }
        jac = rbind(jac, name = tmp, deparse.level = 0)
        if(!is.null(tmp)) rownames(jac)[nrow(jac)] <- paste0(i, "_", j, "_", vlags[l])
      }
    }
  }
  return(- jac)
}

#' Newey-West estimator of the asymptotic covariance matrix of moment conditions
#'
#' @description
#' This function estimates the asymptotic covariance matrix of the moment conditions,
#' following the methodology outlined in Newey-West (1987, 1994). The moment conditions
#' are calculated for lags \eqn{\mathcal{L}=(0,1,2,3,5,20,50)} as is done in the paper.
#'
#' @param theta Vector of coefficients.
#' @param my Numeric matrix containing the time series of the multivariate system as columns.
#' @param bn Bandwidth parameter. If NULL, it is automatically selected.
#' @param verbose Logical. If TRUE, the automatically selected bandwidth is printed.
#' @param type One among `exa`, `asy`, `cau` to determine which version of the model the user wants to estimate
#' among exact, asymptotic, and causal.
#'
#' @returns Numeric matrix.
#'
#' @export
#'
#' @references
#' Newey, Whitney K., and Kenneth D. West. “A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.” Econometrica, vol. 55, no. 3, 1987, pp. 703–08. JSTOR, https://doi.org/10.2307/1913610. Accessed 22 July 2025.
#' Newey, Whitney K., and Kenneth D. West. "Automatic lag selection in covariance matrix estimation." The review of economic studies 61.4 (1994): 631-653.
cov_mom_nw <- function(theta, my, type = "exa", bn = NULL, verbose = TRUE){
  merr = errors_gmm(theta, my, type, ts = TRUE)
  ne = nrow(merr)
  if(is.null(bn)){
    lar = apply(merr, 2, ar, order.max = 1)
    nnw = round(4 * (ne / 100) ^ (2/9))
    vw = rep(0, length(lar[[1]]$resid))
    for(i in 1 : length(lar)) vw = vw + lar[[i]]$resid
    vsj = rep(NA, nnw + 1)
    for(j in 0 : nnw){
      vsj[j + 1] = 1 / (ne - 1) * sum(vw[(j + 2) : ne] * vw[2 : (ne - j)])
    }
    s1 = 2 * sum((1 : nnw) * vsj[- 1])
    s0 = vsj[1] + 2 * sum(vsj[- 1])
    if(s1 > 0) {
      bn = round(1.1447 * (s1 / s0) ^ (2/3) * as.integer(ne ^ (1 / 3)))
    } else if(s1 < 0){
      bn = 1.5 * as.integer(ne ^ (1 / 3))
    }
    if(verbose == TRUE) print(bn)
  }
  vw = 1 - ((1 : bn) / (bn + 1))
  for(i in 0 : bn) {
    tmp = 1 / (ne - 1 - i) * t(merr[(1 + i) : ne, ]) %*% merr[1 : (ne - i), ]
    if(i == 0){
      mS = tmp
    } else if(i > 0){
      mS = mS + vw[i] * (tmp + t(tmp))
    }
  }
  return(mS)
}


#' Moment conditions
#'
#' @description
#' This function computes the moment conditions of the GMM estimator, i.e. the
#' differences between empirical and theoretical cross-covariances, \eqn{\hat\gamma-\gamma(\theta)}.
#' It uses default lags \eqn{\mathcal{L}=(0,1,2,4,5,20,50)}, in line with the paper specification.
#'
#' @param theta Vector of the parameters.
#' @param my Numeric matrix containing the time series of the multivariate system as columns.
#' @param type One among `exa`, `asy`, `cau` to determine which version of the model the user wants to estimate
#' among exact, asymptotic, and causal.
#' @param ts Logical. If TRUE, a time series of length n - 49 is returned, where n
#' is the number of rows in my. If FALSE, average values are returned.
#'
#' @returns Numeric matrix or vector, depending on the input to ts.
#' @export
errors_gmm <- function(theta, my, type = "exa", ts = FALSE,
                       lag_max = 5, lag_add = c(20, 50), delta = 1/252) {
  switch(type,
         "exa" = errors_exa_gmm(theta, my, ts, lag_max, lag_add, delta),
         "asy" = errors_asy_gmm(theta, my, ts, lag_max, lag_add, delta),
         "cau" = errors_cau_gmm(theta, my, ts, lag_max, lag_add, delta))
}


#' Jacobian of the moment conditions
#'
#' @description
#' This function calculates the Jacobian of the moment conditions, given by the difference
#' between empirical and theoretical cross-covariances, \eqn{\hat\gamma-\gamma(\theta)}.
#'
#'
#' @param theta Vector of the parameters.
#' @param my Numeric matrix containing the time series of the multivariate system as columns.
#' @param type One among `exa`, `asy`, `cau` to determine which version of the model the user wants to estimate
#' among exact, asymptotic, and causal.
#'
#' @returns Numeric matrix.
#' @export
jacobian_gmm <- function(theta, my, type = "exa",
                         lag_max = 5, lag_add = c(20, 50), delta = 1/252){
  switch(type,
         "exa" = jacobian_exa_gmm(theta, my, lag_max, lag_add, delta),
         "asy" = jacobian_asy_gmm(theta, my, lag_max, lag_add, delta),
         "cau" = jacobian_cau_gmm(theta, my, lag_max, lag_add, delta))
}


#' Converting parameters between list or vector form
#'
#' @description
#' Functions that convert the way in which parameters are stored, from vector to list
#' or vice-verse.
#'
#' @param lpar List of parameters.
#' @param vpar Vector of parameters.
#' @param d Dimensionality of the system of time-series.
#' @param type One among `exa`, `asy`, `cau` to determine which version of the model the user wants to estimate
#' among exact, asymptotic, and causal.
#'
#' @returns Numeric vector or List of numeric objects.
#' @export
#' @rdname par_transf
par_list2vec = function(lpar, d, type = "exa"){
  dd = d * (d - 1) / 2
  vpar <- switch(type,
                 "exa" = rep(NA, 3 * d + 2 * dd),
                 "asy" = rep(NA, 3 * d + 3 * dd),
                 "cau" = rep(NA, 3 * d + dd))
  if(type %in% c("exa", "cau")) lpar$univ["alpha", ] -> vpar[1 : d]
  if(type == "asy") diag(lpar$cov) -> vpar[1 : d]
  lpar$univ["nu", ] -> vpar[(d + 1) : (2 * d)]
  lpar$univ["H", ] -> vpar[(2 * d + 1) : (3 * d)]
  if(type == "exa"){
    lpar$rho[upper.tri(lpar$rho)] -> vpar[(3 * d + 1) : (3 * d + dd)]
    lpar$eta[upper.tri(lpar$eta)] -> vpar[(3 * d + dd + 1) : (3 * d + 2 * dd)]
  } else if (type == "asy") {
    lpar$cov[upper.tri(lpar$cov)] -> vpar[(3 * d + 1) : (3 * d + dd)]
    lpar$rho[upper.tri(lpar$rho)] -> vpar[(3 * d + dd + 1) : (3 * d + 2 * dd)]
    lpar$eta[upper.tri(lpar$eta)] -> vpar[(3 * d + 2 * dd + 1) : (3 * d + 3 * dd)]
  }
  return(vpar)
}

#' @export
#' @rdname par_transf
par_vec2list = function(vpar, d, type = "exa"){
  dd = d * (d - 1) / 2
  lpar <- list()
  if(type == "exa"){
    lpar$univ = matrix(nr = 3, nc = d)
    rownames(lpar$univ) <- c("H", "nu", "alpha")
    lpar$univ["alpha", ] <- vpar[1 : d]
    lpar$univ["nu", ] <- vpar[(d + 1) : (2 * d)]
    lpar$univ["H", ] <- vpar[(2 * d + 1) : (3 * d)]
    lpar$rho <- lpar$eta <- matrix(nr = d, nc = d)
    lpar$rho[upper.tri(lpar$rho)] <- vpar[(3 * d + 1) : (3 * d + dd)]
    diag(lpar$rho) <- 1
    lpar$rho[lower.tri(lpar$rho)] <- t(lpar$rho)[lower.tri(lpar$rho)]
    lpar$eta[upper.tri(lpar$eta)] <- vpar[(3 * d + dd + 1) : (3 * d + 2 * dd)]
    diag(lpar$eta) <- 0
    lpar$eta[lower.tri(lpar$eta)] <- - t(lpar$eta)[lower.tri(lpar$eta)]
  } else if(type == "asy") {
    lpar$univ = matrix(nr = 2, nc = d)
    rownames(lpar$univ) <- c("H", "nu")
    lpar$univ["nu", ] <- vpar[(d + 1) : (2 * d)]
    lpar$univ["H", ] <- vpar[(2 * d + 1) : (3 * d)]
    lpar$rho <- lpar$eta <- lpar$cov <- matrix(nr = d, nc = d)
    diag(lpar$cov) <- vpar[1 : d]
    lpar$cov[upper.tri(lpar$cov)] <- vpar[(3 * d + 1) : (3 * d + dd)]
    lpar$cov[lower.tri(lpar$cov)] <- t(lpar$cov)[lower.tri(lpar$cov)]
    lpar$rho[upper.tri(lpar$rho)] <- vpar[(3 * d + dd + 1) : (3 * d + 2 * dd)]
    diag(lpar$rho) <- 1
    lpar$rho[lower.tri(lpar$rho)] <- t(lpar$rho)[lower.tri(lpar$rho)]
    lpar$eta[upper.tri(lpar$eta)] <- vpar[(3 * d + 2 * dd + 1) : (3 * d + 3 * dd)]
    diag(lpar$eta) <- 0
    lpar$eta[lower.tri(lpar$eta)] <- - t(lpar$eta)[lower.tri(lpar$eta)]
  } else if(type == "cau"){
    lpar$univ = matrix(nr = 3, nc = d)
    rownames(lpar$univ) <- c("H", "nu", "alpha")
    lpar$univ["alpha", ] <- vpar[1 : d]
    lpar$univ["nu", ] <- vpar[(d + 1) : (2 * d)]
    lpar$univ["H", ] <- vpar[(2 * d + 1) : (3 * d)]
    lpar$rho <- matrix(nr = d, nc = d)
    lpar$rho[upper.tri(lpar$rho)] <- vpar[(3 * d + 1) : (3 * d + dd)]
    diag(lpar$rho) <- 1
    lpar$rho[lower.tri(lpar$rho)] <- t(lpar$rho)[lower.tri(lpar$rho)]
  }
  return(lpar)
}
