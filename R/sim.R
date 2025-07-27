#' @export
#' @rdname sim_mfou
sim_mfou_exa <- function(n, H, rho, eta, a, nu, mu, delta, m = 1) {
  if(any(cov_existence_mfbm(H, rho, eta, out = TRUE) < 0 )) stop()
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions doesn't match.")

  syy = cov_dt_mfou(n = n, a = a, nu = nu, H = H, rho = rho, eta = eta, delta = delta)
  syy_lt = t(chol(syy))
  tmp = syy_lt %*% matrix(rnorm(n = n * d * m), nr = n * d, nc = m) + mu
  if(m == 1) {
    mfou = tmp[seq(1, length(tmp), d), ]
    for(i in 2 : d) {
      mfou = cbind(mfou, tmp[seq(i, length(tmp), d), ])
    }
    colnames(mfou) <- paste0("y", 1 : d)
  } else if (m > 1) {
    mfou = list()
    for(i in 1 : d) {
      mfou[[i]] <- tmp[seq(i, nrow(tmp), d), ]
    }
    names(mfou) <- paste0("y", 1 : d)
  }
  return(mfou)
}

#' @export
#' @rdname sim_mfou
sim_fou_exa <- function(n, H, a, nu, mu, delta, m = 1) {
  syy = cov_dt_fou(n = n, a = a, nu = nu, H = H, delta = delta)
  syy_lt = t(chol(syy))
  tmp = syy_lt %*% matrix(rnorm(n = n * m), nr = n, nc = m) + mu
  return(tmp)
}
