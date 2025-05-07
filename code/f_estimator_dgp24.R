C12 <-
  function(ai, aj, Hi, Hj, delta) {
    H = Hi + Hj
    integrand1 = function(u){u^(H - 2) * exp(- aj * u)}
    integrand2 = function(u){u^(H - 2) * (exp(ai * u) - exp(- aj * u))}
    value = 1 / (ai + aj) * (
      (exp((ai + aj) * delta) - 1) *
        integrate(f = integrand1, lower = delta, upper = Inf)[1]$value +
        integrate(f = integrand2, lower = 0 , upper = delta)[1]$value
    )
    return(value)
  }

C21 <-
  function (ai, aj, Hi, Hj, delta) {
    H = Hi + Hj
    integrand1 = function(u){u^(H - 2) * exp(- ai * u)}
    integrand2 = function(u){u^(H - 2) * (exp(aj * u) - exp(- ai * u))}
    value = 1/(ai + aj) * (
      (exp((ai + aj) * delta) - 1) *
        integrate(f = integrand1, lower = delta, upper = Inf)[1]$value +
        integrate(f = integrand2, lower = 0 , upper = delta)[1]$value
    )
    return(value)
  }

f_coherence <-
  function(rhoij, etaij, Hi, Hj){
    H = Hi + Hj
    if(H != 1){
      value =
        gamma(H + 1)^2 /
        (gamma(2 * Hi + 1) * gamma(2 * Hj + 1) * sin(pi * Hi) * sin(pi * Hj)) *
        (rhoij ^ 2 * sin(0.5 * pi * H) ^ 2 + etaij ^ 2 * cos(0.5 * pi * H) ^ 2)
    }
    else if (H == 1) {
      value =
        (gamma(2 * Hi + 1) * gamma(2 * Hj + 1) * sin(pi * Hi) * sin(pi * Hj)) ^ (- 1) *
        (rhoij ^ 2 + pi ^ 2 / 4 * etaij ^ 2)
    }
    return(value)
  }

f_est_ij <-
  function(xi, xj, ai, aj, nui, nuj, Hi, Hj,
           s = 1, delta = 1/252, constraint.out = FALSE) {
    H = Hi + Hj
    mx = na.omit(cbind(xi, xj))
    mx = mx[is.finite(rowSums(mx)), ]
    n = nrow(mx)
    c12s = C12(ai, aj, Hi, Hj, delta = delta * s)
    c21s = C21(ai, aj, Hi, Hj, delta = delta * s)
    g120 = cov(mx[, 1], mx[, 2])
    g12s = cov(mx[(s + 1) : n, 1], mx[1 : (n - s), 2])
    g21s = cov(mx[1 : (n - s), 1], mx[(s + 1) : n, 2])
    dComm = 1 / (H * (H - 1) * c12s * c21s * nui * nuj)
    rho = - dComm *
      (g120 * c12s - exp(aj * delta * s) * g21s * c12s + g120 * c21s -
         exp(ai * s * delta) * g12s * c21s)
    eta = dComm *
      (g120 * c12s - exp(aj * s * delta) * g21s * c12s -
         g120 * c21s + exp(ai * s * delta) * g12s * c21s)
    constr = f_coherence(rho, eta, Hi, Hj)
    if(constraint.out == TRUE) {
      lOut = list(
        "Estimates" = c("rhoij" = unname(rho), "etaij" = unname(eta)),
        "Constraint" = data.frame(Value = round(constr, 3),
                                  Admiss = constr < 1, row.names = FALSE))
    } else {
      lOut = c("rhoij" = unname(rho), "etaij" = unname(eta))
    }
    return(lOut)
  }

f_est_ij_asy <-
  function(xi, xj, nui, nuj, Hi, Hj,
           s = 1, delta = 1/252, constraint.out = FALSE) {
  H = Hi + Hj
  mx = na.omit(cbind(xi, xj))
  mx = mx[rowSums(mx) != Inf, ]
  mx = mx[rowSums(mx) != -Inf, ]
  n = nrow(mx)
  rho = (2 * cov(mx[, 1], mx[, 2]) -
           cov(mx[(1 + s) : n, 1], mx[1 : (n - s), 2]) -
           cov(mx[(1 + s) : n, 2], mx[1 : (n - s), 1])) /
    (nui * nuj * (s * delta) ^ H)
  eta = (cov(mx[(1 + s) : n, 2], mx[1 : (n - s), 1]) -
           cov(mx[(1 + s) : n, 1], mx[1 : (n - s), 2]) ) /
    (nui * nuj * (s * delta) ^ H)

  constr = f_coherence(Hi, Hj, rho, eta)
  if(constraint.out == TRUE) {
    lOut = list(
      "Estimates" = c("rhoij" = unname(rho), "etaij" = unname(eta)),
      "Constraint" = data.frame(Value = round(constr, 3),
                                Admiss = constr < 1, row.names = FALSE))
  } else if (constraint.out == FALSE) {
    lOut = c("rhoij" = unname(rho), "etaij" = unname(eta))
  }
  return(lOut)
  }

f_est_rho_eta <-
  function(mX, vA, vNu, vH,
           s = 1, delta = 1/252, constraint.out = FALSE){
    mRho <- mEta <-
      matrix(NA, nrow = ncol(mX), ncol = ncol(mX))
    colnames(mRho) <- colnames(mEta) <- colnames(mX)
    rownames(mRho) <- rownames(mEta) <- colnames(mX)
    for(i in 2 : ncol(mX)){
      for(j in 1 : (i - 1)){
        tmp = f_est_ij(mX[, i], mX[, j], vA[i], vA[j], vNu[i], vNu[j], vH[i], vH[j],
                       s = s, delta = delta, constraint.out = constraint.out)
        mRho[i, j] = tmp[1]
        mEta[i, j] = tmp[2]
        }
    }
    diag(mRho) <- 1
    diag(mEta) <- 0
    mRho[upper.tri(mRho, diag = FALSE)] = t(mRho)[upper.tri(mRho, diag = FALSE)]
    mEta[upper.tri(mEta, diag = FALSE)] = - t(mEta)[upper.tri(mEta, diag = FALSE)]
    return(list(
      "rho" = mRho,
      "eta" = mEta))
  }

f_est_rho_eta_asy <-
  function(mX, vNu, vH,
           s = 1, delta = 1/252, constraint.out = FALSE){
    mRho <- mEta <-
      matrix(NA, nrow = ncol(mX), ncol = ncol(mX))
    colnames(mRho) <- colnames(mEta) <- colnames(mX)
    rownames(mRho) <- rownames(mEta) <- colnames(mX)
    for(i in 2 : ncol(mX)){
      for(j in 1 : (i - 1)){
        tmp = f_est_ij_asy(mX[, i], mX[, j], vNu[i], vNu[j], vH[i], vH[j],
                           s = s, delta = delta, constraint.out = constraint.out)
        mRho[i, j] = tmp[1]
        mEta[i, j] = tmp[2]
      }
    }
    diag(mRho) <- 1
    diag(mEta) <- 0
    mRho[upper.tri(mRho, diag = FALSE)] = t(mRho)[upper.tri(mRho, diag = FALSE)]
    mEta[upper.tri(mEta, diag = FALSE)] = - t(mEta)[upper.tri(mEta, diag = FALSE)]
    return(list(
      "rho" = mRho,
      "eta" = mEta))
  }
