source("f_estimator_wxy23.R")
source("f_covariance.R")
source("f_estimator_dgp24.R")
source("f_covariance_derivatives.R")

f_init_cond <-  function(mX, delta = 1/252, type = "exa"){
  # Univariate estimates
  iNS = ncol(mX)
  vNamesUniv = colnames(mX)
  lUniv = list()
  for(i in 1 : iNS){
    lUniv[[i]] = unname(f_point_est_wxy(mX[, i], delta = delta))
  }
  vNu <- vH <- vMu <- vA <- vVar <- numeric()
  for(i in 1 : iNS){
    vVar[i] <- var(mX[is.finite(rowSums(mX)), i], na.rm = TRUE)
    vA[i] <- lUniv[[i]][3]
    vNu[i] <- lUniv[[i]][2]
    vH[i] <- lUniv[[i]][1]
    vMu[i] <- mean(mX[is.finite(rowSums(mX)), i], na.rm = TRUE)
  }
  names(vVar) <- paste0("var_", 1 : iNS)
  names(vNu) <- paste0("nu_", 1 : iNS)
  names(vA) <- paste0("a_", 1 : iNS)
  names(vH) <- paste0("H_", 1 : iNS)
  names(vMu) <- paste0("mu_", 1 : iNS)
  # Bivariate estimates
  mCombn = combn(x = 1 : ncol(mX), m = 2)
  iNC = ncol(mCombn)
  vNamesBiv = apply(X = combn(x = colnames(mX), m = 2), MARGIN = 2, FUN = paste0, collapse = "/")
  lBiv <- list()
  vRho <- vEta <- vCov <- numeric()
  if(type %in% c("exa", "cau")){
    for(i in 1 : iNC){
      j = mCombn[1, i]
      k = mCombn[2, i]
      lBiv[[i]] = 
        f_est_ij(xi = mX[, j], xj = mX[, k],
                 ai = lUniv[[j]][3], aj = lUniv[[k]][3],
                 Hi = lUniv[[j]][1], Hj = lUniv[[k]][1],
                 nui = lUniv[[j]][2], nuj = lUniv[[k]][2],
                 s = 1, delta = delta)
    }
    for(i in 1 : iNC){
      j = mCombn[1, i]
      k = mCombn[2, i]
      vCov[i] = cov(x = mX[is.finite(rowSums(mX)), j], y = mX[is.finite(rowSums(mX)), k], 
                    use = "complete.obs")
        vRho[i] = lBiv[[i]][1]
        vEta[i] = lBiv[[i]][2]
    }
    vTheta0 = 
      c(
        vA,
        vNu,
        vH,
        vRho,
        vEta
      )
    names(vTheta0) <-
      c(
        paste0("a_", vNamesUniv),
        paste0("nu_", vNamesUniv),
        paste0("H_", vNamesUniv),
        paste0("rho_", vNamesBiv),
        paste0("eta_", vNamesBiv)
      )
  } else if 
  (type == "asy"){
    for(i in 1 : iNC){
      j = mCombn[1, i]
      k = mCombn[2, i]
      lBiv[[i]] = 
        f_est_ij_asy(xi = mX[, j], xj = mX[, k],
                     Hi = lUniv[[j]][1], Hj = lUniv[[k]][1],
                     nui = lUniv[[j]][2], nuj = lUniv[[k]][2],
                     s = 1, delta = delta)
    }
    for(i in 1 : iNC){
      j = mCombn[1, i]
      k = mCombn[2, i]
      vCov[i] = cov(x = mX[is.finite(rowSums(mX)), j], y = mX[is.finite(rowSums(mX)), k], 
                    use = "complete.obs")
      vRho[i] = lBiv[[i]][1]
      vEta[i] = lBiv[[i]][2]
    }
    vTheta0 = 
      c(
        vVar,
        vNu,
        vH,
        vCov,
        vRho,
        vEta
      )
    names(vTheta0) <-
      c(
        paste0("var_", vNamesUniv),
        paste0("nu_", vNamesUniv),
        paste0("H_", vNamesUniv),
        paste0("cov_", vNamesBiv),
        paste0("rho_", vNamesBiv),
        paste0("eta_", vNamesBiv)
      )
  }
  
  if(type == "cau"){
    vTheta0 <- vTheta0[1 : (length(vTheta0) - iNC)]
  }
  return(vTheta0)
}
f_errors <- function(tet, mX, iLMax = 5, vLAdd = c(20, 50), delta = 1/252, causal = FALSE, ts = FALSE){
  
  iNS = ncol(mX)
  vA = tet[1 : iNS]
  vNu = tet[(iNS + 1) : (2 * iNS)]
  vH = tet[(2 * iNS + 1) : (3 * iNS)]
  # vMu = tet[(3 * iNS + 1) : (4 * iNS)]
  vMu = apply(mX, 2, mean, na.rm = TRUE)
  vNamesUniv = colnames(mX)
  names(vNu) <- names(vH) <- names(vA) <- names(vMu) <- vNamesUniv
  
  # Univariate moments
  lVar = list()
  lMean = list()
  lAutocov = list()
  for(i in 1 : ncol(mX)){
    n = length(mX[, i])
    lMean[[i]] = mX[, i] - vMu[i]
    lVar[[i]] = mX[, i] ^ 2 - vMu[i] ^ 2 - f_var_fou(a = vA[i], nu = vNu[i], H = vH[i])
    lAutocovTmp = list()
    if(iLMax != 0){
      for(j in 1 : iLMax){
        lAutocovTmp[[j]] = c(rep(NA, j), 
                             mX[(1 + j) : n, i] * mX[1 : (n - j), i] - 
                               vMu[i] ^ 2 - 
                               f_autocov_fou(a = vA[i], nu = vNu[i], H = vH[i], delta = j * delta)
        )
      }
    }
    if(!is.null(vLAdd)){
      for(j in 1 : length(vLAdd)){
        lAutocovTmp[[iLMax + j]] = c(rep(NA, vLAdd[j]), 
                                     mX[(1 + vLAdd[j]) : n, i] * mX[1 : (n - vLAdd[j]), i] - 
                                       vMu[i] ^ 2 - 
                                       f_autocov_fou(a = vA[i], nu = vNu[i], H = vH[i], delta = vLAdd[j] * delta)
        )
      }
    }
    lAutocov[[i]] = lAutocovTmp
  }
  names(lAutocov) <- vNamesUniv
  names(lVar) <- paste0(vNamesUniv, 0)
  names(lMean) <- paste0(vNamesUniv)
  lUnivMoms = c(lVar, do.call(c, lAutocov))
  rm(lMean, lVar, lAutocov)
  tUnivMoms = do.call(cbind, lUnivMoms)
  rm(lUnivMoms)
  if(iNS > 1){
    # Bivariate moments
    vNamesBiv = apply(X = combn(x = colnames(mX), m = 2), MARGIN = 2, FUN = paste0, collapse = "/")
    vNamesBivInv = apply(X = combn(x = colnames(mX), m = 2), MARGIN = 2, FUN = function(x) paste0(rev(x), collapse = "/"))
    mCombn = combn(x = 1 : ncol(mX), m = 2)
    iNC = ncol(mCombn)
    vRho = tet[(3 * iNS + 1) : (3 * iNS + iNC)]
    if(causal == FALSE){
      vEta = tet[(3 * iNS + iNC + 1) : (3 * iNS + 2 * iNC)]
    } else if (causal == TRUE){
      vEta = - vRho * tan(pi/2 * (vH[mCombn[1, ]] + vH[mCombn[2, ]])) * tan(pi/2 * (vH[mCombn[1, ]] - vH[mCombn[2, ]]))
    }
    names(vRho) <- names(vEta) <- vNamesBiv
    lCov = list()
    lCrosscov = list()
    for(i in 1 : iNC){
      s1 = mCombn[1, i]
      s2 = mCombn[2, i]
      x = mX[, c(s1, s2)]
      n = nrow(x)
      fCrosscov1_mom = function(x, l, delta = delta) {
        x[(1 + l) : n, 1] * x[1 : (n - l), 2] -
          vMu[s1] * vMu[s2] -
          f_crosscov_ij_mfou(ai = vA[s1], aj = vA[s2], nui = vNu[s1], nuj = vNu[s2],
                           Hi = vH[s1], Hj = vH[s2], rhoij = vRho[i], etaij = vEta[i],
                           delta = l * delta)
      }
      fCrosscov2_mom = function(x, l, delta = delta) {
        x[(1 + l) : n, 2] * x[1 : (n - l), 1] -
          vMu[s1] * vMu[s2] -
          f_crosscov_ji_mfou(ai = vA[s1], aj = vA[s2], nui = vNu[s1], nuj = vNu[s2],
                           Hi = vH[s1], Hj = vH[s2], rhoij = vRho[i], etaij = vEta[i],
                           delta = l * delta)
      }
      lCrosscov1 = list()
      lCrosscov2 = list()
      if(iLMax != 0){
        for(j in 1 : iLMax){
          lCrosscov1[[j]] = c(rep(NA, j), fCrosscov1_mom(x = x, l = j, delta = delta))
          lCrosscov2[[j]] = c(rep(NA, j), fCrosscov2_mom(x = x, l = j, delta = delta))
        }
      }
      if(!is.null(vLAdd)){
        for(j in 1 : length(vLAdd)){
          lCrosscov1[[iLMax + j]] = c(rep(NA, vLAdd[j]), fCrosscov1_mom(x = x, l = vLAdd[j], delta = delta))
          lCrosscov2[[iLMax + j]] = c(rep(NA, vLAdd[j]), fCrosscov2_mom(x = x, l = vLAdd[j], delta = delta))
        }
      }
      lFoo = c(lCrosscov1, lCrosscov2)
      rm(lCrosscov1, lCrosscov2)
      lCov[[i]] = x[, 1] * x[, 2] - vMu[s1] * vMu[s2] - 
        f_cov_mfou(ai = vA[s1], aj = vA[s2], nui = vNu[s1], nuj = vNu[s2],
                  Hi = vH[s1], Hj = vH[s2], rhoij = vRho[i], etaij = vEta[i])
      lCrosscov[[i]] = lFoo
      rm(lFoo)
    }
    names(lCov) = paste0(vNamesBiv, 0)
    tBivMoms = cbind(do.call(cbind, lCov), do.call(cbind, do.call(c, lCrosscov)))
    rm(lCov)
    mAll = cbind(tUnivMoms, tBivMoms)
    rm(tUnivMoms, tBivMoms)
  } else if(iNS == 1){
    mAll = tUnivMoms
    rm(tUnivMoms)
  }
  
  if(ts == FALSE){
    return(apply(mAll, 2, mean, na.rm = TRUE))
  } else if (ts == TRUE){
    return(mAll)
  }
}
f_errors_asy <- function(tet, mX, iLMax = 5, vLAdd = c(20, 50), delta = 1/252, causal = FALSE){
  
  iNS = ncol(mX)
  vVar = tet[1 : iNS]
  vNu = tet[(iNS + 1) : (2 * iNS)]
  vH = tet[(2 * iNS + 1) : (3 * iNS)]
  # vMu = tet[(3 * iNS + 1) : (4 * iNS)]
  vMu = apply(mX, 2, mean, na.rm = TRUE)
  vNamesUniv = colnames(mX)
  names(vNu) <- names(vH) <- names(vVar) <- names(vMu) <- vNamesUniv
  
  # Univariate moments
  lVar = list()
  lMean = list()
  lAutocov = list()
  for(i in 1 : ncol(mX)){ # think how to collapse params from list to vector
    n = length(mX[, i])
    lMean[[i]] = mX[, i] - vMu[i]
    lVar[[i]] = mX[, i] ^ 2 - vMu[i] ^ 2 - vVar[i]
    lAutocovTmp = list()
    if(iLMax != 0){
      for(j in 1 : iLMax){
        lAutocovTmp[[j]] = c(rep(NA, j), 
                             mX[(1 + j) : n, i] * mX[1 : (n - j), i] - 
                               vMu[i] ^ 2 - 
                               (vVar[i] - 0.5 * vNu[i] ^ 2 * (j * delta) ^ (2 * vH[i]))
        )
      }
    }
    if(!is.null(vLAdd)){
      for(j in 1 : length(vLAdd)){
        lAutocovTmp[[iLMax + j]] = c(rep(NA, vLAdd[j]), 
                                     mX[(1 + vLAdd[j]) : n, i] * mX[1 : (n - vLAdd[j]), i] - 
                                       vMu[i] ^ 2 - 
                                       (vVar[i] - 0.5 * vNu[i] ^ 2 * (vLAdd[j] * delta) ^ (2 * vH[i]))
        )
      }
    }
    names(lAutocovTmp) <- c(1 : iLMax, vLAdd)
    lAutocov[[i]] = lAutocovTmp
  }
  names(lAutocov) <- vNamesUniv
  names(lVar) <- paste0(vNamesUniv, 0)
  names(lMean) <- paste0(vNamesUniv)
  lUnivMoms = c(lVar, do.call(c, lAutocov))
  rm(lMean, lVar, lAutocov)
  tUnivMoms = do.call(cbind, lUnivMoms)
  rm(lUnivMoms)
  
  if(iNS > 1){
    ####################################################
    # Bivariate moments
    
    vNamesBiv = apply(X = combn(x = colnames(mX), m = 2), MARGIN = 2, FUN = paste0, collapse = "/")
    vNamesBivInv = apply(X = combn(x = colnames(mX), m = 2), MARGIN = 2, FUN = function(x) paste0(rev(x), collapse = "/"))
    mCombn = combn(x = 1 : ncol(mX), m = 2)
    iNC = ncol(mCombn)
    vCov = tet[(3 * iNS + 1) : (3 * iNS + iNC)]
    vRho = tet[(3 * iNS + iNC + 1) : (3 * iNS + 2 * iNC)]
    if(causal == FALSE){
      vEta = tet[(3 * iNS + 2 * iNC + 1) : (3 * iNS + 3 * iNC)]
    } else if (causal == TRUE){
      vEta = - vRho * tan(pi/2 * (vH[mCombn[1, ]] + vH[mCombn[2, ]])) * tan(pi/2 * (vH[mCombn[1, ]] - vH[mCombn[2, ]]))
    }
    names(vRho) <- names(vEta) <- vNamesBiv
    
    lCov = list()
    lCrosscov = list()
    for(i in 1 : iNC){
      s1 = mCombn[1, i]
      s2 = mCombn[2, i]
      x = mX[, c(s1, s2)]
      
      n = nrow(x)
      fCrosscov1_mom = function(x, l, delta = delta) {
        x[(1 + l) : n, 1] * x[1 : (n - l), 2] -
          vMu[s1] * vMu[s2] -
          (vCov[i] - (vRho[i] + vEta[i])/2 * vNu[s1] * vNu[s2] * (l * delta) ^ (vH[s1] + vH[s2]))
      }
      fCrosscov2_mom = function(x, l, delta = delta) {
        x[(1 + l) : n, 2] * x[1 : (n - l), 1] -
          vMu[s1] * vMu[s2] -
          (vCov[i] - (vRho[i] - vEta[i])/2 * vNu[s1] * vNu[s2] * (l * delta) ^ (vH[s1] + vH[s2]))
      }
      
      lCrosscov1 = list()
      lCrosscov2 = list()
      if(iLMax != 0){
        for(j in 1 : iLMax){
          lCrosscov1[[j]] = c(rep(NA, j), fCrosscov1_mom(x = x, l = j, delta = delta))
          lCrosscov2[[j]] = c(rep(NA, j), fCrosscov2_mom(x = x, l = j, delta = delta))
        }
      }
      if(!is.null(vLAdd)){
        for(j in 1 : length(vLAdd)){
          lCrosscov1[[iLMax + j]] = c(rep(NA, vLAdd[j]), fCrosscov1_mom(x = x, l = vLAdd[j], delta = delta))
          lCrosscov2[[iLMax + j]] = c(rep(NA, vLAdd[j]), fCrosscov2_mom(x = x, l = vLAdd[j], delta = delta))
        }
      }
      names(lCrosscov1) <- paste0(vNamesBiv[i], c(1 : iLMax, vLAdd))
      names(lCrosscov2) <- paste0(vNamesBivInv[i], c(1 : iLMax, vLAdd))
      lFoo = c(lCrosscov1, lCrosscov2)
      rm(lCrosscov1, lCrosscov2)
      
      lCov[[i]] = x[, 1] * x[, 2] - vMu[s1] * vMu[s2] - vCov[i]
      lCrosscov[[i]] = lFoo
      rm(lFoo)
    }
    
    names(lCov) = paste0(vNamesBiv, 0)
    tBivMoms = cbind(do.call(cbind, lCov), do.call(cbind, do.call(c, lCrosscov)))
    rm(lCov)
    
    mAll = cbind(tUnivMoms, tBivMoms)
    rm(tUnivMoms, tBivMoms)
  } else if(iNS == 1){
    mAll = tUnivMoms
    rm(tUnivMoms)
  }
  return(apply(mAll, 2, mean, na.rm = TRUE))
}
f_loss <-   function(tet, mX, iLMax = 5, vLAdd = c(20, 50), delta = 1/252, mW = NULL, type = "exa"){
  vAll <- 
    if(type == "exa"){
      f_errors(tet = tet, mX = mX, iLMax = iLMax, vLAdd = vLAdd, delta = delta, causal = FALSE)
    } else if(type == "asy"){
      f_errors_asy(tet = tet, mX = mX, iLMax = iLMax, vLAdd = vLAdd, delta = delta, causal = FALSE)
    } else if(type == "cau"){
      f_errors(tet = tet, mX = mX, iLMax = iLMax, vLAdd = vLAdd, delta = delta, causal = TRUE)
    }
  dLoss <- 
    if(is.null(mW)){
      sum(vAll ^ 2)
    } else {
      as.numeric(vAll %*% mW %*% vAll)
    }
  return(dLoss)
}
f_grad <-   function(tet, mX, iLMax = 5, vLAdd = c(20, 50), delta = 1/252, mW = NULL, type = "exa", jac = FALSE){
  iNS = ncol(mX)
  mCombn = combn(x = 1 : iNS, m = 2)
  iNC = ncol(mCombn)
  if(type == "exa"){
    vA = tet[1 : iNS]
    vNu = tet[(iNS + 1) : (2 * iNS)]
    vH = tet[(2 * iNS + 1) : (3 * iNS)]
    # vMu = tet[(3 * iNS + 1) : (4 * iNS)]
    vEmpty = rep(0, 3 * iNS + 2 * iNC)
    if(iNC == 0){
      names(vEmpty) <- 
        c(
          paste0("a", 1 : iNS),
          paste0("nu", 1 : iNS),
          paste0("H", 1 : iNS)
          # paste0("mu", 1 : iNS),
        )
    } else {
      vRho = tet[(3 * iNS + 1) : (3 * iNS + iNC)]
      vEta = tet[(3 * iNS + iNC + 1) : (3 * iNS + 2 * iNC)]
      names(vEmpty) <- 
        c(
          paste0("a", 1 : iNS),
          paste0("nu", 1 : iNS),
          paste0("H", 1 : iNS),
          # paste0("mu", 1 : iNS),
          paste0("rho", 1 : iNC),
          paste0("eta", 1 : iNC)
        ) 
    }
    fVEC = function(vV, iE, c){
      vVC = vV
      vVC[iE] = c
      return(vVC)
    }

    # Univariate moments
    lMean = list()
    lVar = list()
    lAutocov = list()
    for(i in 1 : iNS){ 
      # lMean[[i]] = fVEC(vV = vEmpty, iE = 3 * iNS + i, c = 1)
      lVar[[i]] = fVEC(vV = vEmpty, 
                       iE = c(i, # a
                              iNS + i, # nu
                              2 * iNS + i # H
                              #, 3 * iNS + i # mu
                       ), 
                       c = c(dVda(a = vA[i], v = vNu[i], H = vH[i]),
                             dVdv(a = vA[i], v = vNu[i], H = vH[i]),
                             dVdH(a = vA[i], v = vNu[i], H = vH[i])
                             # ,2 * vMu[i]
                       )
      )
      lAutocovTmp = list()
      if(iLMax != 0){
        for(j in 1 : iLMax){
          lAutocovTmp[[j]] = fVEC(vV = vEmpty, 
                                  iE = c(i, # a
                                         iNS + i, # nu
                                         2 * iNS + i # H
                                         # , 3 * iNS + i
                                  ), # mu 
                                  c = 
                                    c(
                                      dGisda(a = vA[i], v = vNu[i], H = vH[i], delta = j * delta),
                                      dGisdv(a = vA[i], v = vNu[i], H = vH[i], delta = j * delta),
                                      dGisdH(a = vA[i], v = vNu[i], H = vH[i], delta = j * delta)
                                      # , 2 * vMu[i]
                                    )
          )
        }
      }
      if(!is.null(vLAdd)){
        for(j in 1 : length(vLAdd)){
          lAutocovTmp[[iLMax + j]] = 
            fVEC(vV = vEmpty, 
                 iE = c(i, # a
                        iNS + i, # nu 
                        2 * iNS + i # H
                        # ,3 * iNS + i
                 ), # mu
                 c = 
                   c(
                     dGisda(a = vA[i], v = vNu[i], H = vH[i], delta = vLAdd[j] * delta),
                     dGisdv(a = vA[i], v = vNu[i], H = vH[i], delta = vLAdd[j] * delta),
                     dGisdH(a = vA[i], v = vNu[i], H = vH[i], delta = vLAdd[j] * delta)
                     # , 2 * vMu[i]
                   )
            )
        }
      }
      lAutocov[[i]] = lAutocovTmp
    }
    
    # Bivariate moments
    mCombn = combn(x = 1 : iNS, m = 2)
    iNC = ncol(mCombn)
    lCov = list()
    lCrosscov = list()
    if(iNC > 0){
      for(i in 1 : iNC){
        s1 = mCombn[1, i]
        s2 = mCombn[2, i]
        
        lCov[[i]] = 
          fVEC(vV = vEmpty, 
               iE = c(s1, s2, # a
                      iNS + s1, iNS + s2, # nu
                      2 * iNS + s1, 2 * iNS + s2, # H
                      # 3 * iNS + s1, 3 * iNS + s2, # mu
                      3 * iNS + i, # rho
                      3 * iNS + iNC + i), # eta 
               c = 
                 c(
                   dG0da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                         rho = vRho[i], eta = vEta[i]),
                   dG0dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                         rho = vRho[i], eta = vEta[i]),
                   # vMu[s2],
                   # vMu[s1],
                   dG0drho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                           rho = vRho[i], eta = vEta[i]),
                   dG0deta(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                           rho = vRho[i], eta = vEta[i])
                 )
          )
        lCrosscov1 = list()
        lCrosscov2 = list()
        if(iLMax != 0){
          for(j in 1 : iLMax){
            lCrosscov1[[j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i, # rho
                          3 * iNS + iNC + i), # eta,
                   c = 
                     c(
                       dGs12da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = j * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs12drho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12deta(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = j * delta)
                     )
              )
            lCrosscov2[[j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i, # rho
                          3 * iNS + iNC + i), # eta,
                   c = 
                     c(
                       dGs21da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = j * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs21drho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21deta(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = j * delta)
                     )
              )
          }
        }
        if(!is.null(vLAdd)){
          for(j in 1 : length(vLAdd)){
            lCrosscov1[[iLMax + j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i, # rho
                          3 * iNS + iNC + i), # eta,
                   c = 
                     c(
                       dGs12da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs12drho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12deta(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta)
                     )
              )
            lCrosscov2[[iLMax + j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i, # rho
                          3 * iNS + iNC + i), # eta,
                   c = 
                     c(
                       dGs21da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                               rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs21drho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21deta(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta)
                     )
              )
          }
        }
        lFoo = c(lCrosscov1, lCrosscov2)
        rm(lCrosscov1, lCrosscov2)
        lCrosscov[[i]] = lFoo
      }
    }
    
    # Aggregate
    mOut = # minus sign because moment are: empirical - theoretical
      - rbind(
        # means (x1 - mu1)
        # do.call(rbind, lMean),
        # vars
        do.call(rbind, lVar),
        # autocov
        do.call(rbind, do.call(c, lAutocov)),
        # cov
        do.call(rbind, lCov),
        # ccf
        if(iNC > 0){
          do.call(rbind, do.call(c, lCrosscov))
        }
      )
    colnames(mOut) <- names(vEmpty)
    vLoss = f_errors(tet = tet, mX = mX, iLMax = iLMax, vLAdd = vLAdd,
                    delta = delta, causal = FALSE)
  } else if 
  (type == "asy"){
    vVar = tet[1 : iNS]
    vNu = tet[(iNS + 1) : (2 * iNS)]
    vH = tet[(2 * iNS + 1) : (3 * iNS)]
    vEmpty = rep(0, 3 * iNS + 3 * iNC)
    if(iNC == 0){
      names(vEmpty) <- 
        c(
          paste0("var", 1 : iNS),
          paste0("nu", 1 : iNS),
          paste0("H", 1 : iNS)
        )
    } else {
      vCov = tet[(3 * iNS + 1) : (3 * iNS + iNC)]
      vRho = tet[(3 * iNS + iNC + 1) : (3 * iNS + 2 * iNC)]
      vEta = tet[(3 * iNS + 2 * iNC + 1) : (3 * iNS + 3 * iNC)]
      names(vEmpty) <- 
        c(
          paste0("var", 1 : iNS),
          paste0("nu", 1 : iNS),
          paste0("H", 1 : iNS),
          paste0("cov", 1 : iNC),
          paste0("rho", 1 : iNC),
          paste0("eta", 1 : iNC)
        ) 
    }
    fVEC = function(vV, iE, c){
      vVC = vV
      vVC[iE] = c
      return(vVC)
    }
    
    # Univariate moments
    lVar = list()
    lAutocov = list()
    for(i in 1 : iNS){ 
      lVar[[i]] = fVEC(vV = vEmpty, 
                       iE = c(i # var
                       ), 
                       c = c(1))
      lAutocovTmp = list()
      if(iLMax != 0){
        for(j in 1 : iLMax){
          lAutocovTmp[[j]] = fVEC(vV = vEmpty, 
                                  iE = c(i, # var
                                         iNS + i, # nu
                                         2 * iNS + i # H
                                  ),  
                                  c = 
                                    c(
                                      1,
                                      - vNu[i] * (j * delta) ^ (2 * vH[i]),
                                      - vNu[i] ^ 2 * (j * delta) ^ (2 * vH[i]) * log(j * delta)
                                    )
          )
        }
      }
      if(!is.null(vLAdd)){
        for(j in 1 : length(vLAdd)){
          lAutocovTmp[[iLMax + j]] = 
            fVEC(vV = vEmpty, 
                 iE = c(i, # var
                        iNS + i, # nu
                        2 * iNS + i # H
                 ),  
                 c = 
                   c(
                     1,
                     - vNu[i] * (vLAdd[j] * delta) ^ (2 * vH[i]),
                     - vNu[i] ^ 2 * (vLAdd[j] * delta) ^ (2 * vH[i]) * log(vLAdd[j] * delta)
                   )
            )
        }
      }
      lAutocov[[i]] = lAutocovTmp
    }
    
    # Bivariate moments
    mCombn = combn(x = 1 : iNS, m = 2)
    iNC = ncol(mCombn)
    lCov = list()
    lCrosscov = list()
    if(iNC > 0){
      for(i in 1 : iNC){
        s1 = mCombn[1, i]
        s2 = mCombn[2, i]
        lCov[[i]] = 
          fVEC(vV = vEmpty, 
               iE = c(3 * iNS + i), 
               c = 
                 c(1)
          )
        lCrosscov1 = list()
        lCrosscov2 = list()
        if(iLMax != 0){
          for(j in 1 : iLMax){
            lCrosscov1[[j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          3 * iNS + i, # cov
                          3 * iNS + iNC + i, # rho
                          3 * iNS + 2 * iNC + i # eta
                   ),
                   c = 
                     c(
                       - (vRho[i] + vEta[i])/2 * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] + vEta[i])/2 * vNu[s1] * (j * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] + vEta[i])/2 * vNu[s1] * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2]) * log(j * delta),
                       - (vRho[i] + vEta[i])/2 * vNu[s1] * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2]) * log(j * delta),
                       1,
                       - 0.5 * vNu[s1] * vNu[s2] * (j *delta) ^ (vH[s1] + vH[s2]),
                       - 0.5 * vNu[s1] * vNu[s2] * (j *delta) ^ (vH[s1] + vH[s2])))
            
            lCrosscov2[[j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          3 * iNS + i, # cov
                          3 * iNS + iNC + i, # rho
                          3 * iNS + 2 * iNC + i # eta
                   ),
                   c = 
                     c(
                       - (vRho[i] - vEta[i])/2 * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] - vEta[i])/2 * vNu[s1] * (j * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] - vEta[i])/2 * vNu[s1] * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2]) * log(j * delta),
                       - (vRho[i] - vEta[i])/2 * vNu[s1] * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2]) * log(j * delta),
                       1,
                       - 0.5 * vNu[s1] * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2]),
                       0.5 * vNu[s1] * vNu[s2] * (j * delta) ^ (vH[s1] + vH[s2])
                     )
              )
          }
        }
        if(!is.null(vLAdd)){
          for(j in 1 : length(vLAdd)){
            lCrosscov1[[iLMax + j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          3 * iNS + i, # cov
                          3 * iNS + iNC + i, # rho
                          3 * iNS + 2 * iNC + i # eta
                   ),
                   c = 
                     c(
                       - (vRho[i] + vEta[i])/2 * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] + vEta[i])/2 * vNu[s1] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] + vEta[i])/2 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]) * log(vLAdd[j] * delta),
                       - (vRho[i] + vEta[i])/2 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]) * log(vLAdd[j] * delta),
                       1,
                       - 0.5 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]),
                       - 0.5 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2])))
            
            lCrosscov2[[iLMax + j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          3 * iNS + i, # cov
                          3 * iNS + iNC + i, # rho
                          3 * iNS + 2 * iNC + i # eta
                   ),
                   c = 
                     c(
                       - (vRho[i] - vEta[i])/2 * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] - vEta[i])/2 * vNu[s1] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]),
                       - (vRho[i] - vEta[i])/2 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]) * log(vLAdd[j] * delta),
                       - (vRho[i] - vEta[i])/2 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]) * log(vLAdd[j] * delta),
                       1,
                       - 0.5 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2]),
                       0.5 * vNu[s1] * vNu[s2] * (vLAdd[j] * delta) ^ (vH[s1] + vH[s2])
                     )
              )
          }
        }
        lFoo = c(lCrosscov1, lCrosscov2)
        rm(lCrosscov1, lCrosscov2)
        lCrosscov[[i]] = lFoo
      }
    }
    
    # Aggregate
    mOut = # minus sign because moment are: empirical - theoretical
      - rbind(
        # means (x1 - mu1)
        # do.call(rbind, lMean),
        # vars
        do.call(rbind, lVar),
        # autocov
        do.call(rbind, do.call(c, lAutocov)),
        # cov
        do.call(rbind, lCov),
        # ccf
        if(iNC > 0){
          do.call(rbind, do.call(c, lCrosscov))
        }
      )
    vLoss = f_errors_asy(tet = tet, mX = mX, iLMax = iLMax, vLAdd = vLAdd,
                       delta = delta, causal = FALSE)
  } else if 
  (type == "cau"){
    vA = tet[1 : iNS]
    vNu = tet[(iNS + 1) : (2 * iNS)]
    vH = tet[(2 * iNS + 1) : (3 * iNS)]
    # vMu = tet[(3 * iNS + 1) : (4 * iNS)]
    if(iNS == 1){
      vEmpty = rep(0, 3 * iNS)
      names(vEmpty) <- 
        c(
          paste0("a", 1 : iNS),
          paste0("nu", 1 : iNS),
          paste0("H", 1 : iNS)
          # paste0("mu", 1 : iNS),
        )
    } else {
      mCombn = combn(x = 1 : iNS, m = 2)
      iNC = ncol(mCombn)
      vEmpty = rep(0, 3 * iNS + iNC)
      
      vRho = tet[(3 * iNS + 1) : (3 * iNS + iNC)]
      vEta = - vRho *
        tan(pi/2 * (vH[mCombn[1, ]] + vH[mCombn[2, ]])) *
        tan(pi/2 * (vH[mCombn[1, ]] - vH[mCombn[2, ]]))
      
      names(vEmpty) <- 
        c(
          paste0("a", 1 : iNS),
          paste0("nu", 1 : iNS),
          paste0("H", 1 : iNS),
          # paste0("mu", 1 : iNS),
          paste0("rho", 1 : iNC)
        ) 
    }
    fVEC = function(vV, iE, c){
      vVC = vV
      vVC[iE] = c
      return(vVC)
    }
    # Univariate moments
    lMean = list()
    lVar = list()
    lAutocov = list()
    for(i in 1 : iNS){ 
      # lMean[[i]] = fVEC(vV = vEmpty, iE = 3 * iNS + i, c = 1)
      lVar[[i]] = fVEC(vV = vEmpty, 
                       iE = c(i, # a
                              iNS + i, # nu
                              2 * iNS + i # H
                              #, 3 * iNS + i # mu
                       ), 
                       c = c(dVda(a = vA[i], v = vNu[i], H = vH[i]),
                             dVdv(a = vA[i], v = vNu[i], H = vH[i]),
                             dVdH(a = vA[i], v = vNu[i], H = vH[i])
                             # ,2 * vMu[i]
                       )
      )
      lAutocovTmp = list()
      if(iLMax != 0){
        for(j in 1 : iLMax){
          lAutocovTmp[[j]] = fVEC(vV = vEmpty, 
                                  iE = c(i, # a
                                         iNS + i, # nu
                                         2 * iNS + i # H
                                         # , 3 * iNS + i
                                  ), # mu 
                                  c = 
                                    c(
                                      dGisda(a = vA[i], v = vNu[i], H = vH[i], delta = j * delta),
                                      dGisdv(a = vA[i], v = vNu[i], H = vH[i], delta = j * delta),
                                      dGisdH(a = vA[i], v = vNu[i], H = vH[i], delta = j * delta)
                                      # , 2 * vMu[i]
                                    )
          )
        }
      }
      if(!is.null(vLAdd)){
        for(j in 1 : length(vLAdd)){
          lAutocovTmp[[iLMax + j]] = 
            fVEC(vV = vEmpty, 
                 iE = c(i, # a
                        iNS + i, # nu 
                        2 * iNS + i # H
                        # ,3 * iNS + i
                 ), # mu
                 c = 
                   c(
                     dGisda(a = vA[i], v = vNu[i], H = vH[i], delta = vLAdd[j] * delta),
                     dGisdv(a = vA[i], v = vNu[i], H = vH[i], delta = vLAdd[j] * delta),
                     dGisdH(a = vA[i], v = vNu[i], H = vH[i], delta = vLAdd[j] * delta)
                     # , 2 * vMu[i]
                   )
            )
        }
      }
      lAutocov[[i]] = lAutocovTmp
    }
    # Bivariate moments
    lCov = list()
    lCrosscov = list()
    if(iNC > 0){
      for(i in 1 : iNC){
        s1 = mCombn[1, i]
        s2 = mCombn[2, i]
        lCov[[i]] = 
          fVEC(vV = vEmpty, 
               iE = c(s1, s2, # a
                      iNS + s1, iNS + s2, # nu
                      2 * iNS + s1, 2 * iNS + s2, # H
                      # 3 * iNS + s1, 3 * iNS + s2, # mu
                      3 * iNS + i), # rho 
               c = 
                 c(
                   dG0da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                          rho = vRho[i], eta = vEta[i]),
                   dG0dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                         rho = vRho[i], eta = vEta[i]),
                   dG0dH(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                         rho = vRho[i], eta = vEta[i]),
                   # vMu[s2],
                   # vMu[s1],
                   dG0cdrho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2])
                 )
          )
        lCrosscov1 = list()
        lCrosscov2 = list()
        if(iLMax != 0){
          for(j in 1 : iLMax){
            lCrosscov1[[j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i), # rho,
                   c = 
                     c(
                       dGs12da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs12cdH1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = j * delta),
                       dGs12cdH2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = j * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs12cdrho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                  delta = j * delta)
                     )
              )
            lCrosscov2[[j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i), # rho,
                   c = 
                     c(
                       dGs21da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = j * delta),
                       dGs21cdH1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = j * delta),
                       dGs21cdH2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = j * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs21cdrho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                  delta = j * delta)
                     )
              )
          }
        }
        if(!is.null(vLAdd)){
          for(j in 1 : length(vLAdd)){
            lCrosscov1[[iLMax + j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i), # rho,
                   c = 
                     c(
                       dGs12da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs12cdH1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = vLAdd[j] * delta),
                       dGs12cdH2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = vLAdd[j] * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs12cdrho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                  delta = vLAdd[j] * delta)
                     )
              )
            lCrosscov2[[iLMax + j]] = 
              fVEC(vV = vEmpty, 
                   iE = c(s1, s2, # a
                          iNS + s1, iNS + s2, # nu
                          2 * iNS + s1, 2 * iNS + s2, # H
                          # 3 * iNS + s1, 3 * iNS + s2, # mu
                          3 * iNS + i), # rho,
                   c = 
                     c(
                       dGs21da1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21da2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21dv1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21dv2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                rho = vRho[i], eta = vEta[i], delta = vLAdd[j] * delta),
                       dGs21cdH1(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = vLAdd[j] * delta),
                       dGs21cdH2(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                 rho = vRho[i], delta = vLAdd[j] * delta),
                       # vMu[s2],
                       # vMu[s1],
                       dGs21cdrho(a1 = vA[s1], a2 = vA[s2], v1 = vNu[s1], v2 = vNu[s2], H1 = vH[s1], H2 = vH[s2],
                                  delta = vLAdd[j] * delta)
                     )
              )
          }
        }
        lFoo = c(lCrosscov1, lCrosscov2)
        rm(lCrosscov1, lCrosscov2)
        lCrosscov[[i]] = lFoo
      }
    }
    # Aggregate
    mOut = # minus sign because moment are: empirical - theoretical
      - rbind(
        # means (x1 - mu1)
        # do.call(rbind, lMean),
        # vars
        do.call(rbind, lVar),
        # autocov
        do.call(rbind, do.call(c, lAutocov)),
        # cov
        do.call(rbind, lCov),
        # ccf
        if(iNC > 0){
          do.call(rbind, do.call(c, lCrosscov))
        }
      )
    colnames(mOut) <- names(vEmpty)
    vLoss = f_errors(tet = tet, mX = mX, iLMax = iLMax, vLAdd = vLAdd,
                    delta = delta, causal = TRUE)
  }
  if(is.null(mW)){
    mW = diag(nrow = nrow(mOut))
  }
  
  vGrad = 2 * t(mOut) %*% mW %*% vLoss
  if(jac == FALSE){
    return(vGrad)
  } else if (jac == TRUE){
    lOut = list(gradient = vGrad,
                jacobian_gamma = mOut)
    return(lOut)
  }
}
f_mde <-       function(mX, iLMax = 5, vLAdd = c(20, 50), delta = 1/252, mW = NULL, type = "exa"){
  
  if(type != "exa" & type != "asy" & type != "cau"){
    stop("type must be \"exa\", \"asy\", or \"cau\"")
  }
  
  iNS = ncol(mX)
  iNC = iNS * (iNS - 1) / 2
  
  if(is.null(mW)){
    iNL = iLMax + length(vLAdd)
    mW = diag(nrow = iNS * (iNL + 1) + iNS * (iNS - 1) / 2 * (2 * iNL + 1))
  }
  
  if(type == "exa"){
    vLower = c(rep(1e-4, iNS), rep(1e-2, iNS),  rep(1e-3, iNS),     rep(- 1, iNC), rep(- 10, iNC))
    vUpper = c(rep(10, iNS),  rep(5,  iNS),  rep(1 - 1e-3, iNS), rep(1, iNC), rep(  10, iNC))
  } else if (type == "asy"){
    vLower = c(rep(1e-3, iNS), rep(1e-3, iNS),  rep(1e-3, iNS),  rep(- Inf, iNC), rep(- 1, iNC), rep(- 10, iNC))
    vUpper = c(rep(10, iNS),  rep(10,  iNS),  rep(1 - 1e-3, iNS), rep(Inf, iNC), rep(1, iNC), rep(  10, iNC))
  } else if (type == "cau"){
    vLower = c(rep(1e-4, iNS), rep(1e-2, iNS),  rep(1e-3, iNS),     rep(- 1, iNC))
    vUpper = c(rep(10, iNS),  rep(5,  iNS),  rep(1 - 1e-3, iNS), rep(1, iNC))
  }
  
  oOpt =
    optim(
      par = f_init_cond(mX = mX, delta = delta, type = type),
      fn = f_loss,
      gr = f_grad,
      mX = mX,
      type = type,
      delta = delta,
      iL = iLMax,
      vLAdd = vLAdd,
      mW = mW,
      method = "L-BFGS-B",
      lower = vLower,  
      upper = vUpper,
      control = list(maxit = 3000, trace = 6)
    )
  return(oOpt$par)
}
