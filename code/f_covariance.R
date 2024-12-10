# empirical ccf

ccfl = function(x, y, l, mux = NULL, muy = NULL){
  mX = cbind(x, y)
  n = nrow(mX)
  if(is.null(mux) | is.null(muy)){
    mux = mean(x, na.rm = TRUE)
    muy = mean(y, na.rm = TRUE)
  }
  if(l > 0){
    mean(mX[(1 + l) : n, 1] * mX[1 : (n - l), 2] - mux * muy, na.rm = TRUE)
  } else if (l < 0) {
    l = - l
    mean(mX[(1 + l) : n, 2] * mX[1 : (n - l), 1] - mux * muy, na.rm = TRUE)
  } else if(l == 0){
    mean(mX[, 1] * mX[, 2] - mux * muy, na.rm = TRUE)
  }
}

# mfBm

f_var_fbm <- 
  function(s, H, sigma = 1){ 
    gamma = sigma ^ 2 * abs(s) ^ (2 * H)
    return(gamma)
  }

f_autocov_fbm <- 
  function(s, t, H, sigma = 1){ 
    gamma = sigma ^ 2 / 2 * (abs(s) ^ (2 * H) + abs(t) ^ (2 * H) - abs(t - s) ^ (2 * H))
    return(gamma)
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

f_eta_caus <- 
  function(rhoij, Hi, Hj){
  - rhoij * tan(pi/2 * (Hi + Hj)) * tan(pi/2 * (Hi - Hj))
}

f_crosscov_mfbm <- 
  function(s, t, Hi, Hj, rhoij, etaij, dDelta = 1, sigmai = 1, sigmaj = 1){ 
    if(Hi + Hj != 1){
      gamma = (sigmai * sigmaj)/2 * (
        (rhoij - etaij * sign(- s)) * abs(s) ^ (Hi + Hj) +
          (rhoij - etaij * sign(t)) * abs(t) ^ (Hi + Hj) -
          (rhoij - etaij * sign(t - s)) * abs(t - s) ^ (Hi + Hj))
        
    } else if (Hi + Hj == 1){
      gamma = (sigmai * sigmaj)/2 * 
        (rhoij * (abs(s) + abs(t) - abs(s - t)) +
           etaij * (t * log(abs(t)) - s * log(abs(s)) - (t - s) * log(abs(t - s)))
        )
    }
    return(gamma)
  }


# mfOU

f_var_fou <-
  function(a, nu, H){
  value = nu ^ 2 / (2 * a ^ (2 * H)) * gamma(1 + 2 * H)
  return(unname(value))
}

f_autocov_fou <-
  function(a, nu, H, delta){
    integrand = function(y) exp(- abs(y)) * abs(delta * a + y) ^ (2 * H)
    value = nu^2/(2*a^(2*H)) * 
      (0.5 * 
         integrate(f = integrand, lower = - Inf, upper = Inf)$value
       - abs(a * delta)^(2*H)
      )
    return(unname(value))
    
  }

f_autocov_asy_fou <-
  function(nu, H, var, delta){
  value = var - 0.5 * nu ^ 2 * delta ^ (2 * H)
  return(value)
}

f_acf_theo = 
  function(a, nu, H, lag.max = 50, delta = 1/252){
    vAcf =
      sapply(X = c(0 : lag.max) * delta, FUN = f_autocov_fou, 
             a = a, nu = nu,
             H = H)|> unname()
    
    names(vAcf) <- c(0 : lag.max)
    return(vAcf)
  }

f_acf_asy_theo = 
  function(nu, H, lag.max, var, delta = 1/252){
    vAcf =
      sapply(X = c(0 : lag.max) * delta, FUN = f_autocov_asy_fou, 
             nu =  nu, H = H, var = var) |> unname()
    names(vAcf) <- c(0 : lag.max)
    return(vAcf)
  }

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

f_cov_mfou = 
  function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij){
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

f_crosscov_ij_mfou = 
  function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij, delta) {
    H = Hi + Hj
    cov = f_cov_mfou(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij)
    if(H != 1) {
      
      value = exp(- ai * delta) * cov + nui * nuj * exp(- ai * delta) *
        H * (H - 1) * (rhoij + etaij) * 0.5 * C12(ai, aj, Hi, Hj, delta)
      
    } else if (H == 1){
      
      value = exp(- ai * delta) * cov - nui * nuj * exp(- ai * delta) * 
        etaij * 0.5 * C12(ai, aj, Hi, Hj, delta)
    }
    return(unname(value))
  }

f_crosscov_ji_mfou = 
  function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij, delta) {
    H = Hi + Hj
    cov = f_cov_mfou(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij)
    if(H != 1) {
      value = 
        exp(- aj * delta) * cov + nui * nuj * exp(- aj * delta) *
        H * (H - 1) * (rhoij - etaij) * 0.5 * C21(ai, aj, Hi, Hj, delta)
    } else if (H == 1){
      value = 
        exp(- aj * delta) * cov + nui * nuj * exp(- aj * delta) * 
        etaij * 0.5 * C21(ai, aj, Hi, Hj, delta)
    }
    return(unname(value))
  }

f_crosscov_asy_ij_mfou = 
  function(nui, nuj, Hi, Hj, rhoij, etaij, cov, delta) {
    value = 
      cov - 0.5 * nui * nuj * (rhoij + etaij) * delta ^ (Hi + Hj)
    return(value)
  }

f_crosscov_asy_ji_mfou = 
  function(nui, nuj, Hi, Hj, rhoij, etaij, cov, delta) {
    value = cov - 0.5 * nui * nuj * (rhoij - etaij) * delta ^ (Hi + Hj)
    return(value)
  }

f_ccf_theo = 
  function(ai, aj, nui, nuj, Hi, Hj, rhoij, etaij, lag.max = 50, delta = 1/252) {
    vCCF =
      c(
        sapply(X = c(lag.max : 1) * delta, 
               FUN = f_crosscov_ji_mfou, 
               ai = ai, aj = aj, 
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij),
        f_cov_mfou(ai = ai, aj = aj, 
                  nui =  nui, nuj = nuj,
                  Hi = Hi, Hj = Hj,
                  rhoij = rhoij, etaij = etaij),
        sapply(X = c(1 : lag.max) * delta, 
               FUN = f_crosscov_ij_mfou,
               ai = ai, aj = aj, 
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij)
      )|> unname()
    names(vCCF) <- c(- c(lag.max : 1), 0, c(1 : lag.max))
    return(vCCF)
  }

f_ccf_asy_theo = 
  function(nui, nuj, Hi, Hj, rhoij, etaij, cov, lag.max = 50, delta = 1/252) {
    vCCF =
      c(
        sapply(X = c(lag.max : 1) * delta, 
               FUN = f_crosscov_asy_ji_mfou, 
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij,
               cov = cov),
        cov,
        sapply(X = c(1 : lag.max) * delta, 
               FUN = f_crosscov_asy_ij_mfou,
               nui =  nui, nuj = nuj,
               Hi = Hi, Hj = Hj,
               rhoij = rhoij, etaij = etaij,
               cov = cov)
      ) |> unname()
    names(vCCF) <- c(- c(lag.max : 1), 0, c(1 : lag.max))
    return(vCCF)
  }
