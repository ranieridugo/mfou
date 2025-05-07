# point estimates

f_est_H <- 
  function(y, verbose = FALSE){ 
    x = na.omit(y)
    x = x[is.finite(x)]
    n = length(x)
    num = sum((x[5 : n] - 2 * x[3 : (n - 2)] + x[1 : (n - 4)]) ^ 2, na.rm = TRUE) 
    den = sum((x[3 : n] - 2 * x[2 : (n - 1)] + x[1 : (n - 2)]) ^ 2, na.rm = TRUE)
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

f_est_nu <- 
  function(y, H, delta = 1/252){
    x = na.omit(y)
    x = x[is.finite(x)]
    n = length(x)
    num = sum((x[3:n] - 2 * head(x[2:n], -1) + head(x, -2))^2)
    den = (n * (4 - 2 ^ (2 * H)) * delta ^ (2 * H))
    nu = sqrt(num / den)
    return(nu)
  }

f_est_mu <-
  function(y){
  x = na.omit(y)
  x = x[is.finite(x)]
  n = length(x)
  return(sum(x) / n)
}

f_est_alpha <-
  function(y, H, nu){
  x = na.omit(y)
  x = x[is.finite(x)]
  n = length(x)
  num = n * sum(x ^ 2) - sum(x) ^ 2 
  den = n ^ 2 * nu ^ 2 * H * gamma(2 * H)
  alpha = (num / den) ^ (- 1 / (2 * H))
  return(alpha)
}

# standard errors

SS <-
  function(H, j = 1E5){
  rhoj = function(j, H){
    1 / (2 * (4 - 2^(2*H))) *
      (- abs(j + 2) ^ (2 * H) + 4 * abs(j + 1) ^ (2 * H) - 6 *
         abs(j) ^ (2 * H) + 4 * abs(j - 1) ^ (2 * H) - abs(j - 2) ^ (2 * H))
  }
  
  j = c(1:j)
  
  dS11 = 2 + 
    2 ^ (2 - 4 * H) *
    sum((rhoj(H = H, j = j+2) + 4 * rhoj(H = H, j = j+1) + 6 * rhoj(H = H, j = j) + 
           4 * rhoj(H = H, j = abs(j - 1)) + rhoj(H = H, j = abs(j - 2))) ^ 2)
  
  dS12 = 2 ^ (1 - 2 * H) *
    (
      4 * (rhoj(H = H, j = 1) + 1) ^ 2 +
        2 * sum((rhoj(H = H, j = j + 1) + 2 * rhoj(H = H, j = j) +
                   rhoj(H = H, j = j - 1)) ^ 2)
    )
  
  dS22 = 2 + 4 * sum(rhoj(H = H, j = j)^2)
  
  return(c("S11" = dS11, "S12" = dS12, "S22" = dS22))}

f_se_H <- 
  function(x = NULL, H = NULL, n = NULL){
    if(!is.null(x)){
      H = H.hat(y = x)
      n = length(na.omit(x))
    }
    SSx = SS(H = H)
    num = SSx[1] + SSx[3] - 2 * SSx[2]
    den = ((2 * log(2)) ^ 2) * n
    var = num / den
    return(sqrt(var))
  }

f_se_nu <- 
  function(nu, H, n, delta){
    var = H.se(H = H, n = n) ^ 2 * nu ^ 2 * log(1 / delta) ^ 2
    return(sqrt(var))
  }

f_se_mu <- 
  function(nu, alpha, H, delta = 1/252, n){
    var = nu ^ 2 / alpha ^ 2 * (n * delta) ^ (2 * H - 2)
    return(sqrt(var))
  }

f_se_alpha <-
  function(H, alpha, delta, n){
  TT = delta * n
  if(H > 0 & H < 3/4){
    if(H < 1/2){
      phi = 1/(4 * H ^ 2) * ((4*H - 1) + (2*gamma(2 - 4 * H) * gamma(4 * H)) /
                               (gamma(2 * H) * gamma(1 - 2 * H))) 
    } else if (H >= 1/2) {
      phi = (4*H - 1)/(4 * H ^ 2) * (1 + (gamma(3 - 4 * H) * gamma(4 * H - 1))/
                                       (gamma(2 - 2 * H) * gamma(2 * H)))
    }
    var = phi * alpha / TT
  } else if (H == 3/4) {
    var = 16*alpha/(9 * pi) * log(TT)^2/TT
  } else {var = "Error Rosenblatt"}
  return(sqrt(var))
}

# full estimates

f_point_est_wxy <-
  function(y, delta = 1/252, H = NULL){
    n = length(na.omit(y))
    y[is.finite(y)]
    if(!is.null(H)){
      H.p = H
    } else (
      H.p = f_est_H(y = y)
    )
    nu.p = f_est_nu(y = y, H = H.p, delta = delta)
    alpha.p = f_est_alpha(y = y, H = H.p, nu = nu.p)
    mu.p = f_est_mu(y = y)
    out = c("H" = H.p, "v" = nu.p, "a" = alpha.p, "mu" = mu.p)
    return(out)
  }

f_confint <-
  function(mean, se, conf = 0.05){
    lower = mean - se * qnorm(conf / 2)
    upper = mean + se * qnorm(conf / 2)
    return(c(lower, upper))
  } 

f_est_wxy =
  function(y, delta = 1/252){
    n = length(na.omit(y))
    
    H.p = f_est_H(y = y)
    H.stde = f_se_H(x = NULL, H = H.p, n = n)
    H.ci = f_confint(mean = H.p, se = H.stde)
    
    nu.p = f_est_nu(y = y, H = H.p, delta = delta)
    nu.stde = f_se_nu(nu = nu.p, H = H.p, n = n, delta = delta)
    nu.ci = f_confint(mean = nu.p, se = nu.stde)
    alpha.p = f_est_alpha(y = y, H = H.p, nu = nu.p)
    alpha.stde = f_se_alpha(H = H.p, alpha = alpha.p, delta = delta, n = n)
    alpha.ci = f_confint(mean = alpha.p, se = alpha.stde)
    
    mu.p = f_est_mu(y = y)
    mu.stde = f_se_mu(nu = nu.p, alpha = alpha.p, H = H.p, n = n)
    mu.ci = f_confint(mean = mu.p, se = mu.stde)
    
    out = data.frame(H = c(H.p, H.stde, H.ci[1], H.ci[2]),
                     nu = c(nu.p, nu.stde, nu.ci[1], nu.ci[2]),
                     mu = c(mu.p, mu.stde, mu.ci[1], mu.ci[2]),
                     alpha = c(alpha.p, alpha.stde, alpha.ci[1], alpha.ci[2]),
                     row.names = c("p.est", "std.err", "lower", "upper"))
    return(out)
  }

# Forecasting formulae

ksi <-
  function(n, m, v, H, nu, k){
  s = n / 252
  t = (n + m) / 252
  fct1 = sin(pi * (H - 0.5))/pi * v ^ (0.5 - H) * (s - v) ^ (0.5 - H)
  vJ = seq(from = n + 1, to = n + m, by = 1)
  fct2  = nu * (vJ / 252) ^ (H - 0.5) * (vJ/252 - s) ^ (H - 0.5) / (vJ/252 - v) * 1/252 *
    exp(- k * (t - vJ/252))
  return(fct1 * sum(fct2))
}

f_expectedval <-
  function(vx, delta, H, k, nu, mu, h){
    vx = na.omit(vx)
  n = length(vx)
  m = n + h
  s = n * delta
  t = m * delta
  tail(vx, 1) * exp(- k * (t - s)) + mu * (1 - exp(- k * (t - s))) - 
    k * mu / nu * sum(sapply(X = c(1 : n)* delta, FUN = ksi, n = n, m = m, nu = nu, H = H, k = k) * delta) +
    k / nu * sum(sapply(X = c(1 : n) * delta, FUN = ksi, n = n, m = m, nu = nu, H = H, k = k) * vx[1 : n] * delta) + 
    1 / nu * sum(sapply(X = c(1 : (n - 1)) * delta, FUN = ksi, n = n, m = m, nu = nu, H = H, k = k) * diff(vx[1 : n])) 
  
}