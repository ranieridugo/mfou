################################
# VARIANCE

#' @export
fVar_fOU =
  function(a, v, H){
    value = v ^ 2/(2 * a ^ (2 * H)) * gamma(1 + 2 * H)
    return(unname(value))
  }

#' @export
dVda =
  function(a, v, H){
    value = - v ^ 2 * H * gamma(1 + 2 * H) * a ^ (- (1 + 2 * H))
    return(value)
  }

#' @export
dVdv = function(a, v, H){
  value = v * a ^ (- 2 * H) * gamma(1 + 2 * H)
  return(value)
}

#' @export
dVdH = function(a, v, H){
  value =
    v ^ 2 * gamma(1 + 2 * H) * a ^ (- 2 * H) *
    (pracma::psi(0, 1 + 2 * H) - log(a))
  return(value)
}

# a1, a2, v1, v2, rho, eta, H1, H2
#' @export
dV1dtheta = function(a, v, H){
  vVals =
    c(dVda(a = a, v = v, H = H), 0, dVdv(a = a, v = v, H = H), 0, 0, 0,
      dVdH(a = a, v = v, H = H), 0)
  return(vVals)
}

# # check
# calculus::gradient(f = fVar_fOU, var = c(a = 0.1, v = 1, H = 0.4))
# dV1dtheta(0.1, 1, 0.4)

#' @export
dV2dtheta = function(a, v, H){
  vVals =
    c(0, dVda(a = a, v = v, H = H), 0, dVdv(a = a, v = v, H = H), 0, 0,
      0, dVdH(a = a, v = v, H = H))

  return(vVals)
}

################################
# AUTOVARIANCE

#' @export
C11 =
  function(a, H, delta) {
    integrand1 = function(u){u^(2 * (H - 1)) * (exp(a * u) - exp(- a * u))}
    integrand2 = function(u){u^(2 * (H - 1)) * exp(- a * u)}

    value = 1 / (2 * a) * (
      integrate(f = integrand1, lower = 0, upper = delta)$value +
        # calculus::integral(f = integrand1, bounds = list(u = c(0, delta)))$value +
        (exp(2 * a * delta) - 1) *
        integrate(f = integrand2, lower = delta, upper = Inf)$value
      # calculus::integral(f = integrand2, bounds = list(u = c(delta, Inf)))$value
    )

    return(value)
  }

# fAutocov_fOU =
#
#   function(a, v, H, delta){
#
#     integrand = function(y) exp(- abs(y)) * abs(delta * a + y) ^ (2 * H)
#     value = v^2/(2*a^(2*H)) *
#       (0.5 *
#          integrate(f = integrand, lower = - Inf, upper = Inf)$value
#        - abs(a * delta)^(2*H)
#       )
#
#     return(unname(value))
#
#   }

# f_Autocov_fOU = function(a, v, H, delta){
#   cov_iik_mfou(a, v, H, delta)
# }

#' @export
dC11da = function(a, H, delta) {
  integrand1 = function(y){y^(2 * H - 2) * (exp(a * y) - exp(- a * y))}
  integrand2 = function(y){y^(2 * H - 1) * (exp(a * y) + exp(- a * y))}
  integrand3 = function(y){y^(2 * H - 2) * exp(- a * y)}
  integrand4 = function(y){y^(2 * H - 1) * exp(- a * y)}

  value =
    - 1 / (2 * a ^ 2) *
    # calculus::integral(f = integrand1, bounds = list(y = c(0, delta)))$value +
    integrate(f = integrand1, lower = 0, upper = delta)$value +
    1 / (2 * a) *
    # calculus::integral(f = integrand2, bounds = list(y = c(0, delta)))$value +
    integrate(f = integrand2, lower = 0, upper = delta)$value +
    (- (exp(2 * a * delta) - 1)/(2 * a ^ 2) + delta / a * exp(2 * a * delta)) *
    # calculus::integral(f = integrand3, bounds = list(y = c(delta, Inf)))$value -
    integrate(f = integrand3, lower = delta, upper = Inf)$value -
    (exp(2 * a * delta) - 1) / (2 * a) *
    integrate(f = integrand4, lower = delta, upper = Inf)$value
  # calculus::integral(f = integrand4, bounds = list(y = c(delta, Inf)))$value

  return(value)
}

#' @export
dC11dH = function(a, H, delta) {

  integrand1 = function(u){log(u) * u ^ (2 * (H - 1)) * exp(- a * u)}
  integrand2 = function(u){log(u) * u ^ (2 * (H - 1)) * (exp(a * u) - exp(- a * u))}

  value = 1 / a * (
    tryCatch(integrate(f = integrand2, lower = 0, upper = delta)$value,
             error = function(c) {
               calculus::integral(f = integrand2, bounds = list(u = c(0, delta)))$value}) +
      (exp(2 * a * delta) - 1) *
      # calculus::integral(f = integrand1, bounds = list(u = c(delta, Inf)))$value
      integrate(f = integrand1, lower = delta, upper = Inf)$value
  )
  return(value)
}

# # check
# for(i in 1 : 10){
# print(calculus::gradient(f = C11, params = list(delta = i), var = c(a = 0.1, H = 0.4)))
# print(c(dC11da(0.1, 0.4, delta = i), dC11dH(0.1, 0.4, delta = i)))
# }

#' @export
dGisdv = function(a, v, H, delta) {
  value =
    exp(- a * delta) * dVdv(a = a, v = v, H = H) +
    2 * v * exp(- a * delta) * H * (2 * H - 1) *
    C11(a = a, H = H, delta = delta)

  return(value)
}

#' @export
dGisdH = function(a, v, H, delta) {
  value =
    exp(- a * delta) * dVdH(a = a, v = v, H = H) +
    v ^ 2 * exp(- a * delta) *
    ((4 * H - 1) * C11(a = a, H = H, delta = delta) +
       H * (2 * H - 1) * dC11dH(a = a, H = H, delta = delta))
  return(value)
}

#' @export
dGisda = function(a, v, H, delta) {
  value =
    - delta * exp(- a * delta) * fVar_fOU(a = a, v = v, H = H) +
    exp(- a * delta) * dVda(a = a, v = v, H = H) -
    delta * v ^ 2 * exp(- a * delta) * H * (2 * H - 1) *
    C11(a = a, H = H, delta = delta) +
    v ^ 2 * exp(- a * delta) * H * (2 * H - 1) *
    dC11da(a = a, H = H, delta = delta)
  return(value)
}

# # check
# for(i in 1 : 10){
#   delta = 1 / 252
#   print(i)
#   print(calculus::gradient(f = cov_iik_mfou,
#                            var = c(ai = 0.1, nui = 1.2, Hi = 0.2), params = list(k = delta)))
#   print(
#     c(dGisda(a = 0.1, v = 1.2, H = 0.2, delta = delta),
#       dGisdv(a = 0.1, v = 1.2, H = 0.2, delta = delta),
#       dGisdH(a = 0.1, v = 1.2, H = 0.2, delta = delta))
#   )
# }


# Derivatives of I21

#' @export
dC21da1 =
  function(a1, a2, H1, H2, delta) {
    H = H1 + H2

    integrand1 = function(u){u ^ (H - 2) * (exp(a2 * u) - exp(- a1 * u))}
    integrand2 = function(u){u ^ (H - 1) * exp(- a1 * u)}
    integrand3 = function(u){u ^ (H - 2) * exp(- a1 * u)}

    dCom = 1 / (a1 + a2)
    value =
      - dCom ^ 2 *
      # calculus::integral(f = integrand1, bounds = list(u = c(0, delta)))$value +
      integrate(f = integrand1, lower = 0, upper = delta)$value +
      dCom *
      # calculus::integral(f = integrand2, bounds = list(u = c(0, delta)))$value +
      integrate(f = integrand2, lower = 0, upper = delta)$value +
      dCom * delta * exp((a1 + a2) * delta) *
      # calculus::integral(f = integrand3, bounds = list(u = c(delta, Inf)))$value -
      integrate(f = integrand3, lower = delta, upper = Inf)$value -
      dCom ^ 2 * (exp((a1 + a2) * delta) - 1) *
      # calculus::integral(f = integrand3, bounds = list(u = c(delta, Inf)))$value -
      integrate(f = integrand3, lower = delta, upper = Inf)$value -
      dCom * (exp((a1 + a2) * delta) - 1) *
      # calculus::integral(f = integrand2, bounds = list(u = c(delta, Inf)))$value
      integrate(f = integrand2, lower = delta, upper = Inf)$value

    return(value)
  }

#' @export
dC21da2 =
  function(a1, a2, H1, H2, delta) {
    H = H1 + H2
    dCom = 1 / (a1 + a2)

    integrand1 = function(u){u ^ (H - 2) * (exp(a2 * u) - exp(- a1 * u))}
    integrand2 = function(u){u ^ (H - 1) * exp(a2 * u)}
    integrand3 = function(u){u ^ (H - 2) * exp(- a1 * u)}

    value =
      - dCom ^ 2 *
      # calculus::integral(f = integrand1, bounds = list(u = c(0, delta)))$value +
      integrate(f = integrand1, lower = 0, upper = delta)$value +
      dCom *
      # calculus::integral(f = integrand2, bounds = list(u = c(0, delta)))$value +
      integrate(f = integrand2, lower = 0, upper = delta)$value +
      (dCom * delta * exp((a1 + a2) * delta) -
         dCom ^ 2 * (exp((a1 + a2) * delta) - 1)) *
      # calculus::integral(f = integrand3, bounds = list(u = c(delta, Inf)))$value
      integrate(f = integrand3, lower = delta, upper = Inf)$value

    return(value)
  }

#' @export
dC21dH =
  function(a1, a2, H1, H2, delta) {

    dCom = 1 / (a1 + a2)
    H = H1 + H2

    integrand1 = function(u){log(u) * u ^ (H - 2) * (exp(a2 * u) - exp(- a1 * u))}
    integrand2 = function(u){log(u) * u ^ (H - 2) * exp(- a1 * u)}

    value =
      dCom *
      tryCatch(integrate(f = integrand1, lower = 0, upper = delta)$value,
               error = function(c) {
                 calculus::integral(integrand1, bounds = list(u = c(0, delta)))$value}) +
      dCom * (exp((a1 + a2) * delta) - 1) *
      # calculus::integral(f = integrand2, bounds = list(u = c(delta, Inf)))$value
      integrate(f = integrand2, lower = delta, upper = Inf)$value

    return(value)
  }

# # check
# for(i in 1 : 5){
#   print(i)
#   # print(calculus::gradient(f = C21, var = c(a1 = 0.1, a2 = 0.2, H1 = 0.5, H2 = 0.5),
#   #                          params = list(delta = i)))
#   print(
#     c(
#       pracma::grad(I_ji, x0 = 0.1, aj = 0.2, Hi = 0.4, Hj = 0.6, k = i),
#       pracma::grad(I_ji, ai = 0.1, x0 = 0.2, Hi = 0.4, Hj = 0.6, k = i),
#       pracma::grad(I_ji, ai = 0.1, aj = 0.2, x0 = 0.4, Hj = 0.6, k = i))
#   )
#
#   print(
#     c(dC21da1(a1 = 0.1, a2 = 0.2, H1 = 0.4, H2 = 0.6, delta = i),
#       dC21da2(a1 = 0.1, a2 = 0.2, H1 = 0.4, H2 = 0.6, delta = i),
#       dC21dH(a1 = 0.1, a2 = 0.2, H1 = 0.4, H2 = 0.6, delta = i))
#   )
# }

## Derivative of I12

#' @export
dC12da1 =
  function(a1, a2, H1, H2, delta) {
    H = H1 + H2
    dCom = 1 / (a1 + a2)

    integrand1 = function(u){u ^ (H - 2) * (exp(a1 * u) - exp(- a2 * u))}
    integrand2 = function(u){u ^ (H - 1) * exp(a1 * u)}
    integrand3 = function(u){u ^ (H - 2) * exp(- a2 * u)}

    value =
      - dCom ^ 2 *
      # calculus::integral(f = integrand1, bounds = list(u = c(0, delta)))$value +
      integrate(f = integrand1, lower = 0, upper = delta)$value +
      dCom *
      # calculus::integral(f = integrand2, bounds = list(u = c(0, delta)))$value +
      integrate(f = integrand2, lower = 0, upper = delta)$value +
      (dCom * delta * exp((a1 + a2) * delta) -
         dCom ^ 2 * (exp((a1 + a2) * delta) - 1)) *
      # calculus::integral(f = integrand3, bounds = list(u = c(delta, Inf)))$value
      integrate(f  = integrand3, lower = delta, upper = Inf)$value

    return(value)
  }


#' @export
dC12da2 =
  function(a1, a2, H1, H2, delta) {
    H = H1 + H2

    integrand1 = function(u){u ^ (H - 2) * (exp(a1 * u) - exp(- a2 * u))}
    integrand2 = function(u){u ^ (H - 1) * exp(- a2 * u)}
    integrand3 = function(u){u ^ (H - 2) * exp(- a2 * u)}

    dCom = 1 / (a1 + a2)
    value =
      - dCom ^ 2 *
      integrate(f = integrand1, lower = 0, upper = delta)$value +
      # calculus::integral(f = integrand1, bounds = list(u = c(0, delta)))$value +
      dCom *
      integrate(f = integrand2, lower = 0, upper = delta)$value +
      # calculus::integral(f = integrand2, bounds = list(u = c(0, delta)))$value +
      dCom * delta * exp((a1 + a2) * delta) *
      integrate(f = integrand3, lower = delta, upper = Inf)$value -
      # calculus::integral(f = integrand3, bounds = list(u = c(delta, Inf)))$value -
      dCom ^ 2 * (exp((a1 + a2) * delta) - 1) *
      integrate(f = integrand3, lower = delta, upper = Inf)$value -
      # calculus::integral(f = integrand3, bounds = list(u = c(delta, Inf)))$value -
      dCom * (exp((a1 + a2) * delta) - 1) *
      integrate(f = integrand2, lower = delta, upper = Inf)$value
    # calculus::integral(f = integrand2, bounds = list(u = c(delta, Inf)))$value

    return(value)
  }

#' @export
dC12dH =
  function(a1, a2, H1, H2, delta) {

    dCom = 1 / (a1 + a2)
    H = H1 + H2

    integrand1 = function(u){log(u) * u ^ (H - 2) * (exp(a1 * u) - exp(- a2 * u))}
    integrand2 = function(u){log(u) * u ^ (H - 2) * exp(- a2 * u)}

    value =
      dCom *
      # calculus::integral(integrand1, bounds = list(u = c(0, delta)))$value +
      tryCatch( integrate(f = integrand1, lower = 0, upper = delta)$value,
                error = function(c) {
                  calculus::integral(integrand1, bounds = list(u = c(0, delta)))$value}) +

      dCom * (exp((a1 + a2) * delta) - 1) *
      integrate(f = integrand2, lower = delta, upper = Inf)$value
    # calculus::integral(f = integrand2, bounds = list(u = c(delta, Inf)))$value

    return(value)
  }

# # check
# for(i in 1 : 10){
#   print(i)
#   print(calculus::gradient(f = I_ij, var = c(ai = 0.1, aj = 0.2, Hi = 0.3, Hj = 0.4),
#                            params = list(k = i)))
#   print(
#     c(dC12da1(a1 = 0.1, a2 = 0.2, H1 = 0.3, H2 = 0.4, delta = i),
#       dC12da2(a1 = 0.1, a2 = 0.2, H1 = 0.3, H2 = 0.4, delta = i),
#       dC12dH(a1 = 0.1, a2 = 0.2, H1 = 0.3, H2 = 0.4, delta = i))
#   )
# }

################################
# COVARIANCE

#' @export
fCov_bfOU =
  function(a1, a2, v1, v2, H1, H2, rho, eta){

    H = H1 + H2

    if(H != 1){
      value = 1/(2 * (a1 + a2)) *
        gamma(H + 1) * v1 * v2 * ((a1 ^ (1 - H) + a2 ^ (1 - H)) * rho +
                                    (a2 ^ (1 - H) - a1 ^ (1 - H)) * eta)
    } else if (H == 1) {
      value = v1 * v2/(a1 + a2) * (rho + eta * 0.5 * (log(a2) - log(a1)))
    }
    return(unname(value))
  }

#' @export
dG0da1 =
  function(a1, a2, v1, v2, H1, H2, rho, eta) {
    H = H1 + H2
    if(H == 1){
      value =
        - eta * v1 * v2 / (2 * a1 * (a1 + a2)) - v1 * v2 *
        (rho + 0.5 * eta * (- log(a1) + log(a2))) / ((a1 + a2) ^ 2)

    } else {
      dComm1 = gamma(H + 1) * v1 * v2 / (2 * (a1 + a2))
      dComm2 = ((a1 ^ (1 - H) + a2 ^ (1 - H)) * rho +
                  (a2 ^ (1 - H) -  a1 ^ (1 - H)) * eta)
      value =
        dComm1 *
        (
          - 1/(a1 + a2) * dComm2 + (1 - H) * a1 ^ (- H) * (rho - eta)
        )

    }
    return(value)
  }

#' @export
dG0da2 =
  function(a1, a2, v1, v2, H1, H2, rho, eta) {

    H = H1 + H2
    if(H == 1){

      value =
        v1 * v2 *
        (- 1 / (a1 + a2) ^ 2 * (rho + 00.5 * eta * (log(a2) - log(a1))) +
           1 / (a1 + a2) * 0.5 * eta / a2)

    } else {
      dComm1 = gamma(H + 1) * v1 * v2 / (2 * (a1 + a2))
      dComm2 = ((a1 ^ (1 - H) + a2 ^ (1 - H)) * rho +
                  (a2 ^ (1 - H) -  a1 ^ (1 - H)) * eta)

      value =
        dComm1 * (
          - 1 / (a1 + a2) * dComm2 + (1 - H) * a2 ^ (- H) * (rho + eta)
        )
    }
    return(value)
  }

#' @export
dG0dH = function(a1, a2, v1, v2, H1, H2, rho, eta) {

  H = H1 + H2
  if(H == 1){

    value =
      value = 0

  } else {
    dComm1 = gamma(H + 1) * v1 * v2 / (2 * (a1 + a2))
    dComm2 = ((a1 ^ (1 - H) + a2 ^ (1 - H)) * rho +
                (a2 ^ (1 - H) -  a1 ^ (1 - H)) * eta)
    value = dComm1 *
      (pracma::psi(0, H + 1) * dComm2 -
         ((a1 ^ (1 - H) * log(a1) + a2 ^ (1 - H) * log(a2)) * rho +
            (a2 ^ (1 - H) * log(a2) - a1 ^ (1 - H) * log(a1)) * eta))
  }

  return(value)
}

#' @export
dG0cdH1 = # causal
  function(a1, a2, v1, v2, H1, H2, rho){
    H = H1 + H2
    dComm = v1 * v2 * rho /(2 * (a1 + a2))
    t1 = tan(pi/2 * (H1 + H2))
    t2 = tan(pi/2 * (H1 - H2))

    value =
      dComm *
      (
        pracma::psi(0, H + 1) * (
          (a1 ^ (1 - H) + a2 ^ (1 - H)) -
            (a2 ^ (1 - H) - a1 ^ (1 - H)) * t1 * t2) +
          gamma(H + 1) * (
            a2 ^ (1 - H) * log(a2) * (t1 * t2 - 1) -
              a1 ^ (1 - H) * log(a1) * (1 + t1 * t2) -
              (a2 ^ (1 - H) - a1 ^ (1 - H)) * pi/2 * (
                t1 * pracma::sec(pi/2 * (H1 - H2)) ^ 2 + pracma::sec(pi/2 * (H1 + H2)) ^ 2 * t2))
      )

    return(value)
  }

#' @export
dG0cdH2 =
  function(a1, a2, v1, v2, H1, H2, rho){
    H = H1 + H2
    dComm = v1 * v2 * rho /(2 * (a1 + a2))
    t1 = tan(pi/2 * (H1 + H2))
    t2 = tan(pi/2 * (H1 - H2))

    value =
      dComm *
      (
        pracma::psi(0, H + 1) * (
          (a1 ^ (1 - H) + a2 ^ (1 - H)) -
            (a2 ^ (1 - H) - a1 ^ (1 - H)) * t1 * t2) +
          gamma(H + 1) * (
            a2 ^ (1 - H) * log(a2) * (t1 * t2 - 1) -
              a1 ^ (1 - H) * log(a1) * (1 + t1 * t2) +
              (a2 ^ (1 - H) - a1 ^ (1 - H)) * pi/2 * (
                t1 * pracma::sec(pi/2 * (H1 - H2)) ^ 2 - pracma::sec(pi/2 * (H1 + H2)) ^ 2 * t2))
      )

    return(value)
  }

#' @export
dG0cdrho =
  function(a1, a2, v1, v2, H1, H2){

    eta = - tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))
    H = H1 + H2

    value = 1/(2 * (a1 + a2)) *
      gamma(H + 1) * v1 * v2 * ((a1 ^ (1 - H) + a2 ^ (1 - H)) +
                                  (a2 ^ (1 - H) - a1 ^ (1 - H)) * eta)

    return(value)
  }

#' @export
dG0dv1 = function(a1, a2, v1, v2, H1, H2, rho, eta) {

  H = H1 + H2
  if(H == 1){

    value =
      v2  / (a1 + a2) * (rho + 0.5 * eta * (log(a2) - log(a1)))

  } else {

    dComm2 = ((a1 ^ (1 - H) + a2 ^ (1 - H)) * rho +
                (a2 ^ (1 - H) -  a1 ^ (1 - H)) * eta)
    value = gamma(H + 1) * v2 / (2 * (a1 + a2)) * dComm2
  }

  return(value)
}

#' @export
dG0dv2 = function(a1, a2, v1, v2, H1, H2, rho, eta) {

  H = H1 + H2
  if(H == 1){

    value =
      v1 / (a1 + a2) * (rho + 0.5 * eta * (log(a2) - log(a1)))

  } else {

    dComm2 = ((a1 ^ (1 - H) + a2 ^ (1 - H)) * rho +
                (a2 ^ (1 - H) -  a1 ^ (1 - H)) * eta)
    value = gamma(H + 1) * v1 / (2 * (a1 + a2)) * dComm2
  }
  return(value)
}

#' @export
dG0drho = function(a1, a2, v1, v2, H1, H2, rho, eta) {

  H = H1 + H2
  if(H == 1){

    value =
      v1 * v2 / (a1 + a2)

  } else {
    dComm1 = gamma(H + 1) * v1 * v2 / (2 * (a1 + a2))
    value = dComm1 * (a1 ^ (1 - H) + a2 ^ (1 - H))
  }

  return(value)
}

#' @export
dG0deta = function(a1, a2, v1, v2, H1, H2, rho, eta) {

  H = H1 + H2
  if(H == 1){

    value =
      0.5 * v1 * v2 / (a1 + a2) * (log(a2) - log(a1))

  } else {
    dComm1 = gamma(H + 1) * v1 * v2 / (2 * (a1 + a2))
    value = dComm1 * (a2 ^ (1 - H) - a1 ^ (1 - H))
  }

  return(value)

}

# # test
# print(pracma::grad(f = cov_ij0_mfou, x = 1.1, aj = 0.5,
#                    nui = 1, nuj = 3, Hi = 0.4, Hj = 0.6,
#                    rhoij = 0.5, etaij =  0.3))
# print(dG0da1(a1 = 1.1, a2 = 0.5, v1 = 1, v2 = 3, H1 = 0.4, H2 = 0.6,
#               rho = 0.5, eta = 0.3))
#
# print(pracma::grad(f = cov_ij0_mfou, ai = 1.1, x0 = 0.5,
#                    nui = 1, nuj = 3, Hi = 0.4, Hj = 0.6,
#                    rhoij = 0.5, etaij =  0.3))
# print(dG0da2(a1 = 1.1, a2 = 0.5, v1 = 1, v2 = 3, H1 = 0.4, H2 = 0.6,
#               rho = 0.5, eta = 0.3))
#
# print(pracma::grad(f = cov_ij0_mfou, ai = 1.1, aj = 0.5,
#                    nui = 1, nuj = 3, Hi = 0.1, x0 = 0.2,
#                    rhoij = 0.5, etaij =  0.3))
# print(dG0dH(a1 = 1.1, a2 = 0.5, v1 = 1, v2 = 3, H1 = 0.1, H2 = 0.2,
#               rho = 0.5, eta = 0.3))
#
# print(pracma::grad(f = cov_ij0_mfou, ai = 1.1, aj = 0.5,
#                    nui = 1, nuj = 3, Hi = 0.6, Hj = 0.7,
#                    rhoij = 0.5, x0 =  0.3))
# print(dG0deta(a1 = 1.1, a2 = 0.5, v1 = 1, v2 = 3, H1 = 0.6, H2 = 0.7,
#               rho = 0.5, eta = 0.3))
#
# print(pracma::grad(f = cov_ij0_mfou, ai = 1.1, aj = 0.5,
#                    nui = 1, nuj = 3, Hi = 0.4, Hj = 0.6,
#                    x0 = 0.5, etaij =  0.3))
# print(dG0drho(a1 = 1.1, a2 = 0.5, v1 = 1, v2 = 3, H1 = 0.4, H2 = 0.6,
#               rho = 0.5, eta = 0.3))

################################
# CROSS-COVARIANCE


## CCF2

#' @export
dGs21da1 =

  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    if(H == 1){ #  WRONG!
      value =
        exp(- a2 * delta) * dG0da1(a1, a2, v1, v2, H1, H2, rho, eta) +
        0.5 * v1 * v2 * exp(- a2 * delta) * eta *
        dC21da1(a1, a2, H1, H2, delta)

    } else {
      value =
        exp(- a2 * delta) *
        dG0da1(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
               H1 = H1, H2 = H2, rho = rho, eta = eta) +
        v1 * v2 * exp(- a2 * delta) * H * (H - 1) * (rho - eta) / 2 *
        dC21da1(a1 = a1, a2 = a2, H1 = H1, H2 = H2, delta = delta)
    }

    return(value)
  }

# # check
# for(i in 1 : 5){
#   print(i)
#   print(pracma::grad(f = cov_jik_mfou, x0 = 0.1, aj = 0.2,
#                      nui = 1, nuj = 2, Hi = 0.5, Hj = 0.5,
#                      rhoij = 0.5, etaij =  0.3, k = i))
#   print(dGs21da1(a1 = 0.1, a2 = 0.2, v1 = 1, v2 = 2, H1 = 0.5, H2 = 0.5,
#                   rho = 0.5, eta = 0.3, delta = i))
# }

#' @export
dGs21da2 =

  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    if(H == 1){
      value =
        - delta * exp(- a2 * delta) *
        fCov_bfOU(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
                  H1 = H1, H2 = H2, rho = rho, eta = eta) +
        exp(- a2 * delta) * dG0da2(a1, a2, v1, v2, H1, H2, rho, eta) -
        delta * v1 * v2 * exp(- a2 * delta) * eta * 0.5 * I_ji(a1, a2, H1, H2, delta) +
        v1 * v2 * exp(- a2 * delta) * eta * 0.5 * dC21da2(a1, a2, H1, H2, delta)

    } else {
      value =
        - delta * exp(- a2 * delta) *
        fCov_bfOU(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
                  H1 = H1, H2 = H2, rho = rho, eta = eta) +
        exp(- a2 * delta) *
        dG0da2(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
               H1 = H1, H2 = H2, rho = rho, eta = eta) -
        delta * v1 * v2 * exp(- a2 * delta) * H * (H - 1) * (rho - eta)/2 *
        I_ji(a1, a2, H1, H2, delta) +
        v1 * v2 * exp(- a2 * delta) * H * (H - 1) * (rho - eta) / 2 *
        dC21da2(a1 = a1, a2 = a2, H1 = H1, H2 = H2, delta = delta)
    }
    return(value)
  }

# # check
# for(i in 1 : 5){
#   print(i)
#   print(pracma::grad(f = cov_jik_mfou, ai = 0.1, x0 = 0.2,
#                      nui = 1, nuj = 2, Hi = 0.1, Hj = 0.2,
#                      rhoij = 0.5, etaij =  0.3, k = i))
#   print(dGs21da2(a1 = 0.1, a2 = 0.2, v1 = 1, v2 = 2, H1 = 0.1, H2 = 0.2,
#                   rho = 0.5, eta = 0.3, delta = i))
# }

#' @export
dGs21dH =

  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {

    H = H1 + H2
    if(H == 1){
      value =
        v1 * v2 * exp(- a2 * delta) * eta * 0.5 *
        dC21dH(a1, a2, H1, H2, delta)

    } else {
      value =
        exp(- a2 * delta) *
        dG0dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
              H1 = H1, H2 = H2, rho = rho, eta = eta) +
        0.5 * v1 * v2 * exp(- a2 * delta) * (rho - eta) *
        ((2 * H - 1) *
           I_ji(a1, a2, H1, H2, delta) +
           H * (H - 1) *
           dC21dH(a1 = a1, a2 = a2, H1 = H1, H2 = H2, delta = delta))
    }

    return(value)
  }

# # check
# for(i in 1 : 5){
#   print(i)
#   print(pracma::grad(f = cov_jik_mfou, ai = 0.1, aj = 0.2,
#                      nui = 1, nuj = 2, x0 = 0.1, Hj = 0.2,
#                      rhoij = 0.5, etaij =  0.3, k = i))
#   print(dGs21dH(a1 = 0.1, a2 = 0.2, v1 = 1, v2 = 2, H1 = 0.1, H2 = 0.2,
#                   rho = 0.5, eta = 0.3, delta = i))
# }

#' @export
dGs21dv1 =

  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    if(H == 1){
      value =
        exp(- a2 * delta) *
        dG0dv1(a1, a2, v1, v2, H1, H2, rho, eta) +
        v2 * exp(- a2 * delta) * 0.5 * eta * I_ji(a1, a2, H1, H2, delta)

    } else {
      value =
        exp(- a2 * delta) *
        dG0dv1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
               rho = rho, eta = eta) +
        v2 * exp(- a2 * delta) * H * (H - 1) * (rho - eta)/2 *
        I_ji(a1, a2, H1, H2, delta)
    }

    return(value)
  }

#' @export
dGs21dv2 =

  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    if(H == 1){
      value =
        exp(- a2 * delta) *
        dG0dv2(a1, a2, v1, v2, H1, H2, rho, eta) +
        v1 * exp(- a2 * delta) * 0.5 * eta * I_ji(a1, a2, H1, H2, delta)

    } else {
      value =
        exp(- a2 * delta) *
        dG0dv2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
               rho = rho, eta = eta) +
        v1 * exp(- a2 * delta) * H * (H - 1) * (rho - eta)/2 *
        I_ji(a1, a2, H1, H2, delta)
    }

    return(value)
  }


#' @export
dGs21drho =
  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    if(H == 1){
      value =
        exp(- a2 * delta) *
        dG0drho(a1, a2, v1, v2, H1, H2, rho, eta)

    } else {
      value = exp(- a2 * delta) *
        dG0drho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta) +
        v1 * v2 / 2 * exp(- a2 * delta) * H * (H - 1) *
        I_ji(a1, a2, H1, H2, delta)
    }

    return(value)
  }

#' @export
dGs21deta =
  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    if(H == 1){
      value =
        exp(- a2 * delta) * dG0deta(a1, a2, v1, v2, H1, H2, rho, eta) +
        v1 * v2 /2 * exp(- a2 * delta) *
        I_ji(a1, a2, H1, H2, delta)

    } else {
      value = exp(- a2 * delta) *
        dG0deta(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta) -
        v1 * v2 / 2 * exp(- a2 * delta) * H * (H - 1) *
        I_ji(a1, a2, H1, H2, delta)
    }

    return(value)
  }

# # check
# for(i in 1 : 5){
#   print(i)
#   print(pracma::grad(f = cov_jik_mfou, ai = 0.1, aj = 0.2,
#                      nui = 1, nuj = 2, Hi = 0.5, Hj = 0.5,
#                      x0 = 0.5, etaij =  0.3, k = i))
#   print(dGs21drho(a1 = 0.1, a2 = 0.2, v1 = 1, v2 = 2, H1 = 0.5, H2 = 0.5,
#                   rho = 0.5, eta = 0.3, delta = i))
# }


#' @export
dGs21cdrho =  # causal
  function(a1, a2, v1, v2, H1, H2, delta) {

    H = H1 + H2

    value =
      exp(- a2 * delta) *
      dG0cdrho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2) +
      0.5 * v1 * v2 * exp(- a2 * delta) * H * (H - 1) *
      I_ji(a1, a2, H1, H2, delta) *
      (1 + tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2)))

    return(value)
  }

#' @export
dGs21cdH1 =  # causal
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    H = H1 + H2
    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H1 - H2)

    value =
      exp(- a2 * delta) * dG0cdH1(a1, a2, v1, v2, H1, H2, rho) +
      v1 * v2 * exp(- a2 * delta) * rho * 0.5 *
      (
        (2 * H - 1) * I_ji(a1, a2, H1, H2, delta) * (1 + tan(t1) * tan(t2)) +
          H * (H - 1) * (
            dC21dH(a1, a2, H1, H2, delta) * (1 + tan(t1) * tan(t2)) +
              I_ji(a1, a2, H1, H2, delta) * pi/2 * (
                pracma::sec(t1) ^ 2 * tan(t2) + tan(t1) * pracma::sec(t2) ^ 2
              )
          )

      )

    return(value)
  }

#' @export
dGs21cdH2 =  # causal
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    H = H1 + H2
    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H1 - H2)

    value =
      exp(- a2 * delta) * dG0cdH2(a1, a2, v1, v2, H1, H2, rho) +
      v1 * v2 * exp(- a2 * delta) * rho * 0.5 *
      (
        (2 * H - 1) * I_ji(a1, a2, H1, H2, delta) * (1 +  tan(t1) * tan(t2)) +
          H * (H - 1) * (
            dC21dH(a1, a2, H1, H2, delta) * (1 +  tan(t1) * tan(t2)) +
              I_ji(a1, a2, H1, H2, delta) * pi/2 * (
                pracma::sec(t1) ^ 2 * tan(t2) - tan(t1) * pracma::sec(t2) ^ 2
              )
          )

      )

    return(value)
  }

# fCrossCovV2c_bfOU =
#
#   function(a1, a2, v1, v2, H1, H2, rho, delta) {
#
#     H = H1 + H2
#     dCov = fCov_bfOU(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2)))
#     value =
#       exp(- a2 * delta) * dCov + v1 * v2 * exp(- a2 * delta) * H * (H - 1) *
#       (rho - - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))) *
#       0.5 * C21(a1, a2, H1, H2, delta)
#
#     # value =
#     #   exp(- a2 * delta) * fCov_bfOU(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))) +
#     #   v1 * v2 * exp(- a2 * delta) * H * (H - 1) * rho / 2 * (1 + tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))) * C21(a1 = a1, a2 = a2, H1 = H1, H2 = H2, delta = delta)
#     return(unname(value))
#   }
#
# i = 2
# a1 = 0.1; a2 = 0.2; v1 = 1; v2 = 2; H1 = 0.3; H2 = 0.4; rho = 0.5
# eta = - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))
# calculus::gradient(f = fCrossCovV2c_bfOU,
#                    var = c(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho),
#                    params = list(delta = i))
#
# c(
#   dGs21da1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs21da2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs21dv1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs21dv2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs21cdH1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs21cdH2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs21cdrho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, delta = i)
# )
#

#' @export
dGs21acdrho =  # asy causal
  function(a1, a2, v1, v2, H1, H2, delta) {
    return(
      - 0.5 * v1 * v2 * (1 + tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))) * delta ^ (H1 + H2)
    )
  }

dGs21acdH1 =  # asy causal
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H1 - H2)

    return(
      - 0.5 * v1 * v2 * rho * delta ^ (H1 + H2) * (
        pi/2 * (pracma::sec(t1) ^ 2 * tan(t2) + tan(t1) * pracma::sec(t2) ^ 2) + (1 + tan(t1) * tan(t2)) * log(delta))
    )
  }

#' @export
dGs21acdH2 =  # asy causal
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H1 - H2)

    return(
      - 0.5 * v1 * v2 * rho * delta ^ (H1 + H2) * (
        pi/2 * (pracma::sec(t1) ^ 2 * tan(t2) - tan(t1) * pracma::sec(t2) ^ 2) + (1 + tan(t1) * tan(t2)) * log(delta))
    )
  }

# fAsyCrossCovV2_bfOU =
#   function(v1, v2, H1, H2, rho, cov, delta) {
#     eta = - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))
#     value = cov - 0.5 * v1 * v2 * (rho - eta) * delta ^ (H1 + H2)
#     return(value)
#   }
#
# i = 2
# v1 = 1; v2 = 2; H1 = 0.3; H2 = 0.4; rho = 0.5; cov = 0.3
# eta = - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))
# calculus::gradient(f = fAsyCrossCovV2_bfOU,
#                    var = c(v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, cov = cov),
#                    params = list(delta = i))
#
# c(
#   - (rho - eta)/2 * v2 * i ^ (H1 + H2),
#   - (rho - eta)/2 * v1 * i ^ (H1 + H2),
#   dGs21acdH1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs21acdH2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs21acdrho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, delta = i),1
# )

## CCF1

#' @export
dGs12da1 =
  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    value =
      - delta * exp(- a1 * delta) *
      fCov_bfOU(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
                H1 = H1, H2 = H2, rho = rho, eta = eta) +
      exp(- a1 * delta) *
      dG0da1(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
             H1 = H1, H2 = H2, rho = rho, eta = eta) -
      delta * v1 * v2 * exp(- a1 * delta) * H * (H - 1) * 0.5 *
      (rho + eta) * I_ij(a1, a2, H1, H2, delta)+
      v1*v2*exp(-a1 * delta) * H * (H - 1) * 0.5 * (rho + eta) *
      dC12da1(a1 = a1, a2 = a2, H1 = H1, H2 = H2, delta = delta)

    return(value)
  }

# # check
# for(i in 1 : 10){
#   print(i)
#   print(pracma::grad(f = cov_ijk_mfou, x0 = 0.1, aj = 0.1,
#                      nui = 1, nuj = 2, Hi = 0.1, Hj = 0.2,
#                      rhoij = 0.5, etaij = 0.3, k = i))
#   print(dGs12da1(a1 = 0.1, a2 = 0.1, v1 = 1, v2 = 2, H1 = 0.1, H2 = 0.2,
#               rho = 0.5, eta = 0.3, delta = i))
# }

#' @export
dGs12da2 =

  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    value =
      exp(- a1 * delta) *
      dG0da2(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
             H1 = H1, H2 = H2, rho = rho, eta = eta) +
      v1 * v2 * exp(- a1 * delta) * H * (H - 1) * (rho + eta) / 2 *
      dC12da2(a1 = a1, a2 = a2, H1 = H1, H2 = H2, delta = delta)

    return(value)
  }


#' @export
dGs12dH =

  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {

    H = H1 + H2

    value =
      exp(- a1 * delta) *
      dG0dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2,
            H1 = H1, H2 = H2, rho = rho, eta = eta) +

      0.5 * v1 * v2 * exp(- a1 * delta) * (rho + eta) *
      ((2 * H - 1) *
         I_ij(a1, a2, H1, H2, delta) +
         H * (H - 1) *
         dC12dH(a1 = a1, a2 = a2, H1 = H1, H2 = H2, delta = delta))

    return(value)
  }


#' @export
dGs12dv1 =
  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    value =
      exp(- a1 * delta) *
      dG0dv1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
             rho = rho, eta = eta) +
      v2 * exp(- a1 * delta) * H * (H - 1) * (rho + eta)/2 *
      I_ij(a1, a2, H1, H2, delta)

    return(value)
  }

#' @export
dGs12dv2 =
  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    value =
      exp(- a1 * delta) *
      dG0dv2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
             rho = rho, eta = eta) +
      v1 * exp(- a1 * delta) * H * (H - 1) * (rho + eta)/2 *
      I_ij(a1, a2, H1, H2, delta)

    return(value)
  }

# # check
# for(i in 1 : 10){
#   lag = i / 252
#   print(i)
#   print(pracma::grad(f = cov_ijk_mfou, ai = 1, aj = 1.5,
#                      nui = 1, nuj = 2, Hi = 0.1, x0 = 0.2,
#                      rhoij = 0.5, etaij = 0.3, k = lag))
#   print(dGs12dH(a1 = 1, a2 = 1.5, v1 = 1, v2 = 2, H1 = 0.1, H2 = 0.2,
#               rho = 0.5, eta = 0.3, delta = lag))
# }

#' @export
dGs12drho =
  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {
    H = H1 + H2
    value = exp(- a1 * delta) *
      dG0drho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
              rho = rho, eta = eta) +
      0.5 * v1 * v2 * exp(- a1 * delta) * H * (H - 1) *
      I_ij(a1, a2, H1, H2, delta)

    return(value)
  }


#' @export
dGs12deta =
  function(a1, a2, v1, v2, H1, H2, rho, eta, delta) {

    H = H1 + H2
    value = exp(- a1 * delta) *
      dG0deta(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
              rho = rho, eta = eta) +
      v1 * v2 / 2 * exp(- a1 * delta) * H * (H - 1) *
      I_ij(a1, a2, H1, H2, delta)

    return(value)
  }

# # check
# for(i in 8 : 13){
#   lag = i / 252
#   print(i)
#   print(pracma::grad(f = cov_ijk_mfou, ai = 1, aj = 1.5,
#                      nui = 1, nuj = 2, Hi = 0.1, Hj = 0.2,
#                      rhoij = 0.5, x0 = 0.3, k = lag))
#   print(dGs12deta(a1 = 1, a2 = 1.5, v1 = 1, v2 = 2, H1 = 0.1, H2 = 0.2,
#               rho = 0.5, eta = 0.3, delta = lag))
# }

#' @export
dGs12cdrho =  # causal
  function(a1, a2, v1, v2, H1, H2, delta) {

    H = H1 + H2

    value =
      exp(- a1 * delta) *
      dG0cdrho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2) +
      0.5 * v1 * v2 * exp(- a1 * delta) * H * (H - 1) *
      I_ij(a1, a2, H1, H2, delta) *
      (1 + tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H2 - H1)))

    return(value)
  }

#' @export
dGs12cdH1 =
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    H = H1 + H2
    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H2 - H1)

    value =
      exp(- a1 * delta) *
      dG0cdH1(a1, a2, v1, v2, H1, H2, rho) +
      v1 * v2 * exp(- a1 * delta) * rho * 0.5 *
      (
        (2 * H - 1) * I_ij(a1, a2, H1, H2, delta) * (1 + tan(t1) * tan(t2)) +
          H * (H - 1) * (
            dC12dH(a1, a2, H1, H2, delta) * (1 + tan(t1) * tan(t2)) +
              I_ij(a1, a2, H1, H2, delta) * pi/2 * (
                pracma::sec(t1) ^ 2 * tan(t2) - tan(t1) * pracma::sec(t2) ^ 2
              )
          )

      )

    return(value)
  }

#' @export
dGs12cdH2 =  # causal
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    H = H1 + H2
    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H2 - H1)

    value =
      exp(- a1 * delta) *
      dG0cdH2(a1, a2, v1, v2, H1, H2, rho) +
      v1 * v2 * exp(- a1 * delta) * rho * 0.5 *
      (
        (2 * H - 1) * I_ij(a1, a2, H1, H2, delta) * (1 +  tan(t1) * tan(t2)) +
          H * (H - 1) * (
            dC12dH(a1, a2, H1, H2, delta) * (1 +  tan(t1) * tan(t2)) +
              I_ij(a1, a2, H1, H2, delta) * pi/2 * (
                pracma::sec(t1) ^ 2 * tan(t2) + tan(t1) * pracma::sec(t2) ^ 2
              )
          )

      )

    return(value)
  }

# fCrossCovV1c_bfOU =
#   function(a1, a2, v1, v2, H1, H2, rho, delta) {
#
#     H = H1 + H2
#     eta = - rho * tan(pi/2*(H1 + H2)) * tan(pi/2*(H1 - H2))
#     dCov = fCov_bfOU(a1, a2, v1, v2, H1, H2, rho, eta)
#
#     if(H != 1) {
#
#       value = exp(- a1 * delta) * dCov + v1 * v2 * exp(- a1 * delta) *
#         H * (H - 1) * (rho + eta) * 0.5 * C12(a1, a2, H1, H2, delta)
#
#     } else if (H == 1){
#
#       value = exp(- a1 * delta) * dCov + v1 * v2 * exp(- a1 * delta) *
#         eta * 0.5 * C12(a1, a2, H1, H2, delta)
#     }
#     return(unname(value))
#   }
#
# i = 5
# a1 = 0.8; a2 = 0.2; v1 = 1; v2 = 2; H1 = 0.3; H2 = 0.4; rho = 0.5
# calculus::gradient(f = fCrossCovV1c_bfOU,
#                    var = c(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho),
#                    params = list(delta = i))
#
# eta = - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H1 - H2))
# calculus::gradient(f = fCrossCovV1c_bfOU,
#                    var = c(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho),
#                    params = list(delta = i))
#
# c(
#   dGs12da1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs12da2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs12dv1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs12dv2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, eta = eta, delta = i),
#   dGs12cdH1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs12cdH2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs12cdrho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, delta = i)
# )

#' @export
dGs12acdrho =  # asy causal
  function(a1, a2, v1, v2, H1, H2, delta) {
    return(
      - 0.5 * v1 * v2 * (1 + tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H2 - H1))) * delta ^ (H1 + H2)
    )
  }

#' @export
dGs12acdH1 =  # asy causal
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H2 - H1)

    return(
      - 0.5 * v1 * v2 * rho * delta ^ (H1 + H2) * (
        pi/2 * (pracma::sec(t1) ^ 2 * tan(t2) - tan(t1) * pracma::sec(t2) ^ 2) + (1 + tan(t1) * tan(t2)) * log(delta))
    )
  }

#' @export
dGs12acdH2 =  # asy causal
  function(a1, a2, v1, v2, H1, H2, rho, delta) {

    t1 = pi / 2 * (H1 + H2)
    t2 = pi / 2 * (H2 - H1)

    return(
      - 0.5 * v1 * v2 * rho * delta ^ (H1 + H2) * (
        pi/2 * (pracma::sec(t1) ^ 2 * tan(t2) + tan(t1) * pracma::sec(t2) ^ 2) + (1 + tan(t1) * tan(t2)) * log(delta))
    )
  }

# fAsyCrossCovV1_bfOU =
#   function(v1, v2, H1, H2, rho, cov, delta) {
#     eta = rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H2 - H1))
#     value =
#       cov - 0.5 * v1 * v2 * (rho + eta) * delta ^ (H1 + H2)
#     return(value)
#   }
#
# i = 2
# v1 = 1; v2 = 2; H1 = 0.3; H2 = 0.4; rho = 0.5; cov = 0.3
# eta = - rho * tan(pi/2 * (H1 + H2)) * tan(pi/2 * (H2 - H1))
# calculus::gradient(f = fAsyCrossCovV1_bfOU,
#                    var = c(v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, cov = cov),
#                    params = list(delta = i))
#
# c(
#   - (rho - eta)/2 * v2 * i ^ (H1 + H2),
#   - (rho - eta)/2 * v1 * i ^ (H1 + H2),
#   dGs12acdH1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs12acdH2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho, delta = i),
#   dGs12acdrho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, delta = i), 1
# )

#######################################################################
## COMPOSITION

# parameters' order
# a1, a2, v1, v2, rho, eta, H1, H2

#' @export
dG1sdtheta =

  function(a1, v1, H1, delta){

    da1 = dGisda(a = a1, v = v1, H = H1, delta = delta)
    da2 = 0

    dv1 = dGisdv(a = a1, v = v1, H = H1, delta = delta)
    dv2 = 0

    drho = 0
    deta = 0

    dH1 = dGisdH(a = a1, v = v1, H = H1, delta = delta)
    dH2 = 0

    return(c(da1, da2, dv1, dv2, drho, deta, dH1, dH2))
  }

# # check
# for(i in 1 : 10){
#   print(i)
#   print(calculus::gradient(f = fAutocov_fOU,
#                            var = c(a = 0.1, v = 1.2, H = 0.2), params = list(delta = i)))
#   print(dG1sdtheta(0.1, 1.2, 0.2, delta = i))
# }

#' @export
dG2sdtheta =

  function(a2, v2, H2, delta){

    da1 = 0
    da2 = dGisda(a = a2, v = v2, H = H2, delta = delta)

    dv1 = 0
    dv2 = dGisdv(a = a2, v = v2, H = H2, delta = delta)

    drho = 0
    deta = 0

    dH1 = 0
    dH2 = dGisdH(a = a2, v = v2, H = H2, delta = delta)

    return(c(da1, da2, dv1, dv2, drho, deta, dH1, dH2))
  }

# check
# for(i in 1 : 10){
#   print(i)
#   print(calculus::gradient(f = fAutocov_fOU,
#                            var = c(a = 0.1, v = 1.2, H = 0.2), params = list(delta = i)))
#   print(dG2sdtheta(0.1, 1.2, 0.2, delta = i))
# }

#' @export
dG12sdtheta = function(a1, a2, v1, v2, H1, H2, rho, eta, delta){

  da1 = dGs12da1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)
  da2 = dGs12da2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)

  dv1 = dGs12dv1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)
  dv2 = dGs12dv2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)

  drho = dGs12drho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                   rho = rho, eta = eta, delta = delta)
  deta = dGs12deta(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                   rho = rho, eta = eta, delta = delta)

  dH1 = dGs12dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta, delta = delta)
  dH2 = dGs12dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta, delta = delta)

  return(c(da1, da2, dv1, dv2, drho, deta, dH1, dH2))
}

# check
# dG12sdtheta(a1 = 0.1, a2 = 0.2, v1 = 1, v2 = 2, H1 = 0.3, H2 = 0.4, rho = 0.5, eta = 0.3, delta = 10)
# calculus::gradient(f = fCrossCovV1_bfOU,
#                    var = c("a1" = 0.1, "a2" = 0.2, "v1" = 1, "v2" = 2, "rho" = 0.5, "eta" = 0.3,
#                            "H1" = 0.3, "H2" = 0.4),
#                    params = list(delta = 10))


#' @export
dG21sdtheta = function(a1, a2, v1, v2, H1, H2, rho, eta, delta){

  da1 = dGs21da1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)
  da2 = dGs21da2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)

  dv1 = dGs21dv1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)
  dv2 = dGs21dv2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta, delta = delta)

  drho = dGs21drho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                   rho = rho, eta = eta, delta = delta)
  deta = dGs21deta(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                   rho = rho, eta = eta, delta = delta)

  dH1 = dGs21dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta, delta = delta)
  dH2 = dGs21dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta, delta = delta)

  return(c(da1, da2, dv1, dv2, drho, deta, dH1, dH2))
}

# check
# dG21sdtheta(a1 = 0.1, a2 = 0.2, v1 = 1, v2 = 2, H1 = 0.3, H2 = 0.4, rho = 0.5, eta = 0.3, delta = 2)
# calculus::gradient(f = fCrossCovV2_bfOU,
#                    var = c("a1" = 0.1, "a2" = 0.2, "v1" = 1, "v2" = 2, "rho" = 0.5, "eta" = 0.3,
#                            "H1" = 0.3, "H2" = 0.4),
#                    params = list(delta = 2))


#' @export
dG120dtheta =

  function(a1, a2, v1, v2, H1, H2, rho, eta){

    da1 = dG0da1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta)
    da2 = dG0da2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta)

    dv1 = dG0dv1(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta)
    dv2 = dG0dv2(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                 rho = rho, eta = eta)

    drho = dG0drho(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                   rho = rho, eta = eta)
    deta = dG0deta(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                   rho = rho, eta = eta)

    dH1 = dG0dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta)
    dH2 = dG0dH(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2,
                rho = rho, eta = eta)

    return(c(da1, da2, dv1, dv2, drho, deta, dH1, dH2))
  }

# check
# calculus::gradient(f = fCov_bfOU, var = c(a1 = 0.1, a2 = 0.2,
#                                           v1 = 1, v2 = 2, rho = 0.5, eta = 0.1,
#                                           H1 = 0.3, H2 = 0.4))
# dG120dtheta(a1 = 0.1, a2 = 0.2,
#             v1 = 1, v2 = 2, H1 = 0.3, H2 = 0.4, rho = 0.5, eta = 0.1)
#

#' @export
dGdtheta =

  function(a1, a2, v1, v2, H1, H2, rho, eta, iL, delta = 1){

    mOut =
      rbind(
        # gammas21 vector (opposite order)
        t(sapply(X = c(iL : 1) * delta, FUN = dG21sdtheta,
                 a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho,
                 eta = eta)),
        # gamma012
        dG120dtheta(a1 = a1, a2 = a2, v1 = v1, v2 = v2, H1 = H1, H2 = H2, rho = rho,
                    eta = eta),

        # gammas12 vector missing
        # UNSURE: exchanging parameters' index, and derivatives' order
        t(sapply(X = c(1 : iL) * delta, FUN = dG12sdtheta, zero = zero,
                 a1 = a2, a2 = a1, v1 = v2, v2 = v1, H1 = H2, H2 = H1, rho = rho,
                 eta = - eta)),
        # gamma011
        dV1dtheta(a = a1, v = v1, H = H1),
        # gammas11
        t(sapply(X = c(1 : iL) * delta, FUN = dG1sdtheta, zero = zero,
                 a1 = a1, v1 = v1, H1 = H1)),
        # gamma022
        dV2dtheta(a = a2, v = v2, H = H2),
        # gammas22
        t(sapply(X = c(1 : iL) * delta, FUN = dG2sdtheta, zero = zero,
                 a2 = a2, v2 = v2, H2 = H2))
      )

    return(mOut)
  }

# Asymptotic formulae

#' @export
dG1asdtheta =

  function(v, H, dDelta){

    dv110 = 1
    dv220 = 0
    dv120 = 0
    dnu1 = - v * dDelta ^ (2 * H)
    dnu2 = 0
    dH1 = - v ^ 2 * dDelta ^ (2 * H) * log(dDelta)
    dH2 = 0
    drho = 0
    deta = 0

    return(c(dv110, dv220, dv120, dnu1, dnu2, dH1, dH2, drho, deta))
  }

#' @export
dG2asdtheta =

  function(v, H, dDelta){

    dv110 = 0
    dv220 = 1
    dv120 = 0
    dnu1 = 0
    dnu2 = - v * dDelta ^ (2 * H)
    dH1 = 0
    dH2 = - v ^ 2 * dDelta ^ (2 * H) * log(dDelta)
    drho = 0
    deta = 0

    return(c(dv110, dv220, dv120, dnu1, dnu2, dH1, dH2, drho, deta))
  }



#' @export
dG12asdtheta =

  function(v1, v2, H1, H2, rho, eta, dDelta){

    dv110 = 0
    dv220 = 0
    dv120 = 1
    dnu1 = - (rho + eta)/2 * v2 * dDelta ^ (H1 + H2)
    dnu2 = - (rho + eta)/2 * v1 * dDelta ^ (H1 + H2)
    dH1 = - (rho + eta)/2 * v1 * v2 * dDelta ^ (H1 + H2) * log(dDelta)
    dH2 = - (rho + eta)/2 * v1 * v2 * dDelta ^ (H1 + H2) * log(dDelta)
    drho = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)
    deta = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)

    return(c(dv110, dv220, dv120, dnu1, dnu2, dH1, dH2, drho, deta))
  }

#' @export
dG21asdtheta =

  function(v1, v2, H1, H2, rho, eta, dDelta){

    dv110 = 0
    dv220 = 0
    dv120 = 1
    dnu1 = - (rho - eta)/2 * v2 * dDelta ^ (H1 + H2)
    dnu2 = - (rho - eta)/2 * v1 * dDelta ^ (H1 + H2)
    dH1 = - (rho - eta)/2 * v1 * v2 * dDelta ^ (H1 + H2) * log(dDelta)
    dH2 = - (rho - eta)/2 * v1 * v2 * dDelta ^ (H1 + H2) * log(dDelta)
    drho = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)
    deta = 0.5 * v1 * v2 * dDelta ^ (H1 + H2)

    return(c(dv110, dv220, dv120, dnu1, dnu2, dH1, dH2, drho, deta))
  }

# reduced form: H given

#' @export
dG1asdtheta_red =

  function(v, H, dDelta){

    dv110 = 1
    dv220 = 0
    dv120 = 0
    dnu1 = - v * dDelta ^ (2 * H)
    dnu2 = 0
    drho = 0
    deta = 0

    return(c(dv110, dv220, dv120, dnu1, dnu2, drho, deta))
  }

#' @export
dG2asdtheta_red =

  function(v, H, dDelta){

    dv110 = 0
    dv220 = 1
    dv120 = 0
    dnu1 = 0
    dnu2 = - v * dDelta ^ (2 * H)
    drho = 0
    deta = 0

    return(c(dv110, dv220, dv120, dnu1, dnu2, drho, deta))
  }



#' @export
dG12asdtheta_red =

  function(v1, v2, H1, H2, rho, eta, dDelta){

    dv110 = 0
    dv220 = 0
    dv120 = 1
    dnu1 = - (rho + eta)/2 * v2 * dDelta ^ (H1 + H2)
    dnu2 = - (rho + eta)/2 * v1 * dDelta ^ (H1 + H2)
    drho = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)
    deta = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)

    return(c(dv110, dv220, dv120, dnu1, dnu2, drho, deta))
  }

#' @export
dG21asdtheta_red =

  function(v1, v2, H1, H2, rho, eta, dDelta,
           plusinf = Inf, tol = 0.2){

    dv110 = 0
    dv220 = 0
    dv120 = 1
    dnu1 = - (rho - eta)/2 * v2 * dDelta ^ (H1 + H2)
    dnu2 = - (rho - eta)/2 * v1 * dDelta ^ (H1 + H2)
    drho = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)
    deta = 0.5 * v1 * v2 * dDelta ^ (H1 + H2)

    return(c(dv110, dv220, dv120, dnu1, dnu2, drho, deta))
  }

# only rho and eta

#' @export
dG12asdtheta_biv =

  function(v1, v2, H1, H2, rho, eta, dDelta){

    drho = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)
    deta = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)

    return(c(drho, deta))
  }

#' @export
dG21asdtheta_biv =

  function(v1, v2, H1, H2, rho, eta, dDelta,
           plusinf = Inf, tol = 0.2){

    drho = - 0.5 * v1 * v2 * dDelta ^ (H1 + H2)
    deta = 0.5 * v1 * v2 * dDelta ^ (H1 + H2)

    return(c(drho, deta))
  }
