#' Existence condition of the covariance matrix of the mfBm
#'
#' @description
#' Functions to verify the existence of the covariance matrix of a multivariate
#' fractional Brownian motion given parameters' values.
#' - `cov_existence_2fbm` verifies the simplified condition arising in the
#' two dimensional case, referred to in the paper as coherency constraint;
#' - `cov_existence_mfbm` verifies the general existence condition in terms
#' of the positive semi-definiteness of an auxiliary matrix.
#'
#' @param H Vector of Hurst coefficients.
#' @param rho Symmetric matrix of instantaneous correlation coefficients.
#' @param eta Antisymmetric matrix of asymmetry coefficients.
#' @param out Logical. If TRUE, the eigenvalues of the auxiliary matrix are
#' returned.
#' @param Hi Hurst coefficient of first component in the 2-dimensional case.
#' @param Hj Hurst coefficient of second component in the 2-dimensional case.
#' @param rhoij Instantaneous correlation coefficient in the 2-dimensional case.
#' @param etaij Asymmetry coefficient in the 2-dimensional case.
#'
#' @returns
#' Both functions print a message regarding the admissibility of the parameters.
#' `cov_existence_2fbm` returns the value of the constrains, whereas `cov_existence_mfbm` returns the
#' eigenvalues of the auxiliary matrix.
#' @export
#' @rdname cov_existence
#'
#' @references
#' Amblard, Pierre-Olivier, et al. "Basic properties of the multivariate fractional Brownian motion." arXiv preprint arXiv:1007.0828 (2010).
#'
#' @examples
#' H = c(0.0934, 0.1345, 0.1490, 0.1664)
#' rho = matrix(c(1.0000, 0.1471, 0.3143, 0.0739,
#'                0.1471, 1.0000, 0.0599, 0.1218,
#'                0.3143, 0.0599, 1.0000, 0.0299,
#'                0.0739, 0.1218, 0.0299, 1.0000),
#'              ncol = 4)
#' eta = matrix(c(0, 0.0736, -0.1752, 0.0318,
#'                -0.0736, 0, - 0.1550, - 0.0073,
#'                0.1752, 0.1550, 0, 0.0522,
#'                -0.0318, 0.0073, -0.0522, 0),
#'              ncol = 4)
#'
#' # the parameters are admissible for a 4-dimensional mfBm
#' cov_existence_mfbm(H, rho, eta)
#'
#'# the parameters are not always admissible for 2-dimensional mfBms
#'for(i in 1 : 4){
#'  for(j in 1 : i) {
#'    cov_existence_2fbm(rho[i, j], eta[i, j], H[i], H[j])
#'  }
#'}
cov_existence_mfbm <- function(H, rho, eta, out = FALSE){
  d <- length(H)
  G <- matrix(0, nr = d, nc = d)
  for (j in (1 : d)) {
    for (k in (1 : d)) {
      Hjk <- H[j] + H[k]
      if (Hjk != 1) {
        pre <- rho[j, k] * sin(pi / 2 * Hjk)
        pim <- - eta[j, k] * cos(pi / 2 * Hjk)
      }
      else {
        pre <- rho[j, k]
        pim <- - pi / 2 * eta[j, k]
      }
      G[j, k] <- gamma(H[j] + H[k] + 1) *
        complex(real = pre, imaginary = pim)
    }
  }
  eigenvalues <- eigen(G)$values
  tmp <- sum(eigenvalues < 0)
  if (tmp == 0)
    cat("Admissible parameters\n")
  if (tmp > 0)
    cat("Inadmissible parameters\n")
  if(out == TRUE) return(eigenvalues)
}

#' @export
#' @rdname cov_existence
cov_existence_2fbm <- function(rhoij, etaij, Hi, Hj, out = FALSE) {
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

  if (value < 1)
    cat("Admissible parameters\n")
  if (value >= 1)
    cat("Inadmissible parameters\n")

  return(value)
}

#' Expression of \eqn{\eta_{i,j}} in the causal mfBm
#'
#' @description
#' Formulas for the computation of the asymmetry matrix \eqn{\eta},
#' respectively  coefficient \eqn{\eta_{i,j}},
#' as a function of the other parameters, \eqn{\rho} and \eqn{H}, respectively
#'  \eqn{\rho_{i,j},\ H_i, \text{and}\ H_j}, in the casusal multivariate fractional
#' Brownian motion. The causal case corresponds to the process depending only on
#' past innovations in the moving average representation.
#'
#' @param rho Symmetric matrix of instantaneous correlation coefficients.
#' @param H Vector of Hurst coefficients.
#' @param rhoij Correlation coefficient between components i and j.
#' @param Hi Hurst coefficient of component i.
#' @param Hj Hurst coefficient of component j.
#'
#' @returns Numeric value for `etaij_caus` or matrix for `eta_caus`.
#' @export
#' @rdname eta_caus
#'
#' @references
#' Amblard, Pierre-Olivier, et al. "Basic properties of the multivariate fractional Brownian motion." arXiv preprint arXiv:1007.0828 (2010).
etaij_caus <- function(rhoij, Hi, Hj) {
  H = Hi + Hj
  if(H != 1) {
    - rhoij * tan(pi/2 * H) * tan(pi/2 * H)
  } else if (H == 1) {
    2 * rhoij / (pi * tan(pi * Hi))
  }
}

#' @export
#' @rdname eta_caus
eta_caus <- function(rho, H){
  tmp = combn(x = 1 : length(H), m = 2)
  eta <- rho
  eta[lower.tri(eta)] <-
    - rho[lower.tri(rho)] * tan(pi/2 * (H[tmp[1, ]] + H[tmp[2, ]])) * tan(pi/2 * (H[tmp[1, ]] - H[tmp[2, ]]))
  diag(eta) = 0
  eta[upper.tri(eta)] <- - t(eta)[upper.tri(eta)]
  return(eta)
}

#' Variance function of the fBm
#'
#' @description
#' Variance function of the fractional Brownian motion, based on model parameters:
#' \deqn{\text{Var}(B_t^H).}
#'
#' @param t Time.
#' @param H Hurst coefficient.
#' @param sigma Scale coefficient (square root of variance at \eqn{t=1}).
#'
#'
#' @returns Numeric value.
#' @export
var_fbm <- function(t, H, sigma = 1) {
  sigma ^ 2 * abs(t) ^ (2 * H)
}

#' Cross-covariance function of the mfBm
#'
#' @description
#' Cross-covariance function of the multivariate fractional Brownian motion,
#' based on model parameters:
#' \deqn{\text{Cov}\left(B_{t_i}^{H_i},B_{t_j}^{H_j}\right).}
#'
#'
#' @param ti Time of component i;
#' @param tj Time of component j;
#' @param Hi Hurst coefficient of component i.
#' @param Hj Hurst coefficient of component j.
#' @param rhoij Correlation coefficient (covariance at \eqn{t=1}).
#' @param etaij Asymmetry coefficient.
#' @param sigmai Scale coefficient of component i (square root of variance at \eqn{t=1}).
#' @param sigmaj Scale coefficient of component j (square root of variance at \eqn{t=1}).
#'
#' @returns Numeric value.
#' @export
#'
#' @references
#' Amblard, Pierre-Olivier, et al. "Basic properties of the multivariate fractional Brownian motion." arXiv preprint arXiv:1007.0828 (2010).
cov_ijts_mfbm <- function(ti, tj, Hi, Hj, rhoij, etaij, sigmai = 1, sigmaj = 1) {
  if(Hi + Hj != 1){
    (sigmai * sigmaj)/2 * (
      (rhoij + etaij * sign(ti)) * abs(ti) ^ (Hi + Hj) +
        (rhoij - etaij * sign(tj)) * abs(tj) ^ (Hi + Hj) -
        (rhoij - etaij * sign(tj - ti)) * abs(tj - ti) ^ (Hi + Hj))

  } else if (Hi + Hj == 1){
    (sigmai * sigmaj)/2 *
      (rhoij * (abs(ti) + abs(tj) - abs(ti - tj)) +
         etaij * (tj * log(abs(tj)) - ti * log(abs(ti)) - (tj - ti) * log(abs(tj - ti)))
      )
  }
}
