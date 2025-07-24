#' Normalized forecast error variance decomposition
#'
#' @param H Vector of Hurst coefficients.
#' @param rho Matrix of instantaneous correlation coefficients.
#'
#' @returns A numeric matrix
#' @export
#'
#' @references
#' Dugo, Ranieri, Giacomo Giorgio, and Paolo Pigato. "Multivariate Rough Volatility." arXiv preprint arXiv:2412.14353 (2024).
#'
#' Pesaran, H. Hashem, and Yongcheol Shin. "Generalized impulse response analysis in linear multivariate models." Economics letters 58.1 (1998): 17-29.
nfevd <- function(H, rho){
  if(is.null(names(H))) names(H) <- paste0("y", 1 : length(H))
  if(is.null(colnames(rho))) colnames(rho) <- names(H)
  if(is.null(rownames(rho)) | all(rownames(rho) == as.character(1 : nrow(rho)))) rownames(rho) <- colnames(rho)
  if(is.list(H)) H <- unlist(H)

  fGij =
    function(i, j){
      sqrt(beta(H[i] + 0.5, H[i] + 0.5) * beta(H[j] + 0.5, H[j] + 0.5)) /
        (
          beta(H[i] + 0.5, H[j] + 0.5) * sqrt(sin(pi * H[i]) * sin(pi * H[j]))
        ) *
        sin(pi * (H[i] + H[j])) /
        (cos(pi * H[i]) + cos(pi * H[j])) * rho[i, j]
    }


  d = length(H)
  mpsi = matrix(NA, nrow = d, ncol = d)
  dimnames(mpsi) <- dimnames(rho)
  for(i in 1 : d){
    for(j in 1 : d){
      mpsi[i, j] =
        (fGij(i, j) ^ 2 / sqrt(fGij(j, j))) /
        sum(
          sapply(X = 1 : d,
                 FUN = function(x) (fGij(i = i, j = x) ^ 2 / sqrt(fGij(x, x))))
        )
    }
  }
  return(mpsi)
}

#' Spillover indices
#'
#' @description
#' This function computes spillover indices in the
#' multivariate fractional Ornstein-Uhlenbeck process and multivariate fractional
#' Brownian motion models of size \eqn{d} starting
#' from (causal) estimates for the parameters \eqn{H} and \eqn{\rho}.
#'
#' @param H Vector of Hurst coefficients.
#' @param rho Matrix of diffusion correlation coefficients.
#' @param plot Logical. If TRUE, two plots showing directional and net
#' pairwise spillovers are printed
#'
#' @references
#' Dugo, Ranieri, Giacomo Giorgio, and Paolo Pigato. "Multivariate Rough Volatility." arXiv preprint arXiv:2412.14353 (2024).
#'
#' Diebold, Francis X., and Kamil Yilmaz. "Better to give than to receive: Predictive directional measurement of volatility spillovers." International Journal of forecasting 28.1 (2012): 57-66.
#'
#' @examples
#' my = rv[, c("FTSE", "N225", "SPX", "SSEC")]
#'
#' # mfOU
#' est = est_mde_mfou(my, type = "cau")
#' spillovers(est$univ[1, ], est$rho, plot = TRUE)
#'
#' @returns List of numeric objects.
#' @export
#'
#' @examples
#' # mfBm
#' H = apply(my, 2, est_H_fbm)
#' rho = est_rho_eta_mfbm(my)$rho
#' spillovers(H, rho, plot = TRUE)
spillovers <- function(H, rho, plot = FALSE){
  mpsi = nfevd(H, rho)
  tot = sum(mpsi[upper.tri(mpsi)] + mpsi[lower.tri(mpsi)]) / sum(mpsi)
  d = sum(mpsi)
  diag(mpsi) = 0
  received = rowSums(mpsi) / d * 100 # received by others
  transmitted = colSums(mpsi) / d * 100 # transmitted to others
  net = (colSums(mpsi) - rowSums(mpsi)) / d * 100 # net
  mpair = (t(mpsi) - mpsi) / d * 100 # net pairwise transmitted

  if(plot == TRUE){

    # Plot directional
    plot_direct <-
      data.frame(received = received, transmitted = transmitted, net = net) |>
      tibble::rownames_to_column(var = "index") |>
      arrange(match(index, names(H))) |>
      dplyr::select(- index) |>
      dplyr::mutate(received = - received) |>
      (\(df) {
        df$name <- rownames(df)
        df[, c("name", setdiff(names(df), "name"))]
      })()
      tidyr::pivot_longer(!name, names_to = "spill", values_to = "val") |>
      dplyr::group_by(spill) |>
      dplyr::mutate(
        spill = factor(spill, levels = c("received", "transmitted", "net")),
        max_abs = max(abs(val)),
        scaled_val = val / max_abs
      ) |>
      dplyr::mutate(val = dplyr::if_else(spill == "received", - val, val)) |>
      dplyr::ungroup() |>
      ggplot2::ggplot(ggplot2::aes(x = name, y = val, fill = scaled_val)) +
      ggplot2::geom_col() +
      ggplot2::facet_grid(rows = dplyr::vars(spill), switch = "y", scales = "free") +
      ggplot2::scale_fill_gradient2(
        low = "red4", mid = "white", high = "blue4",
        midpoint = 0, limits = c(-1, 1),
        name = "Scaled Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.text.y.left = ggplot2::element_text(angle = 90, size = 14, vjust = 0.5),
        strip.placement = "outside",
        strip.background = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 12 * 1.1),
        axis.title.x = ggplot2::element_blank(),
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 12)) +
        ggtitle("Plot 1 of 2: Directional spillovers")
    print(plot_direct)

    # Plot pairwise
    plot_pair <-
      mpair |>
      reshape2::melt() |>
      ggplot2::ggplot(ggplot2::aes(Var2, Var1, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue4", mid = "white", high = "red4",
                                    midpoint = 0, name = "Net") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "To", y = "From") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12 * 1.1, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 12 * 1.1, vjust = 0.5),
        axis.title.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 12),
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 12)
      ) +
      ggplot2::ggtitle("Plot 2 of 2: Net pairwise spillovers")
    print(plot_pair)
  }

  return(
    list(
      total = tot,
      received = received,
      transmitted = transmitted,
      net = net,
      pairwise = mpair
    )
  )
}

