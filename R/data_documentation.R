#' Log-realized volatility
#'
#' A dataset containing the logarithm of a number of realized volatilities
#' (annualized and tranformed to percentages). This dataset is used in Bibinger,
#' Yu, and Zhang (2025). It is not the same as in the paper Dugo, Giorgio, and
#' Pigato (2025) because of not having permission to disseminate it.
#'
#' @format A data frame with 7309 rows and 8 columns:
#' \describe{
#'   \item{date}{Date of observation.}
#'   \item{AAPL}{log-realized vol for Apple.}
#'   \item{AMGN}{log-realized vol for Amgen.}
#'   \item{ALD}{log-realized vol for Allied Signal.}
#'   \item{AXP}{log-realized vol for American Express.}
#'   \item{BA}{log-realized vol for Boeing.}
#'   \item{CAT}{log-realized vol for Caterpillar.}
#' }
#' @source Data obtained from Risk Lab (https://dachxiu.chicagobooth.edu/#risklab).
"lrv"

#' BOLD signal
#'
#' Per-second observations of a healthy individual's 7T resting-state fMRI, drawn from the Human Connectome Project,
#' based on Tian’s subcortical parcellations. RH and LH stand for right and left hemisphere,
#' respectively.
#'
#' @format A data frame with 3600 rows and 16 columns:
#' \describe{
#'   \item{HIP}{Hippocampus.}
#'   \item{AMY}{Amygdala.}
#'   \item{pTHA}{Posterio Thalamus.}
#'   \item{aTHA}{Anterior Thalamus.}
#'   \item{NAc}{Nucleus Accumbens.}
#'   \item{GP}{Globus Pallidus.}
#'   \item{PUT}{Putamen.}
#'   \item{CAU}{Caudate.}
#'   \item{HIP}{Hippocampus.}
#' }
#' @source Data obtained from the Padova Neuroscience Center.
"bold"
