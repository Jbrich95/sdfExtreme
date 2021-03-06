#' Australia Summer Temperatures outputs
#'
#'Outputs from deformations and model fits for \code{data(Aus_Heat)}.
#'
#'
#' @docType data
#'
#' @usage data(Aus_Heat_Output)
#'
#' @format A list with 5 elements:
#' \describe{
#' \item{likG.IMSP}{Optim output from fitting Inverted Smith model to G-plane sampling locations. See \code{nllIMSPexpSmith}.}
#' \item{likG.MSP}{Optim output from fitting Brown-Resnick model to G-plane sampling locations. See \code{nllMSPexpSmith}.}
#' \item{likD.IMSP}{Optim output from fitting Inverted Smith model to D-plane sampling locations. See \code{nllIMSPexp}.}
#' \item{likD.MSP}{Optim output from fitting Brown-Resnick model to D-plane sampling locations. See \code{nllMSPexp}.}
#' \item{sdf}{D-plane transformation spline parameters. See \code{sdf.heur.EXTR}.}
#' }
#'
#' @keywords output
#'
#' @references Caesar et al. (2006) J. Geophys. Res. 111, D05101,
#' (\href{https://doi.org/10.1029/2005JD006280}{doi})
#'
#' @examples
#' data(Aus_Heat_Output)
"Aus_Heat_Output"
