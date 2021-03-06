#' Highlands Precipitation outputs
#'
#'Outputs from deformations and model fitting for \code{data(High_Precip)}.
#'
#'
#' @docType data
#'
#' @usage data(High_Precip_Output)
#'
#' @format A list with 5 elements:
#' \describe{
#' \item{likG.IMSP}{Optim output from fitting Inverted Brown-Resnick model to G-plane sampling locations. See \code{nllIMSPexp}.}
#' \item{likG.MSP}{Optim output from fitting Brown-Resnick model to G-plane sampling locations. See \code{nllMSPexp}.}
#' \item{likD.IMSP}{Optim output from fitting Inverted Brown-Resnick model to D-plane sampling locations. See \code{nllIMSPexp}.}
#' \item{likD.MSP}{Optim output from fitting Brown-Resnick model to D-plane sampling locations. See \code{nllMSPexp}.}
#' \item{sdf}{D-plane transformation spline parameters. See \code{sdf.heur.EXTR}.}
#' }
#'
#' @keywords output
#'
#' @references Lowe et al. (2018) Met Office, Exter, UK,
#' (\href{https://www.metoffice.gov.uk/pub/data/weather/uk/ukcp18/science-reports/UKCP18-Overview-report.pdf}{pdf})
#'
#' @examples
#' data(High_Precip_Output)
"High_Precip_Output"
