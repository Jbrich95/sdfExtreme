#' Australia Summer Temperatures
#'
#'Data consist of daily summer (DJF) maximum near-surface air temperatures taken from
#'the HadGHCND global gridded dataset and interpolated
#'to 72 grid point locations covering Australia, for the period 1957-2014.
#'
#' @docType data
#'
#' @usage data(Aus_Heat)
#'
#' @format A list with 2 elements:
#' \describe{
#' \item{Temp.}{A matrix with 5234 rows and 72 columns of temperature data.}
#' \item{coords}{A 72 by 2 matrix of lon-lat coordinates.}
#' }
#'
#' @keywords datasets
#'
#' @references Caesar et al. (2006) J. Geophys. Res. 111, D05101,
#' (\href{https://doi.org/10.1029/2005JD006280}{doi})
#'
#' @examples
#' data(Aus_Heat)
"Aus_Heat"
