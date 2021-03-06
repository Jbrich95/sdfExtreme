#' Snowdonia Precipitation Data
#'
#'Data consist of winter (DJF) 12-hour average precipitation rate (mm/day) taken from
#'the UKCP18 climate projection. Data is produced on \eqn{2.2 \times 2.2 km^2} grid-boxes, between
#'the years 1980 and 2000.
#'
#' @docType data
#'
#' @usage data(Snow_Precip)
#'
#' @format A list with 2 elements:
#' \describe{
#' \item{Pr.}{A matrix with 3600 rows and 100 columns of precipitation data.}
#' \item{coords}{A 100 by 2 matrix of lon-lat coordinates, corresponding to the centroid of each grid-box.}
#' }
#'
#' @keywords datasets
#'
#' @references Lowe et al. (2018), Met Office, Exter, UK,
#' (\href{https://www.metoffice.gov.uk/pub/data/weather/uk/ukcp18/science-reports/UKCP18-Overview-report.pdf}{pdf})
#'
#' @examples
#' data(Snow_Precip)
"Snow_Precip"
