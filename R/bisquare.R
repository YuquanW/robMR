#' Title
#'
#' @param x
#' @param k
#'
#' @return
#' @export
#'
#' @examples
rho <- function(x, k) {
  (1-(1-(x/k)^2)^3)*(abs(x)<=k) + 1*(abs(x)>k)
}

psi <- function(x, k) {
  x*(1-(x/k)^2)^2*(abs(x)<=k)
}

wgt_m <- function(x, k) {
  psi(x, k)/x
}

wgt_s <- function(x, k) {
  rho(x, k)/x^2
}
