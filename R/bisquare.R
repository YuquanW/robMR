#' Title
#'
#' @param x
#' @param c
#'
#' @return
#' @export
#'
#' @examples
rho_biw <- function(x, k) {
  (3*(x/k)-3*(x/k)^2+(x/k)^3)*(x<=k) + 1*(x>k)
}

psi_biw <- function(x, k) {
  (3/k-6/k^2*x+3/k^3*x^2)*(x<=k) + 0*(x>k)
}

ppsi_biw <- function(x, k) {
  (-6/k^2+6/k^3*x)*(x<=k) + 0*(x>k)
}

rho <- function(x, k) {
  (1-(1-(x/k)^2)^3)*(abs(x)<=k) + 1*(abs(x)>k)
}

psi <- function(x, k) {
  (x*(1-(x/k)^2)^2)*(abs(x)<=k) + 0*(abs(x)>k)
}

ppsi <- function(x, k) {
  ((1-(x/k)^2)^2-4*(x/k)^2*(1-(x/k)^2))*(abs(x) <= k) + 0*(abs(x)>k)
}

wgt_m <- function(x, c) {
  psi(x, c)/(x+1e-32)*(x!=0) + 1*(x==0)
}
