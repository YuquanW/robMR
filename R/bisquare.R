#' Title
#'
#' @param x
#' @param c
#'
#' @return
#' @export
#'
#' @examples
rho <- function(x, c) {
  (1-(1-(x/c)^2)^3)*(abs(x/c)<=1) + 1*(abs(x/c)>1)
}

psi <- function(x, c) {
  x*(1-(x/c)^2)^2*(abs(x/c)<=1)
}

wgt_m <- function(x, c) {
  psi(x, c)/(x+1e-32)*(x!=0) + 1*(x==0)
}

wgt_s <- function(x, c) {
  rho(x, c)/(x^2+1e-32)*(x!=0) + 6/c^2*(x==0)
}


rho_biw <- function(x, k) {
  (3*(x/k)-3*(x/k)^2+(x/k)^3)*(x<=k) + 1*(x>k)
}

psi_biw <- function(x, k) {
  (3/k-6/k^2*x+3/k^3*x^2)*(x<=k) + 0*(x>k)
}

ppsi_biw <- function(x, k) {
  (-6/k^2+6/k^3*x)*(x<=k) + 0*(x>k)
}
