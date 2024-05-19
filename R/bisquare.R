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
