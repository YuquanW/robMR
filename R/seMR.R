#' Title
#'
#' @return
#' @export
#'
#' @examples
seMR <- function(beta_exp, se_exp, beta_out, se_out, k_m = 1, k_s = 1) {
  beta_ivw <- sum(beta_out*beta_exp/se_out^2)/sum(beta_exp^2/se_out^2)

}
