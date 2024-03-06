#' Title
#'
#' @return
#' @export
#'
#' @examples
seMR <- function(beta_exp,
                 se_exp,
                 beta_out,
                 se_out,
                 mroots = F,
                 c = 1.547,
                 iter_s = 20,
                 tol_s = 1e-5,
                 tol_beta = 1e-5,
                 tol_tau= 1e-3) {
  beta_init <- beta_out/beta_exp
  tau_init <- rep(0, length(beta_init))
  for (i in 1:length(beta_init)) {
    r <- (beta_out - beta_exp*beta_init[i])/sqrt(beta_init[i]^2*se_exp^2 + se_out^2)
    tau <- sqrt(mean(r^2))
    tau_init[i] <- irwls_semr_tau(beta_exp,
                                  se_exp,
                                  beta_out,
                                  se_out,
                                  beta_init[i],
                                  tau,
                                  c,
                                  iter_s,
                                  tol_tau)
  }
  beta1 <- beta_init[which.min(tau_init)]
  tau1 <- min(tau_init)
  for (i in 1:iter_s) {
    beta0 <- beta1
    tau0 <- tau1
    beta1 <- irwls_semr_beta(beta_exp,
                             se_exp,
                             beta_out,
                             se_out,
                             mroots = F,
                             beta0,
                             tau0,
                             c,
                             iter_s,
                             tol_beta)
    tau1 <- irwls_semr_tau(beta_exp,
                           se_exp,
                           beta_out,
                           se_out,
                           beta0,
                           tau0,
                           c,
                           iter_s,
                           tol_tau)
    if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
        abs(tau1 - tau0)/abs(tau1 + 1e-10) <=
        tol_s) {
      break
    }
  }
  if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
      abs(tau1 - tau0)/abs(tau1 + 1e-10) >
      tol_s) {
    warning("S-estimate warning: IRWLS did not converge. Consider to increase iter_s.")
  }
  if (mroots == T) {
    beta1_mroots <- irwls_semr_beta(beta_exp,
                                    se_exp,
                                    beta_out,
                                    se_out,
                                    mroots,
                                    beta1,
                                    tau1,
                                    c,
                                    iter_s,
                                    tol_beta)
    tau1_mroots <- rep(0, length(beta1))
    for (j in 1:length(beta1_mroots)) {
      tau1_mroots[j] <- irwls_semr_tau(beta_exp,
                                       se_exp,
                                       beta_out,
                                       se_out,
                                       beta1_mroots[j],
                                       tau1,
                                       c,
                                       iter_s,
                                       tol_tau)
    }
    beta1 <- beta1_mroots[which.min(tau1_mroots)]
    tau1 <- min(tau1_mroots)
    return(list(beta_s = beta1,
                tau_s = tau1,
                r_hat = (beta_out - beta_exp*beta1)/sqrt(beta1^2*se_exp^2 + se_out^2)))
  } else {
  return(list(beta_s = beta1,
              tau_s = tau1,
              r_hat = (beta_out - beta_exp*beta1)/sqrt(beta1^2*se_exp^2 + se_out^2)))
  }
}
