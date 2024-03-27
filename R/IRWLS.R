irwls_tau <- function(beta_exp,
                      se_exp,
                      beta_out,
                      se_out,
                      beta1,
                      tau1,
                      c,
                      iter_s,
                      tol_tau) {
  for (i in 1:iter_s) {
    tau0 <- tau1
    r <- (beta_out - beta_exp*beta1)/sqrt(beta1^2*se_exp^2 + se_out^2)
    w <- wgt_s(r/tau0, c)
    tau1 <- sqrt(2*mean(w*r^2))
    if (abs(tau1 - tau0)/abs(tau1 + 1e-10) <= tol_tau) {
      break
    }
  }
  if (abs(tau1 - tau0)/abs(tau1 + 1e-10) > tol_tau) {
    warning("S-estimate warning: IRWLS for tau-step did not converge. Consider to increase iter_s.")
  }
  tau1
}
