irwls_semr_beta <- function(beta_exp,
                            se_exp,
                            beta_out,
                            se_out,
                            mroots,
                            beta1,
                            tau1,
                            c,
                            iter_s,
                            tol_beta) {
  if (mroots == T) {
    beta1_mroots <- NULL
  }
  for (i in 1:iter_s) {
    beta0 <- beta1
    r <- (beta_out - beta_exp*beta0)/sqrt(beta0^2*se_exp^2 + se_out^2)
    t <- (beta_out*se_exp^2*beta0 + beta_exp*se_out^2)/(beta0^2*se_exp^2 + se_out^2)^2
    w <- wgt_m(r/tau1, c)*t
    beta1 <- sum(w*beta_out)/sum(w*beta_exp)
    if (mroots == T) {
      if (!(beta1 %in% beta1_mroots)){
        beta1_mroots <- c(beta1_mroots, beta1)
      }
    }
    if (abs(beta1 - beta0)/abs(beta1 + 1e-10) <= tol_beta & mroots == F) {
      break
    }
  }
  if (abs(beta1 - beta0)/abs(beta1 + 1e-10) > tol_beta & mroots == F) {
    warning("S-estimate warning: IRWLS for beta-step did not converge. Consider to increase iter_s or check potential multiple roots.")
    beta1
  } else if (abs(beta1 - beta0)/abs(beta1 + 1e-10) <= tol_beta & mroots == F) {
    beta1
  } else if (mroots == T) {
    if (length(beta1_mroots) == iter_s) {
      warning("S-estimate warning: Still exists multiple roots. Consider to increase iter_s.")
    }
    beta1_mroots
  }
}

irwls_semr_tau <- function(beta_exp,
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
