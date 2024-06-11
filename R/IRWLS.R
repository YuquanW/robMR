irwls_beta1 <- function(beta_exp,
                       se_exp,
                       beta_out,
                       se_out,
                       beta0_new,
                       beta1_new,
                       sigma1_new,
                       c,
                       iter_s,
                       tol_s) {
  for (i in 1:iter_s) {
    beta1_old <- beta1_new
    d <- (beta_out - beta_exp*beta1_old)^2/(beta1_old^2*se_exp^2 + se_out^2)
    s <- (beta_out*se_exp^2*beta1_old + beta_exp*se_out^2)/(beta1_old^2*se_exp^2 + se_out^2)^2
    w <- psi_biw(d/sigma1, c)*s
    beta1_new <- sum(w*beta_out)/sum(w*beta_exp)
    if (abs(beta1_new - beta1_old)/abs(beta1_new + 1e-10) <= tol_s) {
      break
    }
  }
  if (abs(beta1_new - beta1_old)/abs(beta1_new + 1e-10) > tol_s) {
    warning("S-estimate warning: IRWLS for beta-step did not converge. Consider to increase iter_s or check potential multiple roots.")
  }
  beta1_new
}

irwls_sigma <- function(beta_exp,
                        se_exp,
                        beta_out,
                        se_out,
                        beta1,
                        sigma1,
                        c,
                        iter_s,
                        tol_s) {
  for (i in 1:iter_s) {
    sigma0 <- sigma1
    d <- (beta_out - beta_exp*beta1)^2/(beta1^2*se_exp^2 + se_out^2)
    w <- rho_biw(d/sigma0, c)*2/p
    sigma1 <- sum(w*sigma0)
    if (abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) <= tol_s) {
      break
    }
  }
  if (abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) > tol_s) {
    warning("S-estimate warning: IRWLS for sigma-step did not converge. Consider to increase iter_s.")
  }
  sigma1
}
