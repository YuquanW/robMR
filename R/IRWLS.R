# irwls_beta <- function(beta_exp,
#                        se_exp,
#                        beta_out,
#                        se_out,
#                        beta1,
#                        sigma1,
#                        c,
#                        iter_s,
#                        tol_s) {
#   for (i in 1:iter_s) {
#     beta0 <- beta1
#     d <- (beta_out - beta_exp*beta0)^2/(beta0^2*se_exp^2 + se_out^2)
#     s <- (beta_out*se_exp^2*beta0 + beta_exp*se_out^2)/(beta0^2*se_exp^2 + se_out^2)^2
#     w <- psi_biw(d/sigma1, c)*s
#     beta1 <- sum(w*beta_out)/sum(w*beta_exp)
#     if (abs(beta1 - beta0)/abs(beta1 + 1e-10) <= tol_s) {
#       break
#     }
#   }
#   if (abs(beta1 - beta0)/abs(beta1 + 1e-10) > tol_s) {
#     warning("S-estimate warning: IRWLS for beta-step did not converge. Consider to increase iter_s or check potential multiple roots.")
#   }
#   beta1
# }

irwls_beta_init <- function(beta_exp,
                            se_exp,
                            beta_out,
                            se_out,
                            beta1,
                            iter_s,
                            tol_s) {
  for (i in 1:iter_s) {
    beta0 <- beta1
    s <- (beta_out*beta0*se_exp^2+beta_exp*se_out^2)/(se_exp*se_out*sqrt(beta0^2*se_exp^2+se_out^2))
    r <- (beta_out-beta0*beta_exp)/sqrt(beta0^2*se_exp^2+se_out^2)
    beta1 <- sum((sign(abs(s)*r)/(abs(s)*r+1e-32)*(r!=0)+1*(r==0))*s*se_exp*se_out/(beta0^2*se_exp^2+se_out^2)*beta_out/sqrt(beta0^2*se_exp^2+se_out^2))/
      sum((sign(abs(s)*r)/(abs(s)*r+1e-32)*(r!=0)+1*(r==0))*s*se_exp*se_out/(beta0^2*se_exp^2+se_out^2)*beta_exp/sqrt(beta0^2*se_exp^2+se_out^2))
    if (abs(beta1 - beta0)/abs(beta1 + 1e-10) <= tol_s) {
      break
    }
  }
  if (abs(beta1 - beta0)/abs(beta1 + 1e-10) > tol_s) {
    warning("S-estimate warning: IRWLS for beta-step did not converge. Consider to increase iter_s or check potential multiple roots.")
  }
  beta1
}

irwls_beta <- function(beta_exp,
                       se_exp,
                       beta_out,
                       se_out,
                       beta1,
                       sigma1,
                       c,
                       iter_s,
                       tol_s) {
  for (i in 1:iter_s) {
    beta0 <- beta1
    s <- (beta_out*beta0*se_exp^2+beta_exp*se_out^2)/(se_exp*se_out*sqrt(beta0^2*se_exp^2+se_out^2))
    r <- (beta_out-beta0*beta_exp)/sqrt(beta0^2*se_exp^2+se_out^2)
    beta1 <- sum(wgt_m(abs(s)*r/sigma1, c)*s*se_exp*se_out/(beta0^2*se_exp^2+se_out^2)*beta_out/sqrt(beta0^2*se_exp^2+se_out^2))/
      sum(wgt_m(abs(s)*r/sigma1, c)*s*se_exp*se_out/(beta0^2*se_exp^2+se_out^2)*beta_exp/sqrt(beta0^2*se_exp^2+se_out^2))
    if (abs(beta1 - beta0)/abs(beta1 + 1e-10) <= tol_s) {
      break
    }
  }
  if (abs(beta1 - beta0)/abs(beta1 + 1e-10) > tol_s) {
    warning("S-estimate warning: IRWLS for beta-step did not converge. Consider to increase iter_s or check potential multiple roots.")
  }
  beta1
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
  p <- length(beta_exp)
  for (i in 1:iter_s) {
    sigma0 <- sigma1
    s <- (beta_out*beta1*se_exp^2+beta_exp*se_out^2)/(se_exp*se_out*sqrt(beta1^2*se_exp^2+se_out^2))
    r <- (beta_out-beta1*beta_exp)/sqrt(beta1^2*se_exp^2+se_out^2)
    sigma1 <- sum(rho(r/sigma0, c)*sigma0)*2/p
    if (abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) <= tol_s) {
      break
    }
  }
  if (abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) > tol_s) {
    warning("S-estimate warning: IRWLS for sigma-step did not converge. Consider to increase iter_s.")
  }
  sigma1
}

# irwls_sigma <- function(beta_exp,
#                         se_exp,
#                         beta_out,
#                         se_out,
#                         beta1,
#                         sigma1,
#                         c,
#                         iter_s,
#                         tol_s) {
#   for (i in 1:iter_s) {
#     sigma0 <- sigma1
#     d <- (beta_out - beta_exp*beta1)^2/(beta1^2*se_exp^2 + se_out^2)
#     w <- rho_biw(d/sigma0, c)*2/p
#     sigma1 <- sum(w*sigma0)
#     if (abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) <= tol_s) {
#       break
#     }
#   }
#   if (abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) > tol_s) {
#     warning("S-estimate warning: IRWLS for sigma-step did not converge. Consider to increase iter_s.")
#   }
#   sigma1
# }
