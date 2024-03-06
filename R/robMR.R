robMR <- function(beta_exp,
                  se_exp,
                  beta_out,
                  se_out,
                  beta_s,
                  tau_s,
                  mroots = F,
                  c = 4.68,
                  iter_mm = 20,
                  tol_mm = 1e-5
                  ) {
  beta_mm = irwls_semr_beta(beta_exp,
                            se_exp,
                            beta_out,
                            se_out,
                            mroots,
                            beta_s,
                            tau_s,
                            c,
                            iter_mm,
                            tol_mm)
  list(beta_mm = beta_mm,
       r_hat = (beta_out - beta_exp*beta_mm)/sqrt(beta_mm^2*se_exp^2 + se_out^2))
}
