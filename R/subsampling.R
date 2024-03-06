library(dplyr)
subsamp <- function(beta_exp,
                    se_exp,
                    beta_out,
                    se_out,
                    tau_init) {
  dt <- data.frame(beta_exp, se_exp, beta_out, se_out, tau_init)
  p <- length(beta_exp)
  xi <- beta_exp/se_exp
  xi_1 <- min(xi)
  xi_p <- max(xi)
  t <- seq(xi_1 + (xi_p - xi_1)/(2*p), xi_p - (xi_p - xi_1)/(2*p), (xi_p - xi_1)/p)
  prob <- c(pnorm(t[1]), pnorm(t[2:(p-1)]) - pnorm(t[1:(p-2)]), 1 - pnorm(t[p-1]))

  Nsub <- min(count/prob)
  Nsub_i <- floor(Nsub*prob)
  index <- which(xi <= t[1])
  for (i in 1:p) {
    ind <- which(xi)
  }
}
