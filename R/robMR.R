#' Title
#' @import nloptr
#'
#' @return
#' @export
#'
#' @examples

robMR <- function(beta_exp,
                  se_exp,
                  beta_out,
                  se_out,
                  sel_bias = F,
                  c = 2.395,
                  lambda = 1.96,
                  iter_s = 20,
                  tol_s = 1e-5) {
  if (sel_bias == T) {
    lambda <- qnorm(1 - pval/2)
    a_plus <- lambda-beta_exp/se_exp
    a_minus <- -lambda-beta_exp/se_exp
    beta_exp <- beta_exp - se_exp*(dnorm(a_plus)-dnorm(a_minus))/
      (1-pnorm(a_plus)+pnorm(a_minus))
    se_exp <- se_exp*(1 - (a_plus*dnorm(a_plus)-a_minus*dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus))
                      + ((dnorm(a_plus)-dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus)))^2)
  }

  beta1 <- median(beta_out/beta_exp)
  f_sigma <- function(sigma) {
    r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
    (mean(rho_biw(r/sigma, c)) - 0.5)^2
  }
  f_sigma_grad <- function(sigma) {
    r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
    -2*(mean(rho_biw(r/sigma, c)) - 0.5)*mean(psi_biw(r/sigma, c)*r/sigma^2)
  }
  sigma1 <- nloptr(x0 = 1,
                 eval_f = f_sigma,
                 eval_grad_f = f_sigma_grad,
                 opts = list("algorithm" = "NLOPT_LD_MMA",
                             "xtol_rel"=1.0e-8,
                             "print_level" = 0))$solution
  for (i in 1:iter_s) {
    beta0 <- beta1
    sigma0 <- sigma1
    f_beta <- function(beta) {
      r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
      mean(rho_biw(r/sigma0, c))
    }
    f_beta_grad <- function(beta) {
      r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
      -2/sigma0*mean(psi_biw(r/sigma0, c)*(beta_out - beta_exp*beta)*
                     (beta_out*se_exp^2*beta + beta_exp*se_out^2)/(beta^2*se_exp^2 + se_out^2)^2)
    }
    beta1 <- nloptr(x0 = beta0,
                    eval_f = f_beta,
                    eval_grad_f = f_beta_grad,
                    opts = list("algorithm" = "NLOPT_LD_MMA",
                                "xtol_rel"=1.0e-8,
                                "print_level" = 0))$solution

    f_sigma <- function(sigma) {
      r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
      (mean(rho_biw(r/sigma, c)) - 0.5)^2
    }
    f_sigma_grad <- function(sigma) {
      r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
      -2*(mean(rho_biw(r/sigma, c)) - 0.5)*mean(psi_biw(r/sigma, c)*r/sigma^2)
    }
    sigma1 <- nloptr(x0 = 1,
                     eval_f = f_sigma,
                     eval_grad_f = f_sigma_grad,
                     opts = list("algorithm" = "NLOPT_LD_MMA",
                                 "xtol_rel"=1.0e-8,
                                 "print_level" = 0))$solution
    if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
        abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) <=
        tol_s) {
      break
    }
  }

  if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
      abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) >
      tol_s) {
    warning("S-estimate warning: IRWLS did not converge. Consider to increase iter_s.")
  }
  c1 <- integrate(function(x)ppsi_biw(x, c)*x*dchisq(x, df = 1), -Inf, Inf)$value
  c2 <- integrate(function(x)psi_biw(x, c)^2*x*dchisq(x, df = 1), -Inf, Inf)$value
  c3 <- integrate(function(x)psi_biw(x, c)*dchisq(x, df = 1), -Inf, Inf)$value
  c4 <- integrate(function(x)psi_biw(x, c)*x*dchisq(x, df = 1), -Inf, Inf)$value
  c5 <- (beta_out^2*se_exp^4*beta1^2 + beta_exp^2*se_out^4 + 2*beta_out*beta_exp*se_exp^2*se_out^2*beta1)/
    (beta1^2*se_exp^2 + se_out^2)^3
  A <- sum(c2*c5)
  B <- sum(-2*c1*c5 - c3*c5 + c4*se_exp^2*se_out^2/(beta1^2*se_exp^2 + se_out^2)^2)
  se1 <- sqrt(A)/abs(B)

  f_k <- function(k) {
    r <- (beta_out - beta_exp*beta1)^2/(beta1^2*se_exp^2 + se_out^2)
    mean(rho_biw(r/sigma1, k))/mean(rho_biw(beta_exp^2/se_exp^2/sigma1, k))
  }
  k <- optimize(f_k, lower = c, upper = 21.9508)$minimum
  cat("The tuning parameter k is set to", k, "\n")

  f_beta <- function(beta) {
    r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
    mean(rho_biw(r/sigma1, k))
  }
  f_beta_grad <- function(beta) {
    r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
    -2/sigma1*mean(psi_biw(r/sigma1, k)*(beta_out - beta_exp*beta)*
                   (beta_out*se_exp^2*beta + beta_exp*se_out^2)/(beta^2*se_exp^2 + se_out^2)^2)
  }
  f_constrain <- function(beta) {
    r0 <- (beta_out - beta_exp*beta1)^2/(beta1^2*se_exp^2 + se_out^2)
    r1 <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
    mean(rho_biw(r1/sigma1, k)) - mean(rho_biw(r0/sigma1, k))
  }
  beta_mm <- nloptr(x0 = beta1,
                    eval_f = f_beta,
                    eval_grad_f = f_beta_grad,
                    eval_g_ineq = f_constrain,
                    eval_jac_g_ineq = f_beta_grad,
                    opts = list("algorithm" = "NLOPT_LD_MMA",
                                "xtol_rel"=1.0e-8,
                                "print_level" = 0,
                                "check_derivatives" = TRUE,
                                "check_derivatives_print" = "errors"))$solution
  c1 <- integrate(function(x)ppsi_biw(x, k)*x*dchisq(x, df = 1), -Inf, Inf)$value
  c2 <- integrate(function(x)psi_biw(x, k)^2*x*dchisq(x, df = 1), -Inf, Inf)$value
  c3 <- integrate(function(x)psi_biw(x, k)*dchisq(x, df = 1), -Inf, Inf)$value
  c4 <- integrate(function(x)psi_biw(x, k)*x*dchisq(x, df = 1), -Inf, Inf)$value
  c5 <- (beta_out^2*se_exp^4*beta_mm^2 + beta_exp^2*se_out^4 + 2*beta_out*beta_exp*se_exp^2*se_out^2*beta_mm)/
    (beta_mm^2*se_exp^2 + se_out^2)^3

  A <- sum(c2*c5)
  B <- sum(-2*c1*c5 - c3*c5 + c4*se_exp^2*se_out^2/(beta_mm^2*se_exp^2 + se_out^2)^2)
  se_mm <- sqrt(A)/abs(B)
  return(list(beta_mm = beta_mm,
              se_mm = se_mm,
              z_mm = beta_mm/se_mm,
              pval_mm = 2*(1 - pnorm(abs(beta_mm/se_mm))),
              ci_mm = c(beta_mm + qnorm(0.025)*se_mm, beta_mm - qnorm(0.025)*se_mm),
              beta_s = beta1,
              sigma_s = sigma1,
              se_s = se1,
              z_s = beta1/se1,
              pval_s = 2*(1 - pnorm(abs(beta1/se1))),
              ci_s = c(beta1 + qnorm(0.025)*se1, beta1 - qnorm(0.025)*se1),
              r_hat_mm = (beta_out - beta_exp*beta_mm)/sqrt(sigma1*beta_mm^2*se_exp^2 + sigma1*se_out^2),
              r_hat_s = (beta_out - beta_exp*beta1)/sqrt(sigma1*beta1^2*se_exp^2 + sigma1*se_out^2),
              w = psi_biw((beta_out - beta_exp*beta_mm)^2/(sigma1*beta_mm^2*se_exp^2 + sigma1*se_out^2), k)*k/3))
}

library(MRcML)
library(rootSolve)
library(nleqslv)
library(nloptr)
library(mr.raps)
library(mr.divw)
library(MendelianRandomization)


beta_exp <- bmi.bmi$beta.exposure
beta_out <- bmi.bmi$beta.outcome
se_exp <- bmi.bmi$se.exposure
se_out <- bmi.bmi$se.outcome
beta_out2 <- bmi.bmi$beta.outcome+c(rnorm(floor(812*0.48), 0, 1/sqrt(812)),
                                    rep(0, 812 - floor(0.48*812)))

index <- which(bmi.bmi$pval.selection > 5e-5)
beta_exp <- beta_exp[index]
beta_out <- beta_out[index]
se_exp <- se_exp[index]
se_out <- se_out[index]
beta_out2 <- beta_out+c(rnorm(floor(length(index)*0.48), 0, 1/sqrt(length(index))),
                rep(0, length(index) - floor(0.48*length(index))))
robMR(beta_exp,
      se_exp,
      beta_out,
      se_out,
      iter_s = 1000)$beta_mm
robMR2(beta_exp,
      se_exp,
      beta_out,
      se_out,
      iter_s = 1000)$beta_mm
robMR3(beta_exp,
      se_exp,
      beta_out,
      se_out,
      iter_s = 1000,
      iter_mm = 1000)$beta_mm
robMR4(beta_exp,
       se_exp,
       beta_out,
       se_out,
       iter_s = 1000)$beta_cm


mr.raps.mle.all(beta_exp, beta_out, se_exp, se_out)
mr_lasso(mr_input(bx = beta_exp, bxse = se_exp, by = beta_out, byse = se_out))
mr.divw(beta_exp,
        beta_out,
        se_exp,
        se_out)

n <- nrow(bmi.bmi)
res_mm <- matrix(0, 100, 2)
res_raps <- matrix(0, 100, 2)
res_divw <- matrix(0, 100, 2)
for (i in 1:100) {
  beta_out_simu <- rnorm(n, mean = 0+c(rnorm(floor(n*0.4), 0, 1/sqrt(n)),
                                      rep(0, n - floor(0.4*n))), sd = se_out)
  beta_exp_simu <- rnorm(n, mean = beta_exp, sd = se_exp)
  raps <- mr.raps.overdispersed.robust(beta_exp_simu, beta_out_simu, se_exp, se_out)
  res_raps[i, 1] <- raps$beta.hat
  res_raps[i, 2] <- raps$beta.p.value
  mm <- robMR(beta_exp_simu,
              se_exp,
              beta_out_simu,
              se_out,
              iter_s = 1000)
  res_mm[i, 1] <- mm$beta_mm
  res_mm[i, 2] <- mm$pval_mm
  divw <- mr.divw(beta_exp_simu, beta_out_simu, se_exp, se_out)
  res_divw[i, 1] <- raps$beta.hat
  res_divw[i, 2] <- raps$beta.p.value
}

n <- 10000
p <- 200
r <- 100
maf <- runif(p, 0.1, 0.5)
G <- matrix(0, 2*n, p)
for (j in 1:p) {
  G[, j] <- rbinom(2*n, 2, maf[j])
}
G1 <- G[1:n,]
G2 <- G[(n+1):(2*n),]
gamma <- rnorm(p, 0, 1/sqrt(p))
beta <- 0

res_mm <- matrix(0, r, 2)
res_raps <- matrix(0, r, 2)
for (i in 1:r) {
  X1 <- drop(G1%*%gamma) + rnorm(n, sd = 2)
  X2 <- drop(G2%*%gamma) + rnorm(n, sd = 2)
  Y2 <- X2*beta + rnorm(n)
  beta_exp <- beta_out <- se_exp <- se_out <- rep(0, p)
  for (j in 1:p) {
    GWAS1 <- lm(X1~G1[, j])
    GWAS2 <- lm(Y2~G2[, j])
    beta_exp[j] <- coef(summary(GWAS1))[2, 1]
    se_exp[j] <- coef(summary(GWAS1))[2, 2]
    beta_out[j] <- coef(summary(GWAS2))[2, 1]
    se_out[j] <- coef(summary(GWAS2))[2, 2]
  }
  raps <- mr.raps.overdispersed.robust(beta_exp, beta_out, se_exp, se_out)
  res_raps[i, 1] <- raps$beta.hat
  res_raps[i, 2] <- raps$beta.p.value
  mm <- robMR(beta_exp,
              se_exp,
              beta_out,
              se_out,
              iter_s = 1000)
  res_mm[i, 1] <- mm$beta_mm
  res_mm[i, 2] <- mm$pval_mm
}
