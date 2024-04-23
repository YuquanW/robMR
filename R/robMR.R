#' Title
#' @import rootSolve
#'
#' @return
#' @export
#'
#' @examples

robMR3 <- function(beta_exp,
                  se_exp,
                  beta_out,
                  se_out,
                  sel_bias = F,
                  pval = 0.05,
                  c = 1.547,
                  k = 4.68,
                  iter_s = 20,
                  iter_mm = 20,
                  tol_s = 1e-5,
                  tol_tau= 1e-3) {
  if (sel_bias == T) {
    lambda <- qnorm(1 - pval/2)
    a_plus <- lambda-beta_exp/se_exp
    a_minus <- -lambda-beta_exp/se_exp
    beta_exp <- beta_exp - se_exp*(dnorm(a_plus)-dnorm(a_minus))/
      (1-pnorm(a_plus)+pnorm(a_minus))
    se_exp <- se_exp*(1 - (a_plus*dnorm(a_plus)-a_minus*dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus))
                      + ((dnorm(a_plus)-dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus)))^2)
  }
  beta_init <- beta_out/beta_exp
  tau_init <- rep(0, length(beta_init))
  for (i in 1:length(beta_init)) {
    f_tau <- function(tau) {
      r <- (beta_out - beta_exp*beta_init[i])/sqrt(beta_init[i]^2*se_exp^2 + tau*se_out^2)
      mean(rho(r, c)) - 0.5
    }
    f_grad_tau <- function(tau) {
      r <- (beta_out - beta_exp*beta_init[i])/sqrt(beta_init[i]^2*se_exp^2 + tau*se_out^2)
      mean(-3/c^2*psi(r, c)*r*se_out^2/(beta_init[i]^2*se_exp^2 + tau*se_out^2))
    }
    tau_init[i] <- nleqslv(1, f_tau, jac = f_grad_tau)$x
  }
  beta1 <- beta_init[which.min(tau_init)]
  tau1 <- min(tau_init)
  print(c(beta1, tau1))
  for (i in 1:iter_s) {
    beta0 <- beta1
    tau0 <- tau1
    f_beta <- function(beta) {
      r <- (beta_out - beta_exp*beta)/sqrt(beta^2*se_exp^2 + tau0*se_out^2)
      mean(psi(r, c)*
             (beta_out*se_exp^2*beta + beta_exp*tau0*se_out^2)/(beta^2*se_exp^2 + tau0*se_out^2)^(3/2))
    }
    beta1 <- multiroot(f_beta, start = beta0, maxiter = 10000)$root
    f_tau <- function(tau) {
      r <- (beta_out - beta_exp*beta0)/sqrt(beta0^2*se_exp^2 + tau*se_out^2)
      mean(rho(r, c)) - 0.5
    }
    f_grad_tau <- function(tau) {
      r <- (beta_out - beta_exp*beta0)/sqrt(beta0^2*se_exp^2 + tau*se_out^2)
      mean(-3/c^2*psi(r, c)*r*se_out^2/(beta0^2*se_exp^2 + tau*se_out^2))
    }
    tau1 <- nleqslv(tau0, f_tau, jac = f_grad_tau)$x
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
  f_beta <- function(beta) {
    r <- (beta_out - beta_exp*beta)/sqrt(beta^2*se_exp^2 + tau1*se_out^2)
    mean(rho(r, k))
  }
  f_grad_beta <- function(beta) {
    r <- (beta_out - beta_exp*beta)/sqrt(beta^2*se_exp^2 + tau1*se_out^2)
    -6/k^2*mean(psi(r, k)*
           (beta_out*se_exp^2*beta + beta_exp*tau1*se_out^2)/(beta^2*se_exp^2 + tau1*se_out^2)^(3/2))
  }
  f_constrain <- function(beta) {
    r0 <- (beta_out - beta_exp*beta1)/sqrt(beta1^2*se_exp^2 + tau1*se_out^2)
    r1 <- (beta_out - beta_exp*beta)/sqrt(beta^2*se_exp^2 + tau1*se_out^2)
    mean(rho(r1, k)) - mean(rho(r0, k))
  }
  f_jac_constrain <- function(beta) {
    r <- (beta_out - beta_exp*beta)/sqrt(beta^2*se_exp^2 + tau1*se_out^2)
    -6/k^2*mean(psi(r, k)*
           (beta_out*se_exp^2*beta + beta_exp*tau1*se_out^2)/(beta^2*se_exp^2 + tau1*se_out^2)^(3/2))
  }
  beta_mm <- nloptr(x0 = beta1,
                    eval_f = f_beta,
                    eval_grad_f = f_grad_beta,
                    lb = 0,
                    ub = 50,
                    eval_g_ineq = f_constrain,
                    eval_jac_g_ineq = f_jac_constrain,
                    opts = list("algorithm" = "NLOPT_LD_MMA",
                                "xtol_rel"=1.0e-8,
                                "print_level" = 0,
                                "check_derivatives" = TRUE,
                                "check_derivatives_print" = "errors"))$solution
  return(list(beta_mm = beta_mm,
              beta_s = beta1,
              tau_s = tau1,
              r_hat_mm = (beta_out - beta_exp*beta_mm)/(tau1*sqrt(beta_mm*se_exp^2 + se_out^2)),
              r_hat_s = (beta_out - beta_exp*beta1)/(tau1*sqrt(beta1*se_exp^2 + se_out^2))))
}


robMR <- function(beta_exp,
                   se_exp,
                   beta_out,
                   se_out,
                   sel_bias = F,
                   c = 2.3952,
                   k = 21.9508,
                   pval = 0.05,
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
  beta_init <- beta_out/beta_exp
  eta_init <- rep(0, length(beta_init))
  for (i in 1:length(beta_init)) {
    f_eta <- function(eta) {
      r <- (beta_out - beta_exp*beta_init[i])^2/(se_exp^2*beta_init[i]^2 + eta*se_out^2)
      mean(rho_biw(r, c)) - 0.5
    }
    f_eta_grad <- function(eta) {
      r <- (beta_out - beta_exp*beta_init[i])^2/(se_exp^2*beta_init[i]^2 + eta*se_out^2)
      -mean(psi_biw(r, c)*r*se_out^2/(se_exp^2*beta_init[i]^2 + eta*se_out^2))
    }
    eta_init[i] <- nleqslv(1, f_eta, jac = f_eta_grad)$x
  }
  beta1 <- beta_init[which.min(eta_init)]
  eta1 <- min(eta_init)
  print(c(beta1, eta1))
  for (i in 1:iter_s) {
    beta0 <- beta1
    eta0 <- eta1
    f_beta <- function(beta) {
      r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + eta0*se_out^2)
      mean(psi_biw(r, c)*(beta_out - beta_exp*beta)*
             (beta_out*se_exp^2*beta + beta_exp*eta0*se_out^2)/(beta^2*se_exp^2 + eta0*se_out^2)^2)
    }
    beta1 <- multiroot(f_beta, start = beta0, maxiter = 10000)$root
    f_eta <- function(eta) {
      r <- (beta_out - beta_exp*beta0)^2/(beta0^2*se_exp^2 + eta*se_out^2)
      mean(rho_biw(r, c)) - 0.5
    }
    f_eta_grad <- function(eta) {
      r <- (beta_out - beta_exp*beta0)^2/(beta0^2*se_exp^2 + eta*se_out^2)
      -mean(psi_biw(r, c)*r*se_out^2/(se_exp^2*beta0^2 + eta*se_out^2))
    }
    eta1 <- nleqslv(eta0, f_eta, jac = f_eta_grad)$x
    print(c(beta1, eta1))
    if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
        abs(eta1 - eta0)/abs(eta1 + 1e-10) <=
        tol_s) {
      break
    }
  }

  if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
      abs(eta1 - eta0)/abs(eta1 + 1e-10) >
      tol_s) {
    warning("S-estimate warning: IRWLS did not converge. Consider to increase iter_s.")
  }
  f_beta <- function(beta) {
    r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + eta1*se_out^2)
    mean(rho_biw(r, k))
  }
  f_beta_grad <- function(beta) {
    r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + eta1*se_out^2)
    -2*mean(psi_biw(r, k)*(beta_out - beta_exp*beta)*
            (beta_out*se_exp^2*beta + beta_exp*eta1*se_out^2)/(beta^2*se_exp^2 + eta1*se_out^2)^2)
  }
  f_constrain <- function(beta) {
    r0 <- (beta_out - beta_exp*beta1)^2/(beta1^2*se_exp^2 + eta1*se_out^2)
    r1 <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + eta1*se_out^2)
    mean(rho_biw(r1, k)) - mean(rho_biw(r0, k))
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
  return(list(beta_mm = beta_mm,
              beta_s = beta1,
              tau_s = eta1,
              r_hat_mm = (beta_out - beta_exp*beta_mm)/sqrt(beta_mm*se_exp^2 + eta1*se_out^2),
              r_hat_s = (beta_out - beta_exp*beta1)/sqrt(beta1*se_exp^2 + eta1*se_out^2)))
}

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
beta_out2 <- bmi.bmi$beta.outcome+c(rnorm(floor(812*0.3), 0, 10/sqrt(812)),
                                    rep(0, 812 - floor(0.3*812)))

index <- which(bmi.bmi$pval.selection > 5e-6)
beta_exp <- beta_exp[index]
beta_out <- beta_out[index]
se_exp <- se_exp[index]
se_out <- se_out[index]
beta_out2 <- beta_out+c(rnorm(floor(length(index)*0.48), 0, 100/sqrt(length(index))),
                rep(0, length(index) - floor(0.48*length(index))))
robMR(beta_exp,
      se_exp,
      beta_out,
      se_out,
      iter_s = 1000,
      iter_mm = 1000)$beta_mm
robMR2(beta_exp,
      se_exp,
      beta_out,
      se_out,
      iter_s = 1000,
      iter_mm = 1000)$beta_mm
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


mr.raps.mle.all(beta_exp, beta_out2, se_exp, se_out)
mr_lasso(mr_input(bx = beta_exp, bxse = se_exp, by = beta_out2, byse = se_out))
mr.divw(beta_exp,
        beta_out2,
        se_exp,
        se_out)
