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
                  c = 2.395,
                  adjust_eff = T,
                  sel_bias = F,
                  lambda = 1.96,
                  iter_s = 20,
                  tol_s = 1e-3,
                  iter_mm = 20,
                  tol_mm = 1e-5,
                  print = F) {
  p <- length(beta_exp)
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
  # f_sigma <- function(sigma) {
  #   r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
  #   (mean(rho_biw(r/sigma, c)) - 0.5)^2
  # }
  # f_sigma_grad <- function(sigma) {
  #   r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
  #   -2*(mean(rho_biw(r/sigma, c)) - 0.5)*mean(psi_biw(r/sigma, c)*r/sigma^2)
  # }
  sigma1 <- irwls_sigma(beta_exp, se_exp, beta_out, se_out, beta1, 1, c, iter_s, tol_s)
    # nloptr(x0 = 1,
    #                eval_f = f_sigma,
    #                eval_grad_f = f_sigma_grad,
    #                opts = list("algorithm" = "NLOPT_LD_MMA",
    #                            "xtol_rel"=1.0e-8,
    #                            "print_level" = 0))$solution
  for (i in 1:iter_s) {
    beta0 <- beta1
    sigma0 <- sigma1
    # f_beta <- function(beta) {
    #   r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
    #   mean(rho_biw(r/sigma0, c))
    # }
    # f_beta_grad <- function(beta) {
    #   r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
    #   -2/sigma0*mean(psi_biw(r/sigma0, c)*(beta_out - beta_exp*beta)*
    #                    (beta_out*se_exp^2*beta + beta_exp*se_out^2)/(beta^2*se_exp^2 + se_out^2)^2)
    # }
    beta1 <- irwls_beta(beta_exp, se_exp, beta_out, se_out, beta0, sigma0, c, iter_s, tol_s)
      # nloptr(x0 = beta0,
      #               eval_f = f_beta,
      #               eval_grad_f = f_beta_grad,
      #               opts = list("algorithm" = "NLOPT_LD_MMA",
      #                           "xtol_rel"=1.0e-8,
      #                           "print_level" = 0))$solution

    # f_sigma <- function(sigma) {
    #   r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
    #   (mean(rho_biw(r/sigma, c)) - 0.5)^2
    # }
    # f_sigma_grad <- function(sigma) {
    #   r <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
    #   -2*(mean(rho_biw(r/sigma, c)) - 0.5)*mean(psi_biw(r/sigma, c)*r/sigma^2)
    # }
    sigma1 <- irwls_sigma(beta_exp, se_exp, beta_out, se_out, beta0, sigma0, c, iter_s, tol_s)
      # nloptr(x0 = 1,
      #                eval_f = f_sigma,
      #                eval_grad_f = f_sigma_grad,
      #                opts = list("algorithm" = "NLOPT_LD_MMA",
      #                            "xtol_rel"=1.0e-8,
      #                            "print_level" = 0))$solution
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

  if (adjust_eff == T) {
    f_k <- function(k) {
      r <- (beta_out - beta_exp*beta1)^2/(beta1^2*se_exp^2 + se_out^2)
      mean(rho_biw(r/sigma1, k))/mean(rho_biw(beta_exp^2/se_exp^2/sigma1, k))
    }
    k <- optimize(f_k, lower = c, upper = 21.9508)$minimum
    cat("The tuning parameter k is set to", round(k, 4), "after IVstrength-efficiency trade-off.\n")
  } else {
    k <- 21.9508
  }


  # f_beta <- function(beta) {
  #   r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
  #   mean(rho_biw(r/sigma1, k))
  # }
  # f_beta_grad <- function(beta) {
  #   r <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
  #   -2/sigma1*mean(psi_biw(r/sigma1, k)*(beta_out - beta_exp*beta)*
  #                    (beta_out*se_exp^2*beta + beta_exp*se_out^2)/(beta^2*se_exp^2 + se_out^2)^2)
  # }
  # f_constrain <- function(beta) {
  #   r0 <- (beta_out - beta_exp*beta1)^2/(beta1^2*se_exp^2 + se_out^2)
  #   r1 <- (beta_out - beta_exp*beta)^2/(beta^2*se_exp^2 + se_out^2)
  #   mean(rho_biw(r1/sigma1, k)) - mean(rho_biw(r0/sigma1, k))
  # }
  beta_mm <- irwls_beta(beta_exp, se_exp, beta_out, se_out, beta1, sigma1, k, iter_mm, tol_mm)
    # nloptr(x0 = beta1,
    #                 eval_f = f_beta,
    #                 eval_grad_f = f_beta_grad,
    #                 eval_g_ineq = f_constrain,
    #                 eval_jac_g_ineq = f_beta_grad,
    #                 opts = list("algorithm" = "NLOPT_LD_MMA",
    #                             "xtol_rel"=1.0e-8,
    #                             "print_level" = 0,
    #                             "check_derivatives" = print,
    #                             "check_derivatives_print" = ifelse(print, "error", "none")))$solution

  d_s <- (beta_out - beta_exp*beta1)^2/(se_exp^2*beta1^2 + se_out^2)
  t_s <- (beta_out - beta_exp*beta1)/sqrt(se_exp^2*beta1^2 + se_out^2)
  s_s <- (beta1*se_exp^2*beta_out+se_out^2*beta_exp)/(beta1^2*se_exp^2+se_out^2)^(3/2)

  d_mm <- (beta_out - beta_exp*beta_mm)^2/(se_exp^2*beta_mm^2 + se_out^2)
  t_mm <- (beta_out - beta_exp*beta_mm)/sqrt(se_exp^2*beta_mm^2 + se_out^2)
  s_mm <- (beta_mm*se_exp^2*beta_out+se_out^2*beta_exp)/(beta_mm^2*se_exp^2+se_out^2)^(3/2)
  A <- crossprod(cbind(psi_biw(d_s/sigma1, c)*t_s*s_s, rho_biw(d_s/sigma1, c) - 0.5, psi_biw(d_mm/sigma1, k)*t_mm*s_mm))#sum(psi_biw(d/sigma1, c)^2*d*s^2)


  B11 <- sum(2*ppsi_biw(d_s/sigma1, c)*d_s/sigma1*s_s^2+
               psi_biw(d_s/sigma1, c)*s_s^2+
               2*beta1*se_exp^2/(beta1^2*se_exp^2+se_out^2)*psi_biw(d_s/sigma1, c)*t_s*s_s-
               se_exp^2*se_out^2/(beta1^2*se_exp^2+se_out^2)^2*psi_biw(d_s/sigma1, c)*d_s)
  B12 <- sum(ppsi_biw(d_s/sigma1, c)*d_s/sigma1^2*t_s*s_s)
  B13 <- 0
  B21 <- 0
  B22 <- sum(psi_biw(d_s/sigma1, c)*d_s/sigma1^2)
  B23 <- 0
  B31 <- 0
  B32 <- sum(ppsi_biw(d_mm/sigma1, k)*d_mm/sigma1^2*t_mm*s_mm)
  B33 <- sum(2*ppsi_biw(d_mm/sigma1, k)*d_mm/sigma1*s_mm^2+
               psi_biw(d_mm/sigma1, k)*s_mm^2+
               2*beta1*se_exp^2/(beta_mm^2*se_exp^2+se_out^2)*psi_biw(d_mm/sigma1, k)*t_mm*s_mm-
               se_exp^2*se_out^2/(beta_mm^2*se_exp^2+se_out^2)^2*psi_biw(d_mm/sigma1, k)*d_mm)
  B <- matrix(c(B11, B12, B13, B21, B22, B23, B31, B32, B33), 3, 3, byrow = T)
  C <- solve(B, A)%*%t(solve(B))
  se1 <- sqrt(C[1,1])

  # d <- (beta_out - beta_exp*beta_mm)^2/(se_exp^2*beta_mm^2 + se_out^2)
  # t <- (beta_out - beta_exp*beta_mm)/sqrt(se_exp^2*beta_mm^2 + se_out^2)
  # s <- (beta_mm*se_exp^2*beta_out+se_out^2*beta_exp)/(beta_mm^2*se_exp^2+se_out^2)^(3/2)
  # A <- sum(psi_biw(d/sigma1, k)^2*d*s^2)
  # B <- sum(2*ppsi_biw(d/sigma1, k)*d/sigma1*s^2+
  #            psi_biw(d/sigma1, k)*s^2+
  #            2*beta_mm*se_exp^2/(beta_mm^2*se_exp^2+se_out^2)*psi_biw(d/sigma1, k)*t*s-
  #            se_exp^2*se_out^2/(beta_mm^2*se_exp^2+se_out^2)^2*psi_biw(d/sigma1, k)*d)
  se_mm <- sqrt(C[3,3])
  print(sqrt(C[2,2]))
  # beta_mm_b <- rep(0, 10000)
  # sigma_s_b <- rep(0, 10000)
  # beta_s_b <- rep(0, 10000)
  # for (i in 1:10000) {
  #   index <- sample(p, p, replace = T)
  #   beta_exp_b <- beta_exp[index]
  #   se_exp_b <- se_exp[index]
  #   beta_out_b <- beta_out[index]
  #   se_out_b <- se_out[index]
  #   d_b <- (beta_out_b - beta_exp_b*beta_mm)^2/(se_exp_b^2*beta_mm^2 + se_out_b^2)
  #   w_b <- psi_biw(d_b/sigma1, k)
  #   s_b <- (beta_mm*se_exp_b^2*beta_out_b+se_out_b^2*beta_exp_b)/(beta_mm^2*se_exp_b^2+se_out_b^2)^(3/2)
  #   t_b <- (beta_out_b-beta_exp_b*beta_mm)/sqrt(beta_mm^2*se_exp_b^2+se_out_b^2)
  #   beta_mm_tmp <- sum(w_b*s_b*beta_out_b/sqrt(beta_mm^2*se_exp_b^2+se_out_b^2))/
  #     sum(w_b*s_b*beta_exp_b/sqrt(beta_mm^2*se_exp_b^2+se_out_b^2))
  #
  #   d_t <- (beta_out_b - beta_exp_b*beta1)^2/(se_exp_b^2*beta1^2 + se_out_b^2)
  #   w_t <- psi_biw(d_t/sigma1, c)
  #   s_t <- (beta1*se_exp_b^2*beta_out_b+se_out_b^2*beta_exp_b)/(beta1^2*se_exp_b^2+se_out_b^2)^(3/2)
  #   t_t <- (beta_out_b-beta_exp_b*beta1)/sqrt(beta1^2*se_exp_b^2+se_out_b^2)
  #   sigma_s_tmp <- 2*sigma1/p*sum(rho_biw(d_t/sigma1, c))
  #   beta_s_tmp <- sum(w_t*s_t*beta_out_b/sqrt(beta1^2*se_exp_b^2+se_out_b^2))/
  #     sum(w_t*s_t*beta_exp_b/sqrt(beta1^2*se_exp_b^2+se_out_b^2))
  #   C11 <- B33/sum(psi_biw(d_mm/sigma1, k)*s_mm*beta_exp/sqrt(beta_mm^2*se_exp^2+se_out^2))
  #   C12 <- B32/sum(psi_biw(d_mm/sigma1, k)*s_mm*beta_exp/sqrt(beta_mm^2*se_exp^2+se_out^2))
  #   C13 <- 0
  #   C21 <- 0
  #   C22 <- sum(psi_biw(d_s/sigma1, c)*d_s)*2/(p*sigma1)
  #   C23 <- 0
  #   C31 <- 0
  #   C32 <- B12/sum(psi_biw(d_s/sigma1, c)*s_s*beta_exp/sqrt(beta1^2*se_exp^2+se_out^2))
  #   C33 <- B11/sum(psi_biw(d_s/sigma1, c)*s_s*beta_exp/sqrt(beta1^2*se_exp^2+se_out^2))
  #   C <- matrix(c(C11, C12, C13,
  #                 C21, C22, C23,
  #                 C31, C32, C33), 3, 3, byrow = T)
  #   par_b <- c(beta_mm, sigma1, beta1) + drop(solve(C, c(beta_mm_tmp - beta_mm, sigma_s_tmp - sigma1, beta_s_tmp - beta1)))
  #   beta_mm_b[i] <- par_b[1]
  #   sigma_s_b[i] <- par_b[2]
  #   beta_s_b[i] <- par_b[3]
  # }

  return(list(beta_mm = beta_mm,
              se_mm = se_mm,
              z_mm = beta_mm/se_mm,
              pval_mm = 2*(1 - pnorm(abs(beta_mm/se_mm))),
              ci_mm = c(beta_mm + qnorm(0.025)*se_mm, beta_mm - qnorm(0.025)*se_mm),
              #beta_mm_b = beta_mm_b,
              beta_s = beta1,
              sigma_s = sigma1,
              se_s = se1,
              z_s = beta1/se1,
              pval_s = 2*(1 - pnorm(abs(beta1/se1))),
              ci_s = c(beta1 + qnorm(0.025)*se1, beta1 - qnorm(0.025)*se1),
              r_hat_mm = (beta_out - beta_exp*beta_mm)/sqrt(sigma1*beta_mm^2*se_exp^2 + sigma1*se_out^2),
              r_hat_s = (beta_out - beta_exp*beta1)/sqrt(sigma1*beta1^2*se_exp^2 + sigma1*se_out^2),
              w_s = psi_biw(d_s/sigma1, c)*c/3,
              w_mm = psi_biw(d_mm/sigma1, k)*k/3))
}

robMR <- function(beta_exp,
                  se_exp,
                  beta_out,
                  se_out,
                  c = 1.56,
                  sel_bias = F,
                  lambda = 1.96,
                  iter_s = 20,
                  tol_s = 1e-3,
                  iter_mm = 20,
                  tol_mm = 1e-5,
                  print = F) {
  p <- length(beta_exp)
  if (sel_bias == T) {
    lambda <- qnorm(1 - pval/2)
    a_plus <- lambda-beta_exp/se_exp
    a_minus <- -lambda-beta_exp/se_exp
    beta_exp <- beta_exp - se_exp*(dnorm(a_plus)-dnorm(a_minus))/
      (1-pnorm(a_plus)+pnorm(a_minus))
    se_exp <- se_exp*(1 - (a_plus*dnorm(a_plus)-a_minus*dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus))
                      + ((dnorm(a_plus)-dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus)))^2)
  }

  #beta1 <- median(beta_out/beta_exp)
  #beta1 <- irwls_beta_init(beta_exp, se_exp, beta_out, se_out, beta1, iter_s, tol_s)
  f_init <- function(beta) {
    median((beta_out-beta_exp*beta)^2/(beta^2*se_exp^2+se_out^2))
  }
  beta1 <- nloptr(median(beta_out/beta_exp), f_init, opts = list("algorithm" = "NLOPT_LN_BOBYQA",
                                                                 "xtol_rel"=1.0e-8,
                                                                 "print_level" = 0))$x0
  sigma1 <- irwls_sigma(beta_exp, se_exp, beta_out, se_out, beta1, 1, c, iter_s, tol_s)
  # for (i in 1:iter_s) {
  #   beta0 <- beta1
  #   sigma0 <- sigma1
  #   beta1 <- irwls_beta(beta_exp, se_exp, beta_out, se_out, beta0, sigma0, c, iter_s, tol_s)
  #   sigma1 <- irwls_sigma(beta_exp, se_exp, beta_out, se_out, beta0, sigma0, c, iter_s, tol_s)
  #   if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
  #       abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) <=
  #       tol_s) {
  #     break
  #   }
  # }
  #
  # if (abs(beta1 - beta0)/abs(beta1 + 1e-10) +
  #     abs(sigma1 - sigma0)/abs(sigma1 + 1e-10) >
  #     tol_s) {
  #   warning("S-estimate warning: IRWLS did not converge. Consider to increase iter_s.")
  # }

  q_s <- se_exp*se_out/(beta1^2*se_exp^2+se_out^2)
  s_s <- (beta_out*beta1*se_exp^2+beta_exp*se_out^2)/(se_exp*se_out*sqrt(beta1^2*se_exp^2+se_out^2))
  r_s <- (beta_out-beta1*beta_exp)/sqrt(beta1^2*se_exp^2+se_out^2)
  f <- function(k) {
    k*max(q_s)/abs(sum(1/sigma1*ppsi(abs(s_s)*r_s/sigma1, k)*(r_s^2-s_s^2)*q_s^2))
  }
  k <- direct(f, lower = 1.56, upper = 30)$par

  beta_mm <- irwls_beta(beta_exp, se_exp, beta_out, se_out, beta1, sigma1, k, iter_mm, tol_mm)
  s_mm <- (beta_out*beta_mm*se_exp^2+beta_exp*se_out^2)/(se_exp*se_out*sqrt(beta_mm^2*se_exp^2+se_out^2))
  r_mm <- (beta_out-beta_mm*beta_exp)/sqrt(beta_mm^2*se_exp^2+se_out^2)
  q_mm <- se_exp*se_out/(beta_mm^2*se_exp^2+se_out^2)
  w_mm <- wgt_m(abs(s_mm)*r_mm/sigma1, k)
  # A <- crossprod(cbind(psi(abs(s_s)*r_s/sigma1, c)*s_s/abs(s_s)*q_s,
  #                      rho(r_s/sigma1, c)-1/2,
  #                      psi(abs(s_mm)*r_mm/sigma1, k)*s_mm/abs(s_mm)*q_mm))
  # B <- matrix(c(sum(1/sigma1*ppsi(abs(s_s)*r_s/sigma1, c)*(r_s^2-s_s^2)*q_s^2), -sum(1/sigma1^2*ppsi(abs(s_s)*r_s/sigma1, c)*r_s*s_s*q_s), 0,
  #               -sum(psi(r_s/sigma1, c)*s_s*q_s), -sum(psi(r_s/sigma1, c)*r_s/sigma1^2), 0,
  #               0, -sum(1/sigma1^2*ppsi(abs(s_mm)*r_mm/sigma1, k)*r_mm*s_mm*q_mm), sum(1/sigma1*ppsi(abs(s_mm)*r_mm/sigma1, k)*(r_mm^2-s_mm^2)*q_mm^2)), byrow = T, nrow = 3)
  # V <- solve(B, A)%*%t(solve(B))
  beta_mm_b <- rep(0, 10000)
  sigma_s_b <- rep(0, 10000)
  for (i in 1:10000) {
    index <- sample(p, p, replace = T)
    beta_exp_b <- beta_exp[index]
    se_exp_b <- se_exp[index]
    beta_out_b <- beta_out[index]
    se_out_b <- se_out[index]
    r_b <- (beta_out_b-beta_mm*beta_exp_b)/sqrt(beta_mm^2*se_exp_b^2+se_out_b^2)
    s_b <- (beta_out_b*beta_mm*se_exp_b^2+beta_exp_b*se_out_b^2)/(se_exp_b*se_out_b*sqrt(beta_mm^2*se_exp_b^2+se_out_b^2))
    w_b <- wgt_m(abs(s_b)*r_b/sigma1, k)
    q_b <- se_exp_b*se_out_b/(beta_mm^2*se_exp_b^2+se_out_b^2)
    beta_mm_tmp <- sum(w_b*s_b*q_b*beta_out_b/sqrt(beta_mm^2*se_exp_b^2+se_out_b^2))/
      sum(w_b*s_b*q_b*beta_exp_b/sqrt(beta_mm^2*se_exp_b^2+se_out_b^2))

    r_sb <- (beta_out_b-beta1*beta_exp_b)/sqrt(beta1^2*se_exp_b^2+se_out_b^2)
    sigma_s_tmp <- 2*sigma1/p*sum(rho(r_sb/sigma1, c))
    C11 <- sum(1/sigma1*ppsi(abs(s_mm)*r_mm/sigma1, k)*(s_mm^2-r_mm^2)*q_mm^2)/
      sum(w_mm*s_mm*q_mm*beta_exp/sqrt(beta_mm^2*se_exp^2+se_out^2))
    C12 <- sum(1/sigma1^2*ppsi(abs(s_mm)*r_mm/sigma1, k)*r_mm*s_mm*q_mm)/
      sum(w_mm*s_mm*q_mm*beta_exp/sqrt(beta_mm^2*se_exp^2+se_out^2))
    C21 <- 0
    C22 <- sum(psi(r_s/sigma1, c)*r_s)*2/(p*sigma1)
    C <- matrix(c(C11, C12,
                  C21, C22), 2, 2, byrow = T)
    par_b <- c(beta_mm, sigma1) + drop(solve(C, c(beta_mm_tmp - beta_mm, sigma_s_tmp - sigma1)))
    beta_mm_b[i] <- par_b[1]
    sigma_s_b[i] <- par_b[2]
  }
  # se_s <- sqrt(V[1, 1])
  # se_mm <- sqrt(V[3, 3])

  return(list(beta_mm = beta_mm,
              #se_mm = se_mm,
              beta_mm.ci = quantile(beta_mm_b, c(0.025, 0.975)),
              beta_s = beta1,
              #se_s = se_s,
              sigma_s = sigma1,
              w_s = wgt_m(abs(s_s)*r_s/sigma1, c),
              w_mm = wgt_m(abs(s_mm)*r_mm/sigma1, k)
              ))
}


library(MRcML)
library(rootSolve)
library(nleqslv)
library(nloptr)
library(mr.raps)
library(mr.divw)
library(MendelianRandomization)

library(TwoSampleMR)
library(stringr)
ao <- available_outcomes()
BMI <- ao[grepl("Body mass index", ao$trait, fixed = TRUE),]

#selection_data <- extract_instruments("ukb-b-19953", p1 = 1e-2, force_server = T)#"ebi-a-GCST90025996", p1 = 5e-8, force_server = T)#
while (T) {
  exposure_data <- try(extract_instruments("ukb-b-19953", p1 = 1), silent = T)
  if (!is(exposure_data, "try-error")) break
}
#exposure_data <- extract_instruments("ukb-b-19953", p1 = 1)#extract_outcome_data(selection_data$SNP, proxies = F, outcomes = "ieu-a-785")
nms <- colnames(exposure_data)
colnames(exposure_data) <- str_replace(nms, "outcome", "exposure")
outcome_data <- extract_outcome_data(exposure_data$SNP, proxies = F, outcomes = "ieu-a-2")#extract_outcome_data(selection_data$SNP, proxies = F, outcomes = "ieu-a-974")
dat <- harmonise_data(exposure_data, outcome_data)

dat <- bmi.bmi
beta_exp <- dat$beta.exposure
beta_out <- dat$beta.outcome
se_exp <- dat$se.exposure
se_out <- dat$se.outcome
beta_out2 <- bmi.bmi$beta.outcome+c(rnorm(floor(812*0.48), 0, 5*mean(se_out)),
                                    rep(0, 812 - floor(0.48*812)))
p <- length(beta_exp)
beta_hat <- rep(0, 100)
sigma_hat <- rep(0, 100)
eta <- 1
for (i in 1:100) {
  beta_exp_mc <- rnorm(p, beta_exp, eta*se_exp)
  beta_out_mc <- rnorm(p, beta_out, eta*se_out)
  mm <- robMR(beta_exp_mc,
        eta*se_exp,
        beta_out_mc,
        eta*se_out,
        iter_s = 1000)
  beta_hat[i] <- mm$beta_mm
  sigma_hat[i] <- mm$sigma_s
}
beta_hat[abs(sigma_hat-1)<quantile(abs(sigma_hat-1), 0.01)]

index <- which(dat$pval.exposure < 5e-5)
beta_exp_sel <- beta_exp[index]
beta_out_sel <- beta_out[index]
se_exp_sel <- se_exp[index]
se_out_sel <- se_out[index]
beta_out2 <- beta_out+c(rnorm(floor(length(index)*0.48), 0, 1/sqrt(length(index))),
                rep(0, length(index) - floor(0.48*length(index))))

n <- min(dat$samplesize.exposure, dat$samplesize.outcome)
p <- length(beta_exp)
eta <- 0.5
Z <- rnorm(p, sd = eta)
pval <- 5*10^seq(-8, -2)
res <- matrix(0, length(pval), 9)
NIV <- rep(0, length(pval))
for (i in 1:length(pval)) {
  lambda <- qnorm(1-pval[i]/2)
  rr <- abs(beta_exp/se_exp+Z) > lambda
  beta_exp_sel <- beta_exp[rr]
  beta_out_sel <- beta_out[rr]
  se_exp_sel <- se_exp[rr]
  se_out_sel <- se_out[rr]
  a_plus <- lambda - beta_exp_sel/se_exp_sel
  a_minus <- -lambda - beta_exp_sel/se_exp_sel
  beta_exp_rb <- beta_exp_sel - se_exp_sel*(dnorm(a_plus)-dnorm(a_minus))/
    (1-pnorm(a_plus)+pnorm(a_minus))
  se_exp_rb <- se_exp_sel*sqrt(1 - (a_plus*dnorm(a_plus)-a_minus*dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus))
                               + ((dnorm(a_plus)-dnorm(a_minus))/(1-pnorm(a_plus)+pnorm(a_minus)))^2)
  res[i, 5] <- mr.divw(beta_exp_rb,
                       beta_out_sel,
                       se_exp_rb,
                       se_out_sel)$beta.hat
  sel <- abs(beta_exp/se_exp) > lambda
  beta_exp_sel <- beta_exp[sel]
  beta_out_sel <- beta_out[sel]
  se_exp_sel <- se_exp[sel]
  se_out_sel <- se_out[sel]
  res[i, 1] <- mr_lasso(mr_input(bx = beta_exp_sel, bxse = se_exp_sel, by = beta_out_sel, byse = se_out_sel))$Estimate
  res[i, 2] <- MendelianRandomization::mr_ivw(mr_input(bx = beta_exp_sel, bxse = se_exp_sel, by = beta_out_sel, byse = se_out_sel), robust = T)$Estimate
  res[i, 3] <- mr_cML(mr_input(bx = beta_exp_sel, bxse = se_exp_sel, by = beta_out_sel, byse = se_out_sel), DP=F, n = n)$Estimate
  res[i, 4] <- mr.divw(beta_exp_sel,
                       beta_out_sel,
                       se_exp_sel,
                       se_out_sel)$beta.hat
  res[i, 6] <- mr.raps.overdispersed.robust(beta_exp_sel,
                                            beta_out_sel,
                                            se_exp_sel,
                                            se_out_sel)$beta.hat
  res[i, 7] <- robMR(beta_exp_sel,
                     se_exp_sel,
                     beta_out_sel,
                     se_out_sel,
                     sel_bias = T,
                     pval = pval[i],
                     iter_s = 1000)$beta_s
  res[i, 8] <- robMR(beta_exp_sel,
                     se_exp_sel,
                     beta_out_sel,
                     se_out_sel,
                     adjust_eff = F,
                     sel_bias = T,
                     pval = pval[i],
                     iter_s = 1000)$beta_mm
  res[i, 9] <- robMR(beta_exp_sel,
                     se_exp_sel,
                     beta_out_sel,
                     se_out_sel,
                     sel_bias = T,
                     pval = pval[i],
                     iter_s = 1000)$beta_mm

}
for (i in 1:length(pval)) {
  lambda <- qnorm(1-pval[i]/2)
  sel <- abs(beta_exp/se_exp) > lambda
  NIV[i] <- sum(sel)
}
colnames(res) <- c("Lasso", "IVW(rob)", "MR-cML", "dIVW", "RIVW", "RAPS", "MR-S", "MR-MM", "MR-MM(adj)")
res_dt <- data.frame(cbind(Pvalue =  pval, res),check.names = F)
res_dt$Pvalue <- as.character(res_dt$Pvalue)
res_dt <- cbind("Pvalue"=res_dt$Pvalue, "\\#IVs" = NIV, res_dt[, 3:11])
kbl(res_dt, "latex", digits = 3, booktabs = TRUE, escape = FALSE, align = "c") %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))



robMR(beta_exp_rb,
      se_exp_rb,
      beta_out_sel,
      se_out_sel,
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
mr.divw(beta_exp_rb,
        beta_out_sel,
        se_exp_rb,
        se_out_sel)

n <- 200
dat_subset <- dat[1:n, ]
n <- nrow(dat_subset)
res_mm <- matrix(0, 100, 4)
res_raps <- matrix(0, 100, 2)
res_divw <- matrix(0, 100, 2)

beta_exp <- dat_subset$beta.exposure
se_exp <- dat_subset$se.exposure
se_out <- dat_subset$se.outcome
alpha <- c(rnorm(floor(n*0.48), 0, 1/sqrt(n)),
           rep(0, n - floor(0.1*n)))
for (i in 1:100) {
  beta_out_simu <- rnorm(n, mean = beta_exp, sd = se_out)+alpha
  beta_exp_simu <- rnorm(n, mean = beta_exp, sd = se_exp)
  mm <- robMR(beta_exp_simu,
              se_exp,
              beta_out_simu,
              se_out,
              iter_s = 1000)
  res_mm[i, 1] <- mm$beta_mm
  res_mm[i, 2] <- mm$se_mm
  res_mm[i, 3] <- mm$beta_mm - 1.96*mm$se_mm
  res_mm[i, 4] <- mm$beta_mm + 1.96*mm$se_mm
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
