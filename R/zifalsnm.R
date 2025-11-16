###### ZIFA LSN Model ######

#' Zero Inflated Factor Analysis Logistic Skew Normal Multinomial (ZIFA - LSNM) model
#' @description ZIFA-LSNM model's objective is to achieve effective dimension
#' reduction to manage high dimensionality, accounting for the microbiome
#' data compositional nature, zero inflation and, it utilizes
#' skew-normal priors on the latent factors to explicitly model skewness
#' in log-ratio transformation.
#' @export
#' @param X Count matrix of samples and taxa.
#' @param num_fac Number of latent factors.
#' @param NU1P Prior value for the shape parameter of \eqn{\kappa_j}.
#' @param NU2P Prior value for the shape parameter of \eqn{\kappa_j}.
#' @param G1P Prior value for the shape parameter of \eqn{\delta_{jt}}.
#' @param G2P Prior value for the rate parameter of \eqn{\delta_{jt}}.
#' @param AlphaP Prior value for the skewness parameter of latent factors \eqn{F_{it}}.
#'
#' @return Returns optimized values
#'
#' @importFrom sn rsn
#' @import stats
#'
#' @examples
#' \donttest{
#' N <- 50
#' P <- 50
#' K <- 2
#'
#' alpha_prior <- 2
#' alpha_star <- -sqrt(2 / pi) * (alpha_prior / sqrt(1 + alpha_prior^2))
#' f <- matrix(sn::rsn(N * K, alpha_star, 1, alpha_prior), N, K, TRUE)
#'
#' beta <- matrix(runif(P * K, -2.5, 2.5), nrow = P, ncol = K, byrow = TRUE)
#'
#' kappa <- runif(P)
#'
#' beta_0 <- runif(P, -2, 2)
#'
#' z <- matrix(0, N, P, byrow = TRUE)
#' for (i in 1:N) {
#'   z[i, ] <- rbinom(P, size = 1, prob = kappa)
#' }
#' linear_predictor <- matrix(beta_0, N, P, byrow = TRUE) + tcrossprod(f, beta)
#' rho <- (1 - z) * exp(linear_predictor) / (rowSums((1 - z) * exp(linear_predictor)) + 1e-8)
#'
#' cmat <- matrix(0, N, P)
#'
#' for (i in 1:N) cmat[i, ] <- rmultinom(1, as.integer(runif(1, 80000, 100000)), prob = rho[i, ])
#'
#' OUTPUT <- ZIFA_LSNM(X = cmat, num_fac = K, NU1P = 1, NU2P = 3, G1P = 200, G2P = 1 / 200, AlphaP = 2)
#' }

#'

ZIFA_LSNM <- function(X, num_fac = 2, NU1P = 1, NU2P = 3, G1P = 200, G2P = 1 / 200, AlphaP = 2) {
  ###############  Count Matrix  ###############
  x <- X

  ###############  Number of Factors  ###############
  k <- num_fac

  ###############  Number of samples  ###############
  n <- nrow(x)

  ###############  Number of Taxas  ###############
  p <- ncol(x)

  ###############  Some Extra Information ###############
  M <- rowSums(x)

  # shape for kappa
  nu1 <- NU1P

  # Shape for kappa
  nu2 <- NU2P

  # Shape for Gamma
  g1 <- G1P

  # Rate for Gamma
  g2 <- G2P

  # Skewness for Latent Factor
  prior_alpha <- AlphaP

  #--------------------------------------------------------------------------#
  ############### Initial Guess For The Variational Parameters ###############
  #--------------------------------------------------------------------------#

  X_transformation <- scale(log(x + 0.05), scale = T, center = F)
  Decomposition_X <- svd(X_transformation, k, k)

  r <- updated_r <- Decomposition_X$v

  lambda <- updated_lambda <- matrix(0.65, p, k)

  g_1 <- updated_g_1 <- matrix(2, p, k)

  g_2 <- updated_g_2 <- matrix(2, p, k)

  tau_1 <- updated_tau_1 <- rep(2, p)

  tau_2 <- updated_tau_2 <- rep(2, p)

  xi <- updated_xi <- Decomposition_X$u %*% diag(Decomposition_X$d[1:k])

  omega <- updated_omega <- matrix(0.45, n, k)

  alpha <- updated_alpha <- matrix(3, n, k)

  zero_ratio <- apply(x, 2, function(y) {
    sum(y == 0) / n
  })
  pi_ij <- updated_pi_ij <- t(ifelse(t(x) == 0, zero_ratio, 0))
  pi_matrix <- matrix(0, n, p)

  a_0j <- updated_a_0j <- runif(p)

  c_0j <- updated_c_0j <- rep(0.35, p)

  #----------------------------------------------------------#
  ###############  CAVB For ZIFA - LSNM Model  ###############
  #----------------------------------------------------------#

  old_ELBO <- -1e8
  max_iteration <- 500
  difference <- 1e5
  relative_difference <- 1e5
  iter <- 1
  eps <- 5e-5
  delta <- 1e-10
  meanskewnormal_before <- matrix(2, n, k, TRUE)
  r_trace <- meanskewnormal_trace <- elbo_trace <- list()

  for (iter in 1:max_iteration) {
    #------------------------------------------------#
    ###############  Update For Pi_ij  ###############
    #------------------------------------------------#

    TomAndJerry <- matrix((updated_a_0j + 0.5 * updated_c_0j), n, p, byrow = TRUE) + L_ij(new.R = updated_r, new.LAMBDA = updated_lambda, new.XI = updated_xi, new.OMEGA = updated_omega, new.ALPHA = updated_alpha, NROW = n, PCOL = p)

    upsilon <- log(M / (rowSums((1 - updated_pi_ij) * exp(TomAndJerry))))
    for (i in 1:n) {
      for (j in 1:p) {
        pi_matrix[i, j] <- exp(digamma(updated_tau_1[j])) / (exp(digamma(updated_tau_1[j])) + exp(digamma(updated_tau_2[j]) - exp(upsilon[i] + updated_a_0j[j] + 0.5 * updated_c_0j[j] + TomAndJerry[i, j])))
      }
    }
    pi_matrix[x != 0] <- 0
    updated_pi_ij <- ifelse(pi_matrix > 0.5, 1, 0)
    cat("PIIJ updated on iteration step : ", iter, "\n")

    #-------------------------------------------------------#
    ###############  Optimization For Tau_j1  ###############
    #-------------------------------------------------------#

    for (j in 1:p) {
      updated_tau_1[j] <- optim(tau_1[j], UPDATED_TAU_2 = updated_tau_2, PIIJ = updated_pi_ij, NU_1 = nu1, NU_2 = nu2, l = j, N = n, tau_j1_function, tau_j1_gradient, method = "L-BFGS-B", lower = 0.1, upper = +Inf, control = list(fnscale = -1))$par
    }
    cat("Tau_j1 updated on iteration step : ", iter, "\n")

    #-------------------------------------------------------#
    ###############  Optimization For Tau_j2  ###############
    #-------------------------------------------------------#

    for (j in 1:p) {
      updated_tau_2[j] <- optim(tau_2[j], UPDATED_TAU_1 = updated_tau_1, PIIJ = updated_pi_ij, NU_1 = nu1, NU_2 = nu2, l = j, N = n, tau_j2_function, tau_j2_gradient, method = "L-BFGS-B", lower = 0.1, upper = +Inf, control = list(fnscale = -1))$par
    }
    cat("Tau_j2 updated on iteration step : ", iter, "\n")

    #-------------------------------------------------------#
    ###############  Optimization For G1 & G2  ##############
    #-------------------------------------------------------#

    q21 <- try(optim(c(g_1), G2 = updated_g_2, R = updated_r, LAMBDA = updated_lambda, G1_Prior = g1, G2_Prior = g2, N = n, P = p, K = k, method = "L-BFGS-B", lower = matrix(1e-2, p, k), upper = matrix(+Inf, p, k), fn = g1_function, gr = g1_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)
    if ("try-error" %in% class(q21)) {
      updated_g_1 <- g_1

      cat("G1 cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_g_1 <- q21$par

      cat("G1 updated on iteration step : ", iter, "\n")
    }

    q22 <- try(optim(c(g_2), G1 = updated_g_1, R = updated_r, LAMBDA = updated_lambda, G1_Prior = g1, G2_Prior = g2, N = n, P = p, K = k, method = "L-BFGS-B", lower = matrix(1e-2, p, k), upper = matrix(+Inf, p, k), fn = g2_function, gr = g2_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)
    if ("try-error" %in% class(q22)) {
      updated_g_2 <- g_2

      cat("G1 cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_g_2 <- q22$par

      cat("G1 updated on iteration step : ", iter, "\n")
    }

    #--------------------------------------------------#
    ###############  Optimization For R  ###############
    #--------------------------------------------------#

    q1 <- try(optim(c(r), COUNT_MATRIX = x, LAMBDA = updated_lambda, XI = updated_xi, OMEGA = updated_omega, ALPHA = updated_alpha, G1 = updated_g_1, G2 = updated_g_2, A0 = updated_a_0j, C0 = updated_c_0j, PIIJ = updated_pi_ij, N = n, P = p, K = k, M = M, method = "BFGS", fn = r_function, gr = r_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)

    if ("try-error" %in% class(q1)) {
      updated_r <- r

      cat("R cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_r <- matrix(q1$par, p, k)

      cat("R updated on iteration step : ", iter, "\n")
    }

    r_trace[[iter]] <- norm(matrix(q1$par, p, k, TRUE) - updated_r, type = "F") / (norm(updated_r, type = "F") + 1e-8)


    #-------------------------------------------------------#
    ###############  Optimization For Lambda  ###############
    #-------------------------------------------------------#

    q2 <- try(optim(c(lambda), R = updated_r, XI = updated_xi, OMEGA = updated_omega, ALPHA = updated_alpha, G1 = updated_g_1, G2 = updated_g_2, A0 = updated_a_0j, C0 = updated_c_0j, PIIJ = updated_pi_ij, N = n, P = p, K = k, M = M, method = "L-BFGS-B", lower = matrix(1e-2, p, k), upper = matrix(1, p, k), fn = lambda_function, gr = lambda_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)
    if ("try-error" %in% class(q2)) {
      updated_lambda <- lambda

      cat("Lambda cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_lambda <- matrix(q2$par, p, k)

      cat("Lambda updated on iteration step : ", iter, "\n")
    }

    #---------------------------------------------------#
    ###############  Optimization For Xi  ###############
    #---------------------------------------------------#

    q3 <- try(optim(c(xi), COUNT_MATRIX = x, R = updated_r, LAMBDA = updated_lambda, OMEGA = updated_omega, ALPHA = updated_alpha, A0 = updated_a_0j, C0 = updated_c_0j, PIIJ = updated_pi_ij, ALPHA_Prior = prior_alpha, N = n, P = p, K = k, M = M, method = "BFGS", fn = xi_function, gr = xi_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)

    if ("try-error" %in% class(q3)) {
      updated_xi <- xi

      cat("Xi cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_xi <- matrix(q3$par, n, k)

      cat("Xi updated on iteration step : ", iter, "\n")
    }

    #------------------------------------------------------#
    ###############  Optimization For Omega  ###############
    #------------------------------------------------------#

    q4 <- try(optim(c(omega), COUNT_MATRIX = x, R = updated_r, LAMBDA = updated_lambda, XI = updated_xi, ALPHA = updated_alpha, A0 = updated_a_0j, C0 = updated_c_0j, PIIJ = updated_pi_ij, ALPHA_Prior = prior_alpha, N = n, P = p, K = k, M = M, method = "L-BFGS-B", lower = 1e-5, upper = 1, fn = omega_function, gr = omega_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)

    if ("try-error" %in% class(q4)) {
      updated_omega <- omega

      cat("Omega cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_omega <- matrix(q4$par, n, k)

      cat("Omega updated on iteration step : ", iter, "\n")
    }

    #------------------------------------------------------#
    ###############  Optimization For Alpha  ###############
    #------------------------------------------------------#

    q5 <- try(optim(c(alpha), COUNT_MATRIX = x, R = updated_r, LAMBDA = updated_lambda, XI = updated_xi, OMEGA = updated_omega, A0 = updated_a_0j, C0 = updated_c_0j, PIIJ = updated_pi_ij, ALPHA_Prior = prior_alpha, N = n, P = p, K = k, M = M, method = "BFGS", fn = alpha_function, gr = alpha_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)

    if ("try-error" %in% class(q5)) {
      updated_alpha <- alpha

      cat("Alpha cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_alpha <- matrix(q5$par, n, k)

      cat("Alpha updated on iteration step : ", iter, "\n")
    }

    meanskewnormal_after <- updated_xi + updated_omega * (updated_alpha / sqrt(1 + updated_alpha^2)) * sqrt(2 / pi)

    meanskewnormal_trace[[iter]] <- norm(matrix(meanskewnormal_after, n, k, TRUE) - meanskewnormal_before, type = "F") / (norm(meanskewnormal_before, type = "F") + 1e-8)

    meanskewnormal_before <- meanskewnormal_after

    #------------------------------------------------------#
    ###############  Optimization For A0  ##################
    #------------------------------------------------------#

    q6 <- try(optim(c(a_0j), COUNT_MATRIX = x, R = updated_r, LAMBDA = updated_lambda, XI = updated_xi, OMEGA = updated_omega, ALPHA = updated_alpha, PIIJ = updated_pi_ij, C0 = updated_c_0j, N = n, P = p, K = k, M = M, method = "BFGS", fn = a0_function, gr = a0_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)
    if ("try-error" %in% class(q6)) {
      updated_a_0j <- a_0j

      cat("A0 cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_a_0j <- q6$par

      cat("A0 updated on iteration step : ", iter, "\n")
    }

    #------------------------------------------------------#
    ###############  Optimization For C0  ##################
    #------------------------------------------------------#

    q7 <- try(optim(c(c_0j), R = updated_r, LAMBDA = updated_lambda, XI = updated_xi, OMEGA = updated_omega, ALPHA = updated_alpha, PIIJ = updated_pi_ij, A0 = updated_a_0j, N = n, P = p, K = k, M = M, method = "L-BFGS-B", lower = 1e-5, upper = +Inf, fn = c0_function, gr = c0_gradient, control = list(fnscale = -1, maxit = 1000)), silent = TRUE)
    if ("try-error" %in% class(q7)) {
      updated_c_0j <- c_0j

      cat("C0 cannot be updated on iteration step : ", iter, "\n")
    } else {
      updated_c_0j <- q7$par

      cat("C0 updated on iteration step : ", iter, "\n")
    }

    ###############################################
    ###############  Stopping Rule  ###############
    ###############################################


    q <- list(value = ELBO_AP(COUNT_MATRIX = x, R = updated_r, LAMBDA = updated_lambda, XI = updated_xi, OMEGA = updated_omega, ALPHA = updated_alpha, PIIJ = updated_pi_ij, A0 = updated_a_0j, C0 = updated_c_0j, T1 = updated_tau_1, T2 = updated_tau_2, G1 = updated_g_1, G2 = updated_g_2, G1_Prior = g1, G2_Prior = g2, NU1 = nu1, NU2 = nu2, ALPHA_Prior = prior_alpha, N = n, P = p, K = k, M = M))
    new_ELBO <- q$value
    elbo_trace[[iter]] <- q$value
    difference <- abs(new_ELBO - old_ELBO)
    ratio <- abs(new_ELBO / old_ELBO)
    relative_difference <- abs(new_ELBO - old_ELBO) / abs(old_ELBO + delta)
    if (iter == 1) {
      relative_difference_r_trace <- 1
    } else {
      relative_difference_r_trace <- abs(r_trace[[iter]] - r_trace[[iter - 1]]) / abs(r_trace[[iter]] + delta)
    }
    cat("Old ELBO : ", old_ELBO, ",New ELBO : ", new_ELBO, ",Ratio of ELBO's :", ratio, ",Difference of ELBO's :", difference, ",Relative Difference :", relative_difference, "\n", "Relative Difference R trace :", relative_difference_r_trace, "\n")

    old_ELBO <- new_ELBO

    tau_1 <- updated_tau_1
    tau_2 <- updated_tau_2
    g_1 <- updated_g_1
    g_2 <- updated_g_2
    r <- updated_r
    lambda <- updated_lambda
    xi <- updated_xi
    omega <- updated_omega
    alpha <- updated_alpha
    a_0j <- updated_a_0j
    c_0j <- updated_c_0j
    iter <- iter + 1

    if (relative_difference < eps & relative_difference_r_trace < 1e-4) {
      break
    }
    if (iter == max_iteration) {
      cat("Max Iteration Achieved")
      break
    }
  }

  ##################################################
  ###############  Estimates of Rho  ###############
  ##################################################

  DonalDuck <- matrix((updated_a_0j + 0.5 * updated_c_0j), n, p, byrow = TRUE) + L_ij(updated_r, updated_lambda, updated_xi, updated_omega, updated_alpha, n, p)
  estimated_compositions <- exp(DonalDuck) / (rowSums(exp(DonalDuck)) + 1e-8)

  ######################################################
  ###############  Returning the Output  ###############
  ######################################################

  return_list <- list()
  return_list$Pi_ij <- matrix(pi_matrix, n, p)
  return_list$Tau_1 <- updated_tau_1
  return_list$Tau_2 <- updated_tau_2
  return_list$R <- matrix(updated_r, p, k)
  return_list$R_Trace <- r_trace
  return_list$Lambda <- matrix(updated_lambda, p, k)
  return_list$G1 <- updated_g_1
  return_list$G2 <- updated_g_2
  return_list$Xi <- matrix(updated_xi, n, k)
  return_list$Omega <- matrix(updated_omega, n, k)
  return_list$Alpha <- matrix(updated_alpha, n, k)
  return_list$MeanSkewNormal_Trace <- meanskewnormal_trace
  return_list$A0 <- updated_a_0j
  return_list$C0 <- updated_c_0j
  return_list$Estimated_Compositions <- estimated_compositions
  return_list$ELBO_Trace <- elbo_trace
  return(return_list)
}
