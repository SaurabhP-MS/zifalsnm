#---------------------------------------------------#
###############  Function For Tau_j1  ###############
#---------------------------------------------------#

tau_j1_function <- function(tau_1, UPDATED_TAU_2, PIIJ, NU_1, NU_2, l, N) {
  tau_functionj1 <- sapply(1:N, function(i) (PIIJ[i, l] * digamma(tau_1)))
  t1 <- sum(tau_functionj1) + (digamma(tau_1 + UPDATED_TAU_2[l]) * (-N + tau_1 + UPDATED_TAU_2[l] - NU_1 - NU_2)) + (digamma(tau_1) * (NU_1 - tau_1)) + lbeta(tau_1, UPDATED_TAU_2[l])
  return(t1)
}
tau_j1_gradient <- function(tau_1, UPDATED_TAU_2, PIIJ, NU_1, NU_2, l, N) {
  tau_gradientj1 <- sapply(1:N, function(i) {
    (PIIJ[i, l] * trigamma(tau_1))
  })
  grad1 <- sum(tau_gradientj1) + (trigamma(tau_1 + UPDATED_TAU_2[l]) * (-N + tau_1 + UPDATED_TAU_2[l] - NU_1 - NU_2)) + (trigamma(tau_1) * (NU_1 - tau_1))
  return(grad1)
}

#---------------------------------------------------#
###############  Function For Tau_j2  ###############
#---------------------------------------------------#

tau_j2_function <- function(tau_2, UPDATED_TAU_1, PIIJ, NU_1, NU_2, l, N) {
  tau_functionj2 <- sapply(1:N, function(i) {
    -PIIJ[i, l] * digamma(tau_2)
  })
  t2 <- sum(tau_functionj2) + (digamma(UPDATED_TAU_1[l] + tau_2) * (-N + UPDATED_TAU_1[l] + tau_2 - NU_1 - NU_2)) + (digamma(tau_2) * (N + NU_2 - tau_2)) + lbeta(UPDATED_TAU_1[l], tau_2)
  return(t2)
}
tau_j2_gradient <- function(tau_2, UPDATED_TAU_1, PIIJ, NU_1, NU_2, l, N) {
  tau_gradientj2 <- sapply(1:N, function(i) {
    -PIIJ[i, l] * trigamma(tau_2)
  })
  grad2 <- sum(tau_gradientj2) + (trigamma(tau_2) * (N + NU_2 - tau_2)) + (trigamma(UPDATED_TAU_1[l] + tau_2) * (-N - NU_1 - NU_2 + UPDATED_TAU_1[l] + tau_2))
  return(grad2)
}

#---------------------------------------------------#
###############  Function For G1  ###################
#---------------------------------------------------#

g1_function <- function(G1, G2, R, LAMBDA, G1_Prior, G2_Prior, N, P, K) {
  g1_matrix <- matrix(G1, P, K)
  g1_func1 <- (G1_Prior - 0.5) * digamma(g1_matrix) - (g1_matrix / (G2 + 1e-8)) * (0.5 * LAMBDA + 0.5 * R^2 + G2_Prior) + g1_matrix + lgamma(g1_matrix) + (1 - g1_matrix) * digamma(g1_matrix)
  return(sum(g1_func1))
}

g1_gradient <- function(G1, G2, R, LAMBDA, G1_Prior, G2_Prior, N, P, K) {
  g1_matrix <- matrix(G1, P, K)
  g1_grad1 <- (G1_Prior - 0.5) * trigamma(g1_matrix) - (1 / (G2 + 1e-8)) * (0.5 * LAMBDA + 0.5 * R^2 + G2_Prior) + 1 + (1 - g1_matrix) * trigamma(g1_matrix)
  return(c(g1_grad1))
}

#---------------------------------------------------#
###############  Function For G2  ###################
#---------------------------------------------------#

g2_function <- function(G1, G2, R, LAMBDA, G1_Prior, G2_Prior, N, P, K) {
  g2_matrix <- matrix(G2, P, K)
  g2_func1 <- (G1_Prior - 0.5) * log(g2_matrix + 1e-8) + (G1 / (g2_matrix + 1e-8)) * (0.5 * LAMBDA + 0.5 * R^2 + G2_Prior) + log(g2_matrix + 1e-8)
  return(-sum(g2_func1))
}

g2_gradient <- function(G1, G2, R, LAMBDA, G1_Prior, G2_Prior, N, P, K) {
  g2_matrix <- matrix(G2, P, K)
  g2_grad1 <- (1 / (g2_matrix + 1e-8)) * (G1_Prior + 0.5) - (G1 / (g2_matrix^2)) * (0.5 * LAMBDA + 0.5 * R^2 + G2_Prior)
  return(-c(g2_grad1))
}

#----------------------------------------------#
###############  Function For R  ###############
#----------------------------------------------#

r_function <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, G1, G2, A0, C0, PIIJ, N, P, K, M) {
  R_matrix <- matrix(R, P, K)
  xi_omega_alpha_matrix <- XI + (sqrt(2 / pi) * OMEGA * ALPHA / (sqrt(1 + ALPHA^2) + 1e-8))
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R_matrix, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)

  r_func1 <- COUNT_MATRIX * (tcrossprod(xi_omega_alpha_matrix, R_matrix))
  r_func3 <- (R_matrix^2) * (G1 / (G2 + 1e-8))
  r_func4 <- M * log(rowSums((1 - PIIJ) * exp(TomJerry)))

  r1 <- sum(r_func1) - 0.5 * sum(r_func3) - sum(r_func4)
  return(r1)
}

r_gradient <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, G1, G2, A0, C0, PIIJ, N, P, K, M) {
  R_matrix <- matrix(R, P, K)
  r_grad3 <- matrix(0, P, K)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R_matrix, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)
  xi_omega_alpha_matrix <- XI + sqrt(2 / pi) * (OMEGA * ALPHA / (sqrt(1 + ALPHA^2) + 1e-8))
  R_grad_matrix <- R_ij(new.R = R_matrix, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P, KFACT = K)

  r_grad1 <- crossprod(COUNT_MATRIX, xi_omega_alpha_matrix) - R_matrix * (G1 / (G2 + 1e-8))
  r_grad2 <- ((1 - PIIJ) * exp(TomJerry)) / (rowSums((1 - PIIJ) * (exp(TomJerry))) + 1e-8)
  for (t in 1:K) {
    r_grad3[, t] <- colSums(M * r_grad2 * R_grad_matrix[, , t])
  }
  r1 <- r_grad1 - r_grad3
  return(c(r1))
}

#---------------------------------------------------#
###############  Function For Lambda  ###############
#---------------------------------------------------#

lambda_function <- function(R, LAMBDA, XI, OMEGA, ALPHA, G1, G2, A0, C0, PIIJ, N, P, K, M) {
  LAMBDA_Matrix <- matrix(LAMBDA, P, K)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA_Matrix, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)
  l_func1 <- M * log(rowSums((1 - PIIJ) * exp(TomJerry)))

  l_func2 <- (LAMBDA_Matrix) * (G1 / (G2 + 1e-8)) - log(LAMBDA_Matrix + 1e-8)
  l1 <- -sum(l_func1) - 0.5 * sum(l_func2)
  return(l1)
}

lambda_gradient <- function(R, LAMBDA, XI, OMEGA, ALPHA, G1, G2, A0, C0, PIIJ, N, P, K, M) {
  LAMBDA_Matrix <- matrix(LAMBDA, P, K)
  l_grad3 <- matrix(0, P, K)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA_Matrix, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)
  lambda_grad_matrix <- M_ij(new.R = R, new.LAMBDA = LAMBDA_Matrix, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P, KFACT = K)

  l_grad1 <- 0.5 * ((1 / (LAMBDA_Matrix + 1e-8)) - (G1 / (G2 + 1e-8)))
  l_grad2 <- ((1 - PIIJ) * exp(TomJerry)) / (rowSums((1 - PIIJ) * (exp(TomJerry))) + 1e-8)
  for (t in 1:K) {
    l_grad3[, t] <- colSums(M * l_grad2 * lambda_grad_matrix[, , t])
  }

  l1 <- l_grad1 - l_grad3
  return(c(l1))
}

#-----------------------------------------------#
###############  Function For Xi  ###############
#-----------------------------------------------#

xi_function <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, A0, C0, PIIJ, ALPHA_Prior, N, P, K, M) {
  XI_Matrix <- matrix(XI, N, K)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI_Matrix, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)

  xi_func1 <- COUNT_MATRIX * (tcrossprod(XI_Matrix, R))
  xi_func2 <- M * log(rowSums((1 - PIIJ) * exp(TomJerry)))
  xi_func3 <- (0.5) * XI_Matrix^2 + sqrt(2 / pi) * (XI_Matrix * OMEGA * ALPHA / (sqrt(1 + ALPHA^2) + 1e-8)) - (-sqrt(2 / pi) * ALPHA_Prior / (sqrt(1 + ALPHA_Prior^2))) * XI_Matrix - psi_1(XI_Matrix, OMEGA, ALPHA, ALPHA_Prior)

  x1 <- sum(xi_func1) - sum(xi_func2) - sum(xi_func3)
  return(x1)
}

xi_gradient <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, A0, C0, PIIJ, ALPHA_Prior, N, P, K, M) {
  XI_Matrix <- matrix(XI, N, K)
  xi_grad3 <- matrix(0, N, K)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI_Matrix, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)
  xi_grad_matrix <- E_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI_Matrix, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P, KFACT = K)

  xi_grad1 <- COUNT_MATRIX %*% R - XI - (sqrt(2 / pi) * (OMEGA * ALPHA / (sqrt(1 + ALPHA^2) + 1e-8))) + (-sqrt(2 / pi) * ALPHA_Prior / (sqrt(1 + ALPHA_Prior^2))) + psi_1_derivative_xi(XI_Matrix, OMEGA, ALPHA, ALPHA_Prior)
  xi_grad2 <- ((1 - PIIJ) * exp(TomJerry)) / (rowSums((1 - PIIJ) * (exp(TomJerry))) + 1e-8)
  for (t in 1:K) {
    xi_grad3[, t] <- rowSums(M * xi_grad2 * xi_grad_matrix[, , t])
  }

  x1 <- xi_grad1 - xi_grad3
  return(c(x1))
}

#--------------------------------------------------#
###############  Function For Omega  ###############
#--------------------------------------------------#

omega_function <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, A0, C0, PIIJ, ALPHA_Prior, N, P, K, M) {
  OMEGA_Matrix <- matrix(OMEGA, N, K)
  omega_alpha_matrix <- OMEGA_Matrix * ALPHA / (sqrt(1 + ALPHA^2) + 1e-8)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA_Matrix, new.ALPHA = ALPHA, NROW = N, PCOL = P)
  alpha_star <- -sqrt(2 / pi) * (ALPHA_Prior / sqrt(1 + ALPHA_Prior^2))

  o_func1 <- sqrt(2 / pi) * COUNT_MATRIX * (tcrossprod(omega_alpha_matrix, R))
  o_func2 <- M * log(rowSums((1 - PIIJ) * exp(TomJerry)))
  o_func3 <- (0.5) * OMEGA_Matrix^2 + sqrt(2 / pi) * XI * omega_alpha_matrix - sqrt(2 / pi) * alpha_star * omega_alpha_matrix - log(OMEGA_Matrix + 1e-8) - psi_1(XI, OMEGA_Matrix, ALPHA, ALPHA_Prior)

  o1 <- sum(o_func1) - sum(o_func2) - sum(o_func3)
  return(o1)
}

omega_gradient <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, A0, C0, PIIJ, ALPHA_Prior, N, P, K, M) {
  OMEGA_Matrix <- matrix(OMEGA, N, K)
  o_grad3 <- matrix(0, N, K)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA_Matrix, new.ALPHA = ALPHA, NROW = N, PCOL = P)
  omega_grad_matrix <- O_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA_Matrix, new.ALPHA = ALPHA, NROW = N, PCOL = P, KFACT = K)
  alpha_star <- -sqrt(2 / pi) * (ALPHA_Prior / sqrt(1 + ALPHA_Prior^2))

  o_grad1 <- sqrt(2 / pi) * (COUNT_MATRIX %*% R) * (ALPHA / (sqrt(1 + ALPHA^2) + 1e-8)) - OMEGA_Matrix - sqrt(2 / pi) * XI * ALPHA / (sqrt(1 + ALPHA^2) + 1e-8) + sqrt(2 / pi) * alpha_star * ALPHA / (sqrt(1 + ALPHA^2) + 1e-8) + (1 / (OMEGA_Matrix + 1e-8)) + psi_1_derivative_omega(XI, OMEGA_Matrix, ALPHA, ALPHA_Prior)
  o_grad2 <- ((1 - PIIJ) * exp(TomJerry)) / (rowSums((1 - PIIJ) * (exp(TomJerry))) + 1e-8)
  for (t in 1:K) {
    o_grad3[, t] <- rowSums(M * o_grad2 * omega_grad_matrix[, , t])
  }

  o1 <- o_grad1 - o_grad3
  return(c(o1))
}

#--------------------------------------------------#
###############  Function For Alpha  ###############
#--------------------------------------------------#

alpha_function <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, A0, C0, PIIJ, ALPHA_Prior, N, P, K, M) {
  ALPHA_Matrix <- matrix(ALPHA, N, K)
  omega_alpha_matrix <- OMEGA * ALPHA_Matrix / (sqrt(1 + ALPHA_Matrix^2) + 1e-8)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA_Matrix, NROW = N, PCOL = P)
  alpha_star <- -sqrt(2 / pi) * (ALPHA_Prior / sqrt(1 + ALPHA_Prior^2))

  a_func1 <- sqrt(2 / pi) * COUNT_MATRIX * (tcrossprod(omega_alpha_matrix, R))
  a_func2 <- M * log(rowSums((1 - PIIJ) * exp(TomJerry)))
  a_func3 <- sqrt(2 / pi) * XI * omega_alpha_matrix - sqrt(2 / pi) * alpha_star * omega_alpha_matrix - psi_1(XI, OMEGA, ALPHA_Matrix, ALPHA_Prior) + psi_2(ALPHA_Matrix)

  a1 <- sum(a_func1) - sum(a_func2) - sum(a_func3)
  return(a1)
}

alpha_gradient <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, A0, C0, PIIJ, ALPHA_Prior, N, P, K, M) {
  ALPHA_Matrix <- matrix(ALPHA, N, K)
  TomJerry <- matrix((A0 + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA_Matrix, NROW = N, PCOL = P)
  a_grad3 <- matrix(0, N, K)
  alpha_grad_matrix <- A_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA_Matrix, NROW = N, PCOL = P, KFACT = K)
  alpha_star <- -sqrt(2 / pi) * (ALPHA_Prior / sqrt(1 + ALPHA_Prior^2))

  a_grad1 <- sqrt(2 / pi) * OMEGA * (COUNT_MATRIX %*% R) / ((1 + ALPHA_Matrix^2)^1.5 + 1e-8) - sqrt(2 / pi) * OMEGA * XI / ((1 + ALPHA_Matrix^2)^1.5 + 1e-8) + sqrt(2 / pi) * alpha_star * OMEGA / ((1 + ALPHA_Matrix^2)^1.5 + 1e-8) + psi_1_derivative_alpha(XI, OMEGA, ALPHA_Matrix, ALPHA_Prior) - psi_2_derivative(ALPHA_Matrix)
  a_grad2 <- ((1 - PIIJ) * exp(TomJerry)) / (rowSums((1 - PIIJ) * (exp(TomJerry))) + 1e-8)
  for (t in 1:K) {
    a_grad3[, t] <- rowSums(M * a_grad2 * alpha_grad_matrix[, , t])
  }

  a1 <- a_grad1 - a_grad3
  return(c(a1))
}

#--------------------------------------------------#
###############  Function For A0  ##################
#--------------------------------------------------#

a0_function <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, PIIJ, A0, C0, N, P, K, M) {
  a0_vector <- matrix(A0, P, 1)
  TomJerry <- matrix((a0_vector + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)

  a0_func1 <- COUNT_MATRIX %*% a0_vector
  a0_func2 <- M * log(rowSums((1 - PIIJ) * exp(TomJerry)))
  a0_func3 <- 0.5 * a0_vector^2

  a1 <- sum(a0_func1) - sum(a0_func2) - sum(a0_func3)
  sum(a1)
}

a0_gradient <- function(COUNT_MATRIX, R, LAMBDA, XI, OMEGA, ALPHA, PIIJ, A0, C0, N, P, K, M) {
  a0_vector <- matrix(A0, P, 1)
  TomJerry <- matrix((a0_vector + 0.5 * C0), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)
  a0_grad <- COUNT_MATRIX - (M * ((1 - PIIJ) * exp(TomJerry)) / (rowSums((1 - PIIJ) * (exp(TomJerry))) + 1e-8))

  return(c(colSums(a0_grad)) - a0_vector)
}

#--------------------------------------------------#
###############  Function For C0  ##################
#--------------------------------------------------#

c0_function <- function(R, LAMBDA, XI, OMEGA, ALPHA, PIIJ, A0, C0, N, P, K, M) {
  c0_vector <- matrix(C0, P, 1)
  TomJerry <- matrix((A0 + 0.5 * c0_vector), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)

  c0_func1 <- M * log(rowSums((1 - PIIJ) * exp(TomJerry)))
  c0_func2 <- 0.5 * (c0_vector - log(c0_vector + 1e-8))

  c1 <- -sum(c0_func1) - sum(c0_func2)
  sum(c1)
}

c0_gradient <- function(R, LAMBDA, XI, OMEGA, ALPHA, PIIJ, A0, C0, N, P, K, M) {
  c0_vector <- matrix(C0, P, 1)
  TomJerry <- matrix((A0 + 0.5 * c0_vector), N, P, byrow = TRUE) + L_ij(new.R = R, new.LAMBDA = LAMBDA, new.XI = XI, new.OMEGA = OMEGA, new.ALPHA = ALPHA, NROW = N, PCOL = P)

  c0_grad <- -0.5 * (M * ((1 - PIIJ) * exp(TomJerry)) / (rowSums((1 - PIIJ) * (exp(TomJerry))) + 1e-8))

  return(c(colSums(c0_grad)) - 0.5 * (1 - 1 / (c0_vector + 1e-8)))
}
