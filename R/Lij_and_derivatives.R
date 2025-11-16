#---------------------------------------------#
###############  L_ij Function  ###############
#---------------------------------------------#

L_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL) {
  one_minus_lambda_omega_sq <- 1 - new.LAMBDA[rep(1:PCOL, each = NROW), ] * new.OMEGA[rep(1:NROW, PCOL), ]^2
  sqrt_term <- sqrt(one_minus_lambda_omega_sq)

  numerator_1 <- new.LAMBDA[rep(1:PCOL, each = NROW), ] * new.XI[rep(1:NROW, PCOL), ]^2 +
    new.R[rep(1:PCOL, each = NROW), ]^2 * new.OMEGA[rep(1:NROW, PCOL), ]^2 +
    2 * new.R[rep(1:PCOL, each = NROW), ] * new.XI[rep(1:NROW, PCOL), ]

  numerator_2 <- (new.R[rep(1:PCOL, each = NROW), ] + new.LAMBDA[rep(1:PCOL, each = NROW), ] * new.XI[rep(1:NROW, PCOL), ]) *
    new.OMEGA[rep(1:NROW, PCOL), ] * new.ALPHA[rep(1:NROW, PCOL), ]

  denominator_2 <- sqrt_term * sqrt(1 + new.ALPHA[rep(1:NROW, PCOL), ]^2 - new.LAMBDA[rep(1:PCOL, each = NROW), ] * new.OMEGA[rep(1:NROW, PCOL), ]^2) + 1e-8

  log_term <- log(2) - 0.5 * log(one_minus_lambda_omega_sq + 1e-8) + (numerator_1 / (2 * one_minus_lambda_omega_sq + 1e-8))
  pnorm_term <- pnorm(numerator_2 / denominator_2, log.p = TRUE)

  log_term_matrix <- matrix(rowSums(matrix(log_term, nrow = NROW * PCOL)), nrow = NROW, byrow = FALSE)
  pnorm_matrix <- matrix(rowSums(matrix(pnorm_term, nrow = NROW * PCOL)), nrow = NROW, byrow = FALSE)

  result_matrix <- log_term_matrix + pnorm_matrix

  return(result_matrix)
}

#-------------------------------------------------#
###############  A Common Function  ###############
#-------------------------------------------------#

phi_by_Phi <- function(R_value, LAMBDA_value, XI_value, OMEGA_value, ALPHA_value) {
  numerator <- (R_value + LAMBDA_value * XI_value) * OMEGA_value * ALPHA_value
  denominator <- sqrt((1 - LAMBDA_value * OMEGA_value^2) * (1 + ALPHA_value^2 - LAMBDA_value * OMEGA_value^2)) + 1e-8
  ratio <- numerator / denominator
  dnorm(ratio) / pnorm(ratio)
}

#----------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to R  ###############
#----------------------------------------------------------------------#

R_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  # Create 3D arrays by recycling matrices
  R_array <- array(new.R, dim = c(PCOL, KFACT, NROW))
  LAMBDA_array <- array(new.LAMBDA, dim = c(PCOL, KFACT, NROW))
  XI_array <- array(aperm(new.XI, c(2, 1)), dim = c(1, KFACT, NROW))
  OMEGA_array <- array(aperm(new.OMEGA, c(2, 1)), dim = c(1, KFACT, NROW))
  ALPHA_array <- array(aperm(new.ALPHA, c(2, 1)), dim = c(1, KFACT, NROW))

  # Align dimensions for vectorized computations
  XI_array <- XI_array[rep(1, PCOL), , ]
  OMEGA_array <- OMEGA_array[rep(1, PCOL), , ]
  ALPHA_array <- ALPHA_array[rep(1, PCOL), , ]

  numerator1 <- R_array * OMEGA_array^2 + XI_array
  denominator1 <- (1 - LAMBDA_array * OMEGA_array^2) + 1e-8
  term1 <- numerator1 / denominator1

  numerator2 <- OMEGA_array * ALPHA_array
  denominator2 <- sqrt(1 - LAMBDA_array * OMEGA_array^2) * sqrt(1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2) + 1e-8
  term2 <- phi_by_Phi(R_array, LAMBDA_array, XI_array, OMEGA_array, ALPHA_array) * (numerator2 / denominator2)

  R_result <- term1 + term2

  # Rearrange dimensions back to NROW x PCOL x KFACT
  return(aperm(R_result, c(3, 1, 2)))
}

#---------------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Lambda  ###############
#---------------------------------------------------------------------------#

M_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  # Create 3D arrays by recycling matrices
  R_array <- array(new.R, dim = c(PCOL, KFACT, NROW))
  LAMBDA_array <- array(new.LAMBDA, dim = c(PCOL, KFACT, NROW))
  XI_array <- array(aperm(new.XI, c(2, 1)), dim = c(1, KFACT, NROW))
  OMEGA_array <- array(aperm(new.OMEGA, c(2, 1)), dim = c(1, KFACT, NROW))
  ALPHA_array <- array(aperm(new.ALPHA, c(2, 1)), dim = c(1, KFACT, NROW))

  # Align dimensions for vectorized computations
  XI_array <- XI_array[rep(1, PCOL), , ]
  OMEGA_array <- OMEGA_array[rep(1, PCOL), , ]
  ALPHA_array <- ALPHA_array[rep(1, PCOL), , ]

  denominator1 <- (1 - LAMBDA_array * OMEGA_array^2)
  term1 <- 0.5 * (OMEGA_array^2 / (denominator1 + 1e-8))
  term2 <- 0.5 * ((XI_array^2 + R_array^2 * OMEGA_array^4 + 2 * R_array * XI_array * OMEGA_array^2) / (denominator1^2 + 1e-8))

  denominator2 <- sqrt(1 - LAMBDA_array * OMEGA_array^2) * sqrt(1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2)
  denominator3 <- 2 * ((1 - LAMBDA_array * OMEGA_array^2) * (1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2))^1.5
  term3 <- (XI_array * OMEGA_array * ALPHA_array) / (denominator2 + 1e-8)
  term4 <- ((R_array + LAMBDA_array * XI_array) * OMEGA_array * ALPHA_array * (-2 * OMEGA_array^2 - OMEGA_array^2 * ALPHA_array^2 + 2 * LAMBDA_array * OMEGA_array^4)) / (denominator3 + 1e-8)

  term5 <- phi_by_Phi(R_array, LAMBDA_array, XI_array, OMEGA_array, ALPHA_array) * (term3 - term4)

  M_result <- term1 + term2 + term5

  # Rearrange dimensions back to NROW x PCOL x KFACT
  return(aperm(M_result, c(3, 1, 2)))
}

#-----------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Xi  ###############
#-----------------------------------------------------------------------#

E_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  # Create 3D arrays by recycling matrices
  R_array <- array(new.R, dim = c(PCOL, KFACT, NROW))
  LAMBDA_array <- array(new.LAMBDA, dim = c(PCOL, KFACT, NROW))
  XI_array <- array(aperm(new.XI, c(2, 1)), dim = c(1, KFACT, NROW))
  OMEGA_array <- array(aperm(new.OMEGA, c(2, 1)), dim = c(1, KFACT, NROW))
  ALPHA_array <- array(aperm(new.ALPHA, c(2, 1)), dim = c(1, KFACT, NROW))

  # Align dimensions for vectorized computations
  XI_array <- XI_array[rep(1, PCOL), , ]
  OMEGA_array <- OMEGA_array[rep(1, PCOL), , ]
  ALPHA_array <- ALPHA_array[rep(1, PCOL), , ]

  numerator1 <- LAMBDA_array * XI_array + R_array
  denominator1 <- (1 - LAMBDA_array * OMEGA_array^2) + 1e-8
  term1 <- numerator1 / denominator1

  numerator2 <- LAMBDA_array * OMEGA_array * ALPHA_array
  denominator2 <- sqrt(1 - LAMBDA_array * OMEGA_array^2) * sqrt(1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2) + 1e-8
  term2 <- phi_by_Phi(R_array, LAMBDA_array, XI_array, OMEGA_array, ALPHA_array) * (numerator2 / denominator2)

  E_result <- term1 + term2

  # Rearrange dimensions back to NROW x PCOL x KFACT
  return(aperm(E_result, c(3, 1, 2)))
}

#--------------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Omega  ###############
#--------------------------------------------------------------------------#

O_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  # Create 3D arrays by recycling matrices
  R_array <- array(new.R, dim = c(PCOL, KFACT, NROW))
  LAMBDA_array <- array(new.LAMBDA, dim = c(PCOL, KFACT, NROW))
  XI_array <- array(aperm(new.XI, c(2, 1)), dim = c(1, KFACT, NROW))
  OMEGA_array <- array(aperm(new.OMEGA, c(2, 1)), dim = c(1, KFACT, NROW))
  ALPHA_array <- array(aperm(new.ALPHA, c(2, 1)), dim = c(1, KFACT, NROW))

  # Align dimensions for vectorized computations
  XI_array <- XI_array[rep(1, PCOL), , ]
  OMEGA_array <- OMEGA_array[rep(1, PCOL), , ]
  ALPHA_array <- ALPHA_array[rep(1, PCOL), , ]

  denominator1 <- (1 - LAMBDA_array * OMEGA_array^2)
  term1 <- (LAMBDA_array * OMEGA_array / (denominator1 + 1e-8))
  term2 <- ((R_array^2 * OMEGA_array + LAMBDA_array^2 * OMEGA_array * XI_array^2 + 2 * LAMBDA_array * OMEGA_array * R_array * XI_array) / (denominator1^2 + 1e-8))

  denominator2 <- sqrt(1 - LAMBDA_array * OMEGA_array^2) * sqrt(1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2)
  denominator3 <- ((1 - LAMBDA_array * OMEGA_array^2) * (1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2))^1.5
  term3 <- ((R_array + LAMBDA_array * XI_array) * ALPHA_array) / (denominator2 + 1e-8)
  term4 <- ((R_array + LAMBDA_array * XI_array) * OMEGA_array * ALPHA_array * (-2 * LAMBDA_array * OMEGA_array - LAMBDA_array * OMEGA_array * ALPHA_array^2 + 2 * LAMBDA_array^2 * OMEGA_array^3)) / (denominator3 + 1e-8)

  term5 <- phi_by_Phi(R_array, LAMBDA_array, XI_array, OMEGA_array, ALPHA_array) * (term3 - term4)

  O_result <- term1 + term2 + term5

  # Rearrange dimensions back to NROW x PCOL x KFACT
  return(aperm(O_result, c(3, 1, 2)))
}

#--------------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Alpha  ###############
#--------------------------------------------------------------------------#

A_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  # Create 3D arrays by recycling matrices
  R_array <- array(new.R, dim = c(PCOL, KFACT, NROW))
  LAMBDA_array <- array(new.LAMBDA, dim = c(PCOL, KFACT, NROW))
  XI_array <- array(aperm(new.XI, c(2, 1)), dim = c(1, KFACT, NROW))
  OMEGA_array <- array(aperm(new.OMEGA, c(2, 1)), dim = c(1, KFACT, NROW))
  ALPHA_array <- array(aperm(new.ALPHA, c(2, 1)), dim = c(1, KFACT, NROW))

  # Align dimensions for vectorized computations
  XI_array <- XI_array[rep(1, PCOL), , ]
  OMEGA_array <- OMEGA_array[rep(1, PCOL), , ]
  ALPHA_array <- ALPHA_array[rep(1, PCOL), , ]

  denominator2 <- sqrt(1 - LAMBDA_array * OMEGA_array^2) * sqrt(1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2)
  denominator3 <- ((1 - LAMBDA_array * OMEGA_array^2) * (1 + ALPHA_array^2 - LAMBDA_array * OMEGA_array^2))^1.5
  term3 <- ((R_array + LAMBDA_array * XI_array) * OMEGA_array) / (denominator2 + 1e-8)
  term4 <- ((R_array + LAMBDA_array * XI_array) * OMEGA_array * ALPHA_array * (ALPHA_array - LAMBDA_array * OMEGA_array^2 * ALPHA_array)) / (denominator3 + 1e-8)

  A_result <- phi_by_Phi(R_array, LAMBDA_array, XI_array, OMEGA_array, ALPHA_array) * (term3 - term4)

  # Rearrange dimensions back to NROW x PCOL x KFACT
  return(aperm(A_result, c(3, 1, 2)))
}
