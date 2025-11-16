#---------------------------------------------#
###############  Psi Functions  ###############
#---------------------------------------------#

psi_1 <- function(xi, omega, alpha, alpha_prior) {
  alpha_star <- -sqrt(2 / pi) * (alpha_prior / sqrt(1 + alpha_prior^2))

  A1 <- alpha_prior * (xi + (omega * alpha / (sqrt(1 + alpha^2) + 1e-8)) * sqrt(2 / pi) - alpha_star)

  A2 <- (2 / pi) * (alpha^2 / (1 + alpha^2 + 1e-8))

  M1 <- sqrt(0.5 / pi) * A1 * exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)

  M2 <- (0.5 / pi) * exp(-A1^2) / (pnorm(A1)^2 + 1e-8)

  Taylor_Approximated_Value <- pnorm(A1, log.p = TRUE) + 0.5 * alpha_prior^2 * omega^2 * (1 - A2) * (-M1 - M2)

  return(Taylor_Approximated_Value)
}

psi_1_derivative_xi <- function(xi, omega, alpha, alpha_prior) {
  alpha_star <- -sqrt(2 / pi) * (alpha_prior / sqrt(1 + alpha_prior^2))

  A1 <- alpha_prior * (xi + (omega * alpha / (sqrt(1 + alpha^2) + 1e-8)) * sqrt(2 / pi) - alpha_star)

  A2 <- (2 / pi) * (alpha^2 / (1 + alpha^2 + 1e-8))

  dA1 <- alpha_prior

  dM1 <- dA1 * (sqrt(0.5 / pi) * (exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8))) * (1 - A1^2 - A1 * (dnorm(A1) / (pnorm(A1) + 1e-8)))

  dM2 <- dA1 * (-1 / pi) * (exp(-A1^2) / (pnorm(A1)^2 + 1e-8)) * (A1 + (dnorm(A1) / (pnorm(A1) + 1e-8)))

  psi_deriv_xi <- dA1 * (dnorm(A1) / (pnorm(A1) + 1e-8)) + 0.5 * alpha_prior^2 * omega^2 * (1 - A2) * (-dM1 - dM2)

  return(psi_deriv_xi)
}

psi_1_derivative_omega <- function(xi, omega, alpha, alpha_prior) {
  alpha_star <- -sqrt(2 / pi) * (alpha_prior / sqrt(1 + alpha_prior^2))

  A1 <- alpha_prior * (xi + (omega * alpha / (sqrt(1 + alpha^2) + 1e-8)) * sqrt(2 / pi) - alpha_star)

  A2 <- (2 / pi) * (alpha^2 / (1 + alpha^2 + 1e-8))

  dA1 <- sqrt(2 / pi) * (alpha_prior * alpha / (sqrt(1 + alpha^2) + 1e-8))

  M1 <- sqrt(0.5 / pi) * A1 * exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)

  M2 <- (0.5 / pi) * exp(-A1^2) / (pnorm(A1)^2 + 1e-8)

  dM1 <- dA1 * sqrt(0.5 / pi) * (exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)) * (1 - A1^2 - A1 * (dnorm(A1) / (pnorm(A1) + 1e-8)))

  dM2 <- dA1 * (-1 / pi) * (exp(-A1^2) / (pnorm(A1)^2 + 1e-8)) * (A1 + (dnorm(A1) / (pnorm(A1) + 1e-8)))

  psi_deriv_omega <- dA1 * (dnorm(A1) / (pnorm(A1) + 1e-8)) + alpha_prior^2 * omega * (1 - A2) * (-M1 - M2) + 0.5 * alpha_prior^2 * omega^2 * (1 - A2) * (-dM1 - dM2)

  return(psi_deriv_omega)
}

psi_1_derivative_alpha <- function(xi, omega, alpha, alpha_prior) {
  alpha_star <- -sqrt(2 / pi) * (alpha_prior / sqrt(1 + alpha_prior^2))

  A1 <- alpha_prior * (xi + (omega * alpha / (sqrt(1 + alpha^2) + 1e-8)) * sqrt(2 / pi) - alpha_star)

  A2 <- (2 / pi) * (alpha^2 / (1 + alpha^2 + 1e-8))

  dA1 <- alpha_prior * sqrt(2 / pi) * omega / ((1 + alpha^2)^1.5 + 1e-8)

  dA2 <- (4 / pi) * (alpha / ((1 + alpha^2)^2 + 1e-8))

  M1 <- sqrt(0.5 / pi) * A1 * exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)

  M2 <- (0.5 / pi) * exp(-A1^2) / (pnorm(A1)^2 + 1e-8)

  dM1 <- dA1 * sqrt(0.5 / pi) * (exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)) * (1 - A1^2 - A1 * (dnorm(A1) / (pnorm(A1) + 1e-8)))

  dM2 <- dA1 * (-1 / pi) * (exp(-A1^2) / (pnorm(A1)^2 + 1e-8)) * (A1 + (dnorm(A1) / (pnorm(A1) + 1e-8)))

  psi_deriv_alpha <- dA1 * (dnorm(A1) / (pnorm(A1) + 1e-8)) + 0.5 * alpha_prior^2 * omega^2 * (-dA2) * (-M1 - M2) + 0.5 * alpha_prior^2 * omega^2 * (1 - A2) * (-dM1 - dM2)

  return(psi_deriv_alpha)
}

psi_2 <- function(alpha) {
  A1 <- sqrt(2 / pi) * ((alpha^2) / (sqrt(1 + alpha^2) + 1e-8))

  N1 <- sqrt(0.5 / pi) * A1 * exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)

  N2 <- (0.5 / pi) * exp(-A1^2) / (pnorm(A1)^2 + 1e-8)

  Taylor_Approximated_Value <- pnorm(A1, log.p = TRUE) + (alpha^2 / 2 - A1^2 / 2) * (-N1 - N2)

  return(Taylor_Approximated_Value)
}

psi_2_derivative <- function(alpha) {
  A1 <- sqrt(2 / pi) * ((alpha^2) / (sqrt(1 + alpha^2) + 1e-8))

  I1 <- sqrt(2 / pi) * ((2 * alpha + alpha^3) / ((1 + alpha^2)^1.5 + 1e-8))

  N1 <- sqrt(0.5 / pi) * A1 * exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)

  N2 <- (0.5 / pi) * exp(-A1^2) / (pnorm(A1)^2 + 1e-8)

  dN1 <- I1 * sqrt(0.5 / pi) * (exp(-0.5 * A1^2) / (pnorm(A1) + 1e-8)) * (1 - A1^2 - A1 * (dnorm(A1) / (pnorm(A1) + 1e-8)))

  dN2 <- I1 * (-1 / pi) * (exp(-A1^2) / (pnorm(A1)^2 + 1e-8)) * (A1 + (dnorm(A1) / (pnorm(A1) + 1e-8)))

  psi_deriv <- I1 * (dnorm(A1) / (pnorm(A1) + 1e-8)) + (alpha - A1 * I1) * (-N1 - N2) + (alpha^2 / 2 - A1^2 / 2) * (-dN1 - dN2)

  return(psi_deriv)
}
