# mills() is defined in Lij_and_derivatives.R

mills_pieces <- function(x) {
  h  <- mills(x)
  M1 <- x * h
  M2 <- h * h
  list(h   = h,
       M1  = M1,
       M2  = M2,
       dM1 = h * (1 - x^2 - x * h),
       dM2 = -2 * h^2 * (x + h))
}

psi_1 <- function(xi, omega, alpha, alpha_prior) {
  
  alpha_star <- -sqrt(2/pi) * (alpha_prior / sqrt(1 + alpha_prior^2))
  
  mu   <- alpha_prior * (xi + sqrt(2/pi) * omega * alpha / sqrt(1 + alpha^2) -
                           alpha_star)
  sig2 <- alpha_prior^2 * omega^2 * (1 - (2/pi) * alpha^2 / (1 + alpha^2))
  
  mp <- mills_pieces(mu)
  
  Taylor_Approximated_Value <-
    pnorm(mu, log.p = TRUE) + 0.5 * sig2 * (-mp$M1 - mp$M2)
  
  return(Taylor_Approximated_Value)
}


psi_1_derivative_xi <- function(xi, omega, alpha, alpha_prior) {
  
  alpha_star <- -sqrt(2/pi) * (alpha_prior / sqrt(1 + alpha_prior^2))
  
  mu   <- alpha_prior * (xi + sqrt(2/pi) * omega * alpha / sqrt(1 + alpha^2) -
                           alpha_star)
  sig2 <- alpha_prior^2 * omega^2 * (1 - (2/pi) * alpha^2 / (1 + alpha^2))
  
  dmu <- alpha_prior
  
  mp <- mills_pieces(mu)
  
  psi_deriv_xi <- dmu * mp$h +
    0.5 * sig2 * (-mp$dM1 - mp$dM2) * dmu
  
  return(psi_deriv_xi)
}


psi_1_derivative_omega <- function(xi, omega, alpha, alpha_prior) {
  
  alpha_star <- -sqrt(2/pi) * (alpha_prior / sqrt(1 + alpha_prior^2))
  
  A2   <- (2/pi) * alpha^2 / (1 + alpha^2)
  mu   <- alpha_prior * (xi + sqrt(2/pi) * omega * alpha / sqrt(1 + alpha^2) -
                           alpha_star)
  sig2 <- alpha_prior^2 * omega^2 * (1 - A2)
  
  dmu   <- alpha_prior * sqrt(2/pi) * alpha / sqrt(1 + alpha^2)
  dsig2 <- 2 * alpha_prior^2 * omega * (1 - A2)
  
  mp <- mills_pieces(mu)
  
  psi_deriv_omega <- dmu * mp$h +
    0.5 * dsig2 * (-mp$M1  - mp$M2) +
    0.5 * sig2  * (-mp$dM1 - mp$dM2) * dmu
  
  return(psi_deriv_omega)
}


psi_1_derivative_alpha <- function(xi, omega, alpha, alpha_prior) {
  
  alpha_star <- -sqrt(2/pi) * (alpha_prior / sqrt(1 + alpha_prior^2))
  
  A2   <- (2/pi) * alpha^2 / (1 + alpha^2)
  mu   <- alpha_prior * (xi + sqrt(2/pi) * omega * alpha / sqrt(1 + alpha^2) -
                           alpha_star)
  sig2 <- alpha_prior^2 * omega^2 * (1 - A2)
  
  dmu   <- alpha_prior * sqrt(2/pi) * omega / (1 + alpha^2)^1.5
  dA2   <- (4/pi) * alpha / (1 + alpha^2)^2
  dsig2 <- -alpha_prior^2 * omega^2 * dA2
  
  mp <- mills_pieces(mu)
  
  psi_deriv_alpha <- dmu * mp$h +
    0.5 * dsig2 * (-mp$M1  - mp$M2) +
    0.5 * sig2  * (-mp$dM1 - mp$dM2) * dmu
  
  return(psi_deriv_alpha)
}


psi_2 <- function(alpha) {
  
  muV <- sqrt(2/pi) * alpha^2 / sqrt(1 + alpha^2)
  s2V <- alpha^2 - muV^2
  
  mp <- mills_pieces(muV)
  
  Taylor_Approximated_Value <-
    pnorm(muV, log.p = TRUE) + 0.5 * s2V * (-mp$M1 - mp$M2)
  
  return(Taylor_Approximated_Value)
}


psi_2_derivative <- function(alpha) {
  
  muV <- sqrt(2/pi) * alpha^2 / sqrt(1 + alpha^2)
  s2V <- alpha^2 - muV^2
  
  dmuV <- sqrt(2/pi) * (2 * alpha + alpha^3) / (1 + alpha^2)^1.5
  
  mp <- mills_pieces(muV)
  
  psi_deriv <- dmuV * mp$h +
    (alpha - muV * dmuV) * (-mp$M1  - mp$M2) +
    0.5 * s2V            * (-mp$dM1 - mp$dM2) * dmuV
  
  return(psi_deriv)
}