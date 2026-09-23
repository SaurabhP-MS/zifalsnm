mills <- function(x) {
  h <- exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
  i <- x < -30
  if (any(i)) h[i] <- -x[i] - 1/x[i] + 2/x[i]^3     
  h
}


#-------------------------------------------------#
###############  L_ij Function  ###################
#-------------------------------------------------#

L_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL){
  
  jj <- rep(1:PCOL, each = NROW)
  ii <- rep(1:NROW, PCOL)
  
  lam2 <- new.LAMBDA[jj, , drop = FALSE]
  r    <- new.R[jj, , drop = FALSE]
  xi   <- new.XI[ii, , drop = FALSE]
  om   <- new.OMEGA[ii, , drop = FALSE]
  al   <- new.ALPHA[ii, , drop = FALSE]
  
  A <- 1 - lam2 * om^2
  B <- A + al^2                       
  
  num1 <- lam2 * xi^2 + r^2 * om^2 + 2 * r * xi
  g    <- (r + lam2 * xi) * om * al / sqrt(A * B)
  
  log_term   <- log(2) - 0.5 * log(A) + num1 / (2 * A)
  pnorm_term <- pnorm(g, log.p = TRUE)
  
  matrix(rowSums(log_term),   nrow = NROW, byrow = FALSE) +
    matrix(rowSums(pnorm_term), nrow = NROW, byrow = FALSE)
}


.lij_core <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA,
                      NROW, PCOL, KFACT) {
  
  r    <- array(new.R,      dim = c(PCOL, KFACT, NROW))
  lam2 <- array(new.LAMBDA, dim = c(PCOL, KFACT, NROW))
  
  xi <- aperm(array(new.XI,    dim = c(NROW, KFACT, PCOL)), c(3, 2, 1))
  om <- aperm(array(new.OMEGA, dim = c(NROW, KFACT, PCOL)), c(3, 2, 1))
  al <- aperm(array(new.ALPHA, dim = c(NROW, KFACT, PCOL)), c(3, 2, 1))
  
  A      <- 1 - lam2 * om^2
  B      <- A + al^2
  sqrtAB <- sqrt(A * B)
  g      <- (r + lam2 * xi) * om * al / sqrtAB
  
  list(r = r, lam2 = lam2, xi = xi, om = om, al = al,
       A = A, B = B, sqrtAB = sqrtAB, AB15 = sqrtAB^3,
       h = mills(g))                     
}


#-------------------------------------------------#
###############  A Common Function  ###############
#-------------------------------------------------#

phi_by_Phi <- function(R_value, LAMBDA_value, XI_value, OMEGA_value, ALPHA_value) {
  A <- 1 - LAMBDA_value * OMEGA_value^2
  B <- A + ALPHA_value^2
  mills((R_value + LAMBDA_value * XI_value) * OMEGA_value * ALPHA_value / sqrt(A * B))
}


#----------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to R  ###############
#----------------------------------------------------------------------#

R_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  
  z <- .lij_core(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT)
  
  term1 <- (z$r * z$om^2 + z$xi) / z$A
  term2 <- z$h * (z$om * z$al / z$sqrtAB)
  
  aperm(term1 + term2, c(3, 1, 2))       
}


#---------------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Lambda  ###############
#---------------------------------------------------------------------------#

M_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  
  z <- .lij_core(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT)
  
  term1 <- 0.5 * z$om^2 / z$A
  term2 <- 0.5 * (z$xi^2 + z$r^2 * z$om^4 + 2 * z$r * z$xi * z$om^2) / z$A^2
  
  term3 <- z$xi * z$om * z$al / z$sqrtAB
  
  term4 <- (z$r + z$lam2 * z$xi) * z$om * z$al *
    (-z$om^2 * (z$A + z$B)) / (2 * z$AB15)
  term5 <- z$h * (term3 - term4)
  
  aperm(term1 + term2 + term5, c(3, 1, 2))
}


#-----------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Xi  ###############
#-----------------------------------------------------------------------#

E_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  
  z <- .lij_core(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT)
  
  term1 <- (z$lam2 * z$xi + z$r) / z$A
  term2 <- z$h * (z$lam2 * z$om * z$al / z$sqrtAB)
  
  aperm(term1 + term2, c(3, 1, 2))
}


#--------------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Omega  ###############
#--------------------------------------------------------------------------#

O_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  
  z <- .lij_core(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT)
  
  term1 <- z$lam2 * z$om / z$A
  term2 <- (z$r^2 * z$om + z$lam2^2 * z$om * z$xi^2 +
              2 * z$lam2 * z$om * z$r * z$xi) / z$A^2
  
  term3 <- (z$r + z$lam2 * z$xi) * z$al / z$sqrtAB

  term4 <- (z$r + z$lam2 * z$xi) * z$om * z$al *
    (-z$lam2 * z$om * (z$A + z$B)) / z$AB15
  term5 <- z$h * (term3 - term4)
  
  aperm(term1 + term2 + term5, c(3, 1, 2))
}


#--------------------------------------------------------------------------#
###############  Derivative of L_ij Function w.r.t to Alpha  ###############
#--------------------------------------------------------------------------#

A_ij <- function(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT) {
  
  z <- .lij_core(new.R, new.LAMBDA, new.XI, new.OMEGA, new.ALPHA, NROW, PCOL, KFACT)
  
  term3 <- (z$r + z$lam2 * z$xi) * z$om / z$sqrtAB
  
  term4 <- (z$r + z$lam2 * z$xi) * z$om * z$al * (z$al * z$A) / z$AB15
  
  aperm(z$h * (term3 - term4), c(3, 1, 2))
}
