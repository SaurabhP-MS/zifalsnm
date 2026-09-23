xlogx <- function(u) {
  out <- numeric(length(u))
  i <- u > 0
  out[i] <- u[i] * log(u[i])
  dim(out) <- dim(u)
  out
}



ELBO_AP <- function(COUNT_MATRIX,R,LAMBDA,XI,OMEGA,ALPHA,PIIJ,A0,C0,T1,T2,G1,G2,
                    G1_Prior,G2_Prior,NU1,NU2,ALPHA_Prior,N,P,K,M){
  
  new.tau1 <- T1
  new.tau2 <- T2
  A0_vector <- matrix(A0,P,1)
  alpha_star <- -sqrt(2/pi)*(ALPHA_Prior/sqrt(1 + ALPHA_Prior^2))
  
  # Term 11
  term11 <- COUNT_MATRIX %*% A0_vector
  
  # Term 12
  xi_omega_alpha_matrix <- XI + (OMEGA*ALPHA / sqrt(1+ALPHA^2))*sqrt(2/pi)
  term12 <- COUNT_MATRIX*(tcrossprod(xi_omega_alpha_matrix, R))
  
  # Term 2
  TomJerry <- matrix((A0 + 0.5*C0),N,P,byrow=TRUE) +
    L_ij(new.R=R, new.LAMBDA=LAMBDA, new.XI=XI, new.OMEGA=OMEGA,
         new.ALPHA=ALPHA, NROW=N, PCOL=P)
  term2 <- M*log(1+rowSums((1-PIIJ)*exp(TomJerry)))
  
  # Term 3
  dg1  <- digamma(new.tau1)
  dg2  <- digamma(new.tau2)
  dg12 <- digamma(new.tau1 + new.tau2)
  term3 <- sum((NU1 - new.tau1)*(dg1 - dg12) +
                 (NU2 - new.tau2)*(dg2 - dg12) +
                 lbeta(new.tau1, new.tau2))
  
  # Term 4
  E1 <- matrix(dg1 - dg12, N, P, byrow = TRUE)
  E2 <- matrix(dg2 - dg12, N, P, byrow = TRUE)
  term4 <- sum(PIIJ*E1 + (1 - PIIJ)*E2 - xlogx(PIIJ) - xlogx(1 - PIIJ))
  
  # Term 5
  term5 <- 0.5*OMEGA^2 + 0.5*XI^2 +
    sqrt(2/pi)*((XI*OMEGA*ALPHA)/sqrt(1+ALPHA^2)) -
    alpha_star*XI -
    sqrt(2/pi)*((alpha_star*OMEGA*ALPHA)/sqrt(1+ALPHA^2)) -
    psi_1(XI, OMEGA, ALPHA, ALPHA_Prior) + psi_2(ALPHA) - log(OMEGA)
  
  # Term 6
  term6 <- A0^2 + C0 - log(C0)
  
  # Term 7
  term7 <- (digamma(G1) - log(G2))*(G1_Prior - 0.5) -
    (G1/G2)*(0.5*LAMBDA + 0.5*R^2 + G2_Prior) +
    0.5*log(LAMBDA) + G1 - log(G2) + lgamma(G1) + (1 - G1)*digamma(G1)
  
  ELBO_value <- sum(term11) + sum(term12) - sum(term2) + term3 + term4 -
    sum(term5) - 0.5*sum(term6) + sum(term7)
  return(ELBO_value)
}
