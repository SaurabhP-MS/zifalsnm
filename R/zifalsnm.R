###### ZIFA-LSN Model ######

#' Zero Inflated Factor Analysis Logistic Skew Normal Multinomial (ZIFA-LSNM) model
#'
#' @description ZIFA-LSNM model's objective is to achieve effective dimension
#' reduction to manage high dimensionality, accounting for the microbiome
#' data compositional nature, zero inflation and, it utilizes
#' skew-normal priors on the latent factors to explicitly model skewness
#' in log-ratio transformation. The model is fitted by coordinate ascent
#' variational inference (CAVI).
#'
#' @param X Count matrix with samples in rows and taxa in columns. The last
#'   column is used as the reference taxon for the log-ratio transformation.
#'   Taxa with only zero counts should be removed before fitting the model.
#' @param number_of_factors Number of latent factors \eqn{K}.
#' @param epsilon Convergence tolerance: the algorithm stops when the absolute
#'   relative change in the ELBO between two iterations is below `epsilon`.
#' @param NU1Prior Prior value for the first shape parameter of \eqn{\kappa_j}.
#' @param NU2Prior Prior value for the second shape parameter of \eqn{\kappa_j}.
#' @param G1Prior Prior value for the shape parameter of \eqn{\delta_{jt}}.
#' @param G2Prior Prior value for the rate parameter of \eqn{\delta_{jt}}.
#' @param AlphaPrior Prior value for the skewness parameter of the latent
#'   factors \eqn{F_{it}}.
#' @param Max_Iter Maximum number of CAVI iterations.
#' @param verbose If `TRUE`, progress messages are printed at every iteration.
#'
#' @return A list with the optimized variational parameters:
#' \describe{
#'   \item{Pi_ij}{Estimated zero-inflation indicators (\eqn{N \times (P-1)}).}
#'   \item{Tau_1, Tau_2}{Variational Beta parameters for \eqn{\kappa_j}.}
#'   \item{R, Lambda}{Variational mean and variance of the factor loadings.}
#'   \item{G1, G2}{Variational Gamma parameters for \eqn{\delta_{jt}}.}
#'   \item{Xi, Omega, Alpha}{Variational skew-normal location, scale and
#'     skewness of the latent factors.}
#'   \item{A0, C0}{Variational mean and variance of the intercepts \eqn{\beta_{0j}}.}
#'   \item{Estimated_Compositions}{Estimated compositions (\eqn{N \times P});
#'     the last column is the reference taxon.}
#'   \item{ELBO_Trace}{ELBO value at each iteration.}
#'   \item{*_Trace}{Values of each variational parameter at each iteration.}
#' }
#'
#' @importFrom sn rsn
#' @import stats
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 50 # Number of samples
#' P <- 50 # Number of taxa
#' K <- 2  # Number of factors
#'
#' # Latent factors
#' alpha_prior <- 2
#' alpha_star <- -sqrt(2 / pi) * (alpha_prior / sqrt(1 + alpha_prior^2))
#' f <- matrix(sn::rsn(N * K, alpha_star, 1, alpha_prior), N, K, TRUE)
#'
#' # Factor loadings
#' beta <- matrix(runif(P * K, -2.5, 2.5), nrow = P, ncol = K, byrow = TRUE)
#' kappa <- runif(P)
#' beta_0 <- runif(P, -2, 2)
#'
#' z <- matrix(0, N, P, byrow = TRUE)
#' for (i in 1:N) {
#'   z[i, ] <- rbinom(P, size = 1, prob = kappa)
#' }
#' linear_predictor <- matrix(beta_0, N, P, byrow = TRUE) + tcrossprod(f, beta)
#' rho <- (1 - z) * exp(linear_predictor) /
#'   (rowSums((1 - z) * exp(linear_predictor)) + 1e-8)
#'
#' cmat <- matrix(0, N, P)
#' for (i in 1:N) {
#'   cmat[i, ] <- rmultinom(1, as.integer(runif(1, 80000, 100000)), prob = rho[i, ])
#' }
#'
#' # Remove taxa with only zero counts before fitting the model
#' cmat <- cmat[, colSums(cmat) > 0]
#'
#' OUTPUT <- ZIFA_LSNM(X = cmat, number_of_factors = K, epsilon = 5e-5,
#'                     NU1Prior = 1, NU2Prior = 3, G1Prior = 200,
#'                     G2Prior = 1 / 200, AlphaPrior = 2, Max_Iter = 200)
#' }

ZIFA_LSNM <- function(X, number_of_factors = 2, epsilon = 5e-5, NU1Prior = 1,
                      NU2Prior = 3, G1Prior = 200, G2Prior = 1 / 200,
                      AlphaPrior = 2, Max_Iter = 200, verbose = TRUE){

  P <- ncol(X)
  x<-X[,-P]
  n <- nrow(x)
  p <- ncol(x)
  k <- number_of_factors
  M <- rowSums(X)
  nu1 <- NU1Prior
  nu2 <- NU2Prior
  g1 <- G1Prior
  g2 <- G2Prior
  prior_alpha <- AlphaPrior

  # Initial Guess For The Variational Parameters

  X_transformation <- scale(log2(x+0.05), scale=T, center=T)
  Decomposition_X <- svd(X_transformation, k, k)
  
  r <- Decomposition_X$v
  lambda <- matrix(0.65,p,k)
  g_1 <- matrix(2,p,k)
  g_2 <- matrix(2,p,k)
  tau_1 <-  rep(2,p)
  tau_2 <- rep(2,p)
  xi <- Decomposition_X$u %*% diag(Decomposition_X$d[1:k])
  omega <- matrix(0.45,n,k)
  alpha <- matrix(prior_alpha,n,k)
  zero_ratio <- apply(x, 2, function(y) {sum(y==0)/n})
  pi_ij <- updated_pi_ij <- t(ifelse(t(x)==0, zero_ratio, 0))
  pi_matrix <- matrix(0,n,p)
  a_0j <- runif(p)
  c_0j <- rep(0.35,p)


  old_ELBO <- NA_real_
  elbo_trace <- pi_list <- tau1_list <- tau2_list <- list()
  r_list <- lambda_list <- g1_list <- g2_list <- xi_list <- omega_list <- list()
  alpha_list <- a0_list <- c0_list <- list()

  for(iter in 1:Max_Iter) {

    TomAndJerry <- matrix((a_0j + 0.5*c_0j),n,p,byrow=TRUE) +
      L_ij(new.R=r, new.LAMBDA=lambda, new.XI=xi,
           new.OMEGA=omega, new.ALPHA=alpha, NROW=n, PCOL=p)

    upsilon <- log(M/(rowSums((1-updated_pi_ij)*exp(TomAndJerry))))
    for (i in 1:n){
      for (j in 1:p){
        pi_matrix[i,j] <- exp(digamma(tau_1[j])) / (exp(digamma(tau_1[j])) +
                                                              exp(digamma(tau_2[j]) - exp(upsilon[i] + TomAndJerry[i,j])))
      }
    }
    pi_matrix[x!=0] <- 0
    updated_pi_ij <- ifelse(pi_matrix > 0.5, 1, 0)
    pi_list[[iter]] <- updated_pi_ij
    
    
    if (verbose) cat("Pi_ij updated on iteration step : ",iter,"\n")

    # Optimization For Tau_j1 & Tau_j2

    for(j in 1:p) {
      tau_1[j]<- optim(tau_1[j], UPDATED_TAU_2=tau_2,
                               PIIJ=updated_pi_ij, NU_1=nu1, NU_2=nu2, l=j, N=n,
                               tau_j1_function, tau_j1_gradient,
                               method = "L-BFGS-B",lower=0.1,
                               upper= +Inf,control=list(fnscale=-1))$par
    }
    
    tau1_list[[iter]] <- tau_1
    
    if (verbose) cat("Tau_j1 updated on iteration step : ",iter,"\n")

    for(j in 1:p) {
      tau_2[j]<- optim(tau_2[j], UPDATED_TAU_1=tau_1,
                               PIIJ=updated_pi_ij, NU_1=nu1, NU_2=nu2, l=j, N=n,
                               tau_j2_function, tau_j2_gradient,
                               method = "L-BFGS-B",lower=0.1, upper= +Inf,
                               control=list(fnscale=-1))$par
    }
    
    tau2_list[[iter]] <- tau_2
    
    if (verbose) cat("Tau_j2 updated on iteration step : ",iter,"\n")

    # Optimization For G1 & G2
    
    q21 <- try(optim(c(g_1), G2=g_2, R=r, LAMBDA=lambda,
                     G1_Prior=g1, G2_Prior=g2, N=n, P=p, K=k,method="L-BFGS-B",
                     lower=matrix(1e-2,p,k), upper=matrix(+Inf,p,k),
                     fn=g1_function, gr=g1_gradient,
                     control=list(fnscale=-1,maxit=1000)), silent = TRUE)
    if("try-error" %in% class(q21)){

      if (verbose) cat("G1 cannot be updated on iteration step : ",iter,"\n")
    }else{

      g_1 <- matrix(q21$par, p, k)

      if (verbose) cat("G1 updated on iteration step : ",iter,"\n")
    }
    g1_list[[iter]] <- g_1

    q22 <- try(optim(c(g_2), G1=g_1, R=r, LAMBDA=lambda,
                     G1_Prior=g1, G2_Prior=g2, N=n, P=p, K=k,method="L-BFGS-B",
                     lower=matrix(1e-2,p,k), upper=matrix(+Inf,p,k),
                     fn=g2_function, gr=g2_gradient,
                     control=list(fnscale=-1,maxit=1000)), silent = TRUE)
    if("try-error" %in% class(q22)){

      if (verbose) cat("G2 cannot be updated on iteration step : ",iter,"\n")
    }else{

      g_2 <- matrix(q22$par, p, k)

      if (verbose) cat("G2 updated on iteration step : ",iter,"\n")
    }
    g2_list[[iter]] <- g_2

    # Optimization For R

    q1 <- try(optim(c(r), COUNT_MATRIX=x, LAMBDA=lambda, XI=xi,
                    OMEGA=omega, ALPHA=alpha, G1=g_1,
                    G2=g_2, A0=a_0j, C0=c_0j,
                    PIIJ=updated_pi_ij, N=n, P=p, K=k, M=M, method="BFGS",
                    fn=r_function, gr=r_gradient,
                    control=list(fnscale=-1,maxit=1000)), silent=TRUE)

    if("try-error" %in% class(q1)){

      if (verbose) cat("R cannot be updated on iteration step : ",iter,"\n")
    }else{

      r <- matrix(q1$par, p, k)

      if (verbose) cat("R updated on iteration step : ",iter,"\n")
    }
    r_list[[iter]] <- r

    # Optimization For Lambda

    q2 <- try(optim(c(lambda), R=r, XI=xi, OMEGA=omega,
                    ALPHA=alpha, G1=g_1, G2=g_2,
                    A0=a_0j, C0=c_0j, PIIJ=updated_pi_ij, N=n,
                    P=p, K=k, M=M, method="L-BFGS-B", lower=matrix(1e-8,p,k),
                    upper=matrix(1,p,k), fn=lambda_function, gr=lambda_gradient,
                    control=list(fnscale=-1,maxit=1000)), silent = TRUE)

    if("try-error" %in% class(q2)){

      if (verbose) cat("Lambda cannot be updated on iteration step : ",iter,"\n")
    }else{

      lambda <- matrix(q2$par, p, k)

      if (verbose) cat("Lambda updated on iteration step : ",iter,"\n")
    }
    lambda_list[[iter]] <- lambda

    # Optimization For Xi

    q3 <- try(optim(c(xi), COUNT_MATRIX=x, R=r, LAMBDA=lambda,
                    OMEGA=omega, ALPHA=alpha, A0=a_0j,
                    C0=c_0j, PIIJ=updated_pi_ij, ALPHA_Prior=prior_alpha,
                    N=n, P=p, K=k, M=M, method="BFGS", fn=xi_function,
                    gr=xi_gradient, control=list(fnscale=-1,maxit=1000)), silent=TRUE)

    if("try-error" %in% class(q3)){

      if (verbose) cat("Xi cannot be updated on iteration step : ",iter,"\n")
    }else{

      xi <- matrix(q3$par, n, k)

      if (verbose) cat("Xi updated on iteration step : ",iter,"\n")
    }
    xi_list[[iter]] <- xi
    # Optimization For Omega

    q4 <- try(optim(c(omega), COUNT_MATRIX=x, R=r, LAMBDA=lambda,
                    XI=xi, ALPHA=alpha, A0=a_0j,
                    C0=c_0j, PIIJ=updated_pi_ij, ALPHA_Prior=prior_alpha,
                    N=n, P=p, K=k, M=M, method="L-BFGS-B", lower=1e-8, upper=1,
                    fn=omega_function, gr=omega_gradient,
                    control=list(fnscale=-1,maxit=1000)), silent = TRUE)

    if("try-error" %in% class(q4)){

      if (verbose) cat("Omega cannot be updated on iteration step : ",iter,"\n")

    }else{

      omega <- matrix(q4$par, n, k)

      if (verbose) cat("Omega updated on iteration step : ",iter,"\n")
    }
    omega_list[[iter]] <- omega
    # Optimization For Alpha

    q5 <- try(optim(c(alpha), COUNT_MATRIX=x, R=r, LAMBDA=lambda,
                    XI=xi, OMEGA=omega, A0=a_0j,
                    C0=c_0j, PIIJ=updated_pi_ij, ALPHA_Prior=prior_alpha,
                    N=n, P=p, K=k, M=M, method="BFGS", fn=alpha_function,
                    gr=alpha_gradient, control=list(fnscale=-1,maxit=1000)), silent = TRUE)

    if("try-error" %in% class(q5)){

      if (verbose) cat("Alpha cannot be updated on iteration step : ",iter,"\n")
    }else{

      alpha <- matrix(q5$par, n, k)

      if (verbose) cat("Alpha updated on iteration step : ",iter,"\n")
    }
    alpha_list[[iter]] <- alpha

    # Optimization For A0

    q6 <- try(optim(c(a_0j), COUNT_MATRIX=x,R=r, LAMBDA=lambda,
                    XI=xi, OMEGA=omega, ALPHA=alpha,
                    PIIJ=updated_pi_ij, C0=c_0j,N=n,P=p,K=k,M=M,
                    method="BFGS", fn=a0_function, gr=a0_gradient,
                    control=list(fnscale=-1, maxit=1000)), silent = TRUE)
    if("try-error" %in% class(q6)){

      if (verbose) cat("A0 cannot be updated on iteration step : ",iter,"\n")
    }else{

      a_0j <- q6$par

      if (verbose) cat("A0 updated on iteration step : ",iter,"\n")
    }
    a0_list[[iter]] <- a_0j

    # Optimization For C0

    q7 <- try(optim(c(c_0j),R=r, LAMBDA=lambda, XI=xi,
                    OMEGA=omega, ALPHA=alpha, PIIJ=updated_pi_ij,
                    A0=a_0j,N=n,P=p,K=k,M=M,method="L-BFGS-B", lower=1e-5,
                    upper=+Inf ,fn=c0_function, gr=c0_gradient,
                    control=list(fnscale=-1, maxit=1000)), silent = TRUE)
    if("try-error" %in% class(q7)){

      if (verbose) cat("C0 cannot be updated on iteration step : ",iter,"\n")
    }else{

      c_0j <- q7$par

      if (verbose) cat("C0 updated on iteration step : ",iter,"\n")
    }
    c0_list[[iter]] <- c_0j

    # Stopping Rule
    
    new_ELBO <- ELBO_AP(COUNT_MATRIX = x, R = r, LAMBDA = lambda, XI = xi,
                        OMEGA = omega, ALPHA = alpha, PIIJ = updated_pi_ij,
                        A0 = a_0j, C0 = c_0j, T1 = tau_1, T2 = tau_2,
                        G1 = g_1, G2 = g_2, G1_Prior = g1, G2_Prior = g2,
                        NU1 = nu1, NU2 = nu2, ALPHA_Prior = prior_alpha,
                        N = n, P = p, K = k, M = M)
    
    elbo_trace[[iter]] <- new_ELBO
    
    if (iter == 1) {
      if (verbose) cat("Iteration", iter, ": ELBO =", format(new_ELBO), "\n")
    } else {
      increment <- new_ELBO - old_ELBO            
      relative  <- increment / abs(old_ELBO)
      if (verbose) cat("Old ELBO :",old_ELBO,",New ELBO :",new_ELBO,
          ",Difference of ELBO's :",increment,",Relative Difference :",relative,"\n")
      if (abs(relative) < epsilon) { if (verbose) cat("Converged at iteration", iter, "\n"); break }
    }
    old_ELBO <- new_ELBO
    
    if (iter == Max_Iter) if (verbose) cat("Max iteration reached")
  }

  # Estimated compositions

  DonalDuck <-matrix((a_0j + 0.5*c_0j),n,p,byrow=TRUE) +
    L_ij(new.R=r, new.LAMBDA=lambda, new.XI=xi,
         new.OMEGA=omega, new.ALPHA=alpha, NROW=n, PCOL=p)
  denom <- 1 + rowSums(exp(DonalDuck))
  comp_main <- exp(DonalDuck) / denom
  comp_ref <- 1 / denom
  estimated_compositions <- cbind(comp_main, comp_ref)

  
  list(Pi_ij = updated_pi_ij,
       Tau_1 = tau_1, Tau_2 = tau_2,
       R = r, Lambda = lambda, G1 = g_1, G2 = g_2,
       Xi = xi, Omega = omega, Alpha = alpha,
       A0 = a_0j, C0 = c_0j,
       Estimated_Compositions = estimated_compositions,
       ELBO_Trace = elbo_trace,
       G1_Trace = g1_list, G2_Trace = g2_list,
       R_Trace = r_list, Lambda_Trace = lambda_list,
       Xi_Trace = xi_list, Omega_Trace = omega_list,
       Alpha_Trace = alpha_list, A0_Trace = a0_list,
       C0_Trace = c0_list, Tau1_Trace = tau1_list, Tau2_Trace = tau2_list)

}
