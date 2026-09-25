data_dir <- "C:/Work Files/Simulation_Study_For_Proj_1/K=5/Alpha_true = 0/Sim_100_50_5_0"
setwd("C:/Work Files/Simulation_Study_For_Proj_1/K=5/Alpha_true = 0/Sim_100_50_5_0")
N <- 100
P <- 50
K <- 5
files1 <- list.files(
  path       = data_dir,
  pattern    = "^Output_ZIFA_LSN_100_50_5_0_0_\\d+\\.RData$",
  full.names = TRUE
)

files2 <- list.files(
  path       = data_dir,
  pattern    = "^Output_ZIFA_LSN_100_50_5_0_2_\\d+\\.RData$",
  full.names = TRUE
)

files3 <- list.files(
  path       = data_dir,
  pattern    = "^Output_ZIFA_LSN_100_50_5_0_10_\\d+\\.RData$",
  full.names = TRUE
)

files4 <- list.files(
  path       = data_dir,
  pattern    = "^Output_ZIPPCA_LPNM_100_50_5_0_\\d+\\.RData$",
  full.names = TRUE
)

files6 <- list.files(
  path       = data_dir,
  pattern    = "^True_Compositions_100_50_5_0_\\d+\\.RData$",
  full.names = TRUE
)

idx1 <- as.integer(gsub(".*_(\\d+)\\.RData$", "\\1", files1))
files1 <- files1[order(idx1)]

idx2 <- as.integer(gsub(".*_(\\d+)\\.RData$", "\\1", files2))
files2 <- files2[order(idx2)]

idx3 <- as.integer(gsub(".*_(\\d+)\\.RData$", "\\1", files3))
files3 <- files3[order(idx3)]

idx4 <- as.integer(gsub(".*_(\\d+)\\.RData$", "\\1", files4))
files4 <- files4[order(idx4)]

idx6 <- as.integer(gsub(".*_(\\d+)\\.RData$", "\\1", files6))
files6 <- files6[order(idx6)]

zifa_list_100_50_5_0_0 <- lapply(files1, function(fp) {
  e <- new.env()
  load(fp, envir = e)
  objs <- ls(envir = e)
  if (length(objs) != 1) {
    stop("Expected exactly one object in ", basename(fp))
  }
  e[[objs]]
})

zifa_list_100_50_5_0_2 <- lapply(files2, function(fp) {
  e <- new.env()
  load(fp, envir = e)
  objs <- ls(envir = e)
  if (length(objs) != 1) {
    stop("Expected exactly one object in ", basename(fp))
  }
  e[[objs]]
})

zifa_list_100_50_5_0_10 <- lapply(files3, function(fp) {
  e <- new.env()
  load(fp, envir = e)
  objs <- ls(envir = e)
  if (length(objs) != 1) {
    stop("Expected exactly one object in ", basename(fp))
  }
  e[[objs]]
})


zippcalpnm_list_100_50_5_0 <- lapply(files4, function(fp) {
  e <- new.env()
  load(fp, envir = e)
  objs <- ls(envir = e)
  if (length(objs) != 1) {
    stop("Expected exactly one object in ", basename(fp))
  }
  e[[objs]]
})


compositions_list_100_50_5_0 <- lapply(files6, function(fp) {
  e <- new.env()
  load(fp, envir = e)
  objs <- ls(envir = e)
  if (length(objs) != 1) {
    stop("Expected exactly one object in ", basename(fp))
  }
  e[[objs]]
})


get(load("True_Factor_Loadings_100_50_5_0.RData"))
get(load("True_Latent_factors_100_50_5_0.RData"))
get(load("True_Kappa_100_50_5_0.RData"))
get(load("True_Beta_0_100_50_5_0.RData"))


meanforskewnormal_100_50_5_0_0 <- meanforskewnormal_100_50_5_0_2 <- list()
meanforskewnormal_100_50_5_0_10 <- list()

n_rep <- 250

zifa_100_50_5_0_0_time <- rep(0,n_rep)
zifa_100_50_5_0_2_time <- rep(0,n_rep)
zifa_100_50_5_0_10_time <- rep(0,n_rep)
zippcalpnm_100_50_5_0_time <- rep(0,n_rep)


for(i in 1:n_rep){
  meanforskewnormal_100_50_5_0_0[[i]] <- zifa_list_100_50_5_0_0[[i]]$Xi + 
    zifa_list_100_50_5_0_0[[i]]$Omega*(zifa_list_100_50_5_0_0[[i]]$Alpha/sqrt(1+zifa_list_100_50_5_0_0[[i]]$Alpha^2))*sqrt(2/pi)
  meanforskewnormal_100_50_5_0_2[[i]] <- zifa_list_100_50_5_0_2[[i]]$Xi + 
    zifa_list_100_50_5_0_2[[i]]$Omega*(zifa_list_100_50_5_0_2[[i]]$Alpha/sqrt(1+zifa_list_100_50_5_0_2[[i]]$Alpha^2))*sqrt(2/pi)
  meanforskewnormal_100_50_5_0_10[[i]] <- zifa_list_100_50_5_0_10[[i]]$Xi + 
    zifa_list_100_50_5_0_10[[i]]$Omega*(zifa_list_100_50_5_0_10[[i]]$Alpha/sqrt(1+zifa_list_100_50_5_0_10[[i]]$Alpha^2))*sqrt(2/pi)
  
  zifa_100_50_5_0_0_time[i] <- as.numeric(zifa_list_100_50_5_0_0[[i]]$Time_Taken_Mins, units = "mins")  
  zifa_100_50_5_0_2_time[i] <- as.numeric(zifa_list_100_50_5_0_2[[i]]$Time_Taken_Mins, units="mins")
  zifa_100_50_5_0_10_time[i] <- as.numeric(zifa_list_100_50_5_0_10[[i]]$Time_Taken_Mins, units="mins")
  zippcalpnm_100_50_5_0_time[i] <- as.numeric(zippcalpnm_list_100_50_5_0[[i]]$Time_Taken_Mins, units="mins")
}

corr_fl_zifa_100_50_5_0_0 <- corr_fs_zifa_100_50_5_0_0 <- rep(0, n_rep)
corr_fl_zifa_100_50_5_0_2 <- corr_fs_zifa_100_50_5_0_2 <- rep(0, n_rep)
corr_fl_zifa_100_50_5_0_10 <- corr_fs_zifa_100_50_5_0_10 <- rep(0, n_rep)

corr_fl_zippcalpnm_100_50_5_0 <- corr_fs_zippcalpnm_100_50_5_0 <- rep(0,n_rep)


for(i in 1:n_rep){
  corr_fl_zifa_100_50_5_0_0[i] <- cancor(beta, zifa_list_100_50_5_0_0[[i]]$R)$cor[5]
  corr_fl_zifa_100_50_5_0_2[i] <- cancor(beta, zifa_list_100_50_5_0_2[[i]]$R)$cor[5]
  corr_fl_zifa_100_50_5_0_10[i] <- cancor(beta, zifa_list_100_50_5_0_10[[i]]$R)$cor[5]
  
  corr_fl_zippcalpnm_100_50_5_0[i] <- cancor(rbind(beta,0),  zippcalpnm_list_100_50_5_0[[i]]$lvs$factor_coefs_j)$cor[5]
  
  corr_fs_zifa_100_50_5_0_0[i] <- cancor(f,  meanforskewnormal_100_50_5_0_0[[i]])$cor[5]
  corr_fs_zifa_100_50_5_0_2[i] <- cancor(f,  meanforskewnormal_100_50_5_0_2[[i]])$cor[5]
  corr_fs_zifa_100_50_5_0_10[i] <- cancor(f,  meanforskewnormal_100_50_5_0_10[[i]])$cor[5]
  
  corr_fs_zippcalpnm_100_50_5_0[i] <- cancor(f,  zippcalpnm_list_100_50_5_0[[i]]$lvs$factor_scores)$cor[5]
  
}

rmse_compositions_zifa_100_50_5_0_0 <- rep(0,n_rep)
rmse_compositions_zifa_100_50_5_0_2 <- rep(0,n_rep)
rmse_compositions_zifa_100_50_5_0_10 <- rep(0,n_rep)
rmse_compositions_zippcalpnm_100_50_5_0 <- rep(0,n_rep)

rmse_kappa_zifa_100_50_5_0_0 <- rep(0, n_rep)
rmse_kappa_zifa_100_50_5_0_2 <- rep(0, n_rep)
rmse_kappa_zifa_100_50_5_0_10 <- rep(0, n_rep)

rmse_kappa_zippcalpnm_100_50_5_0 <- rep(0,n_rep)

rmse_prod_zifa_100_50_5_0_0 <- rep(0,n_rep)
rmse_prod_zifa_100_50_5_0_2 <- rep(0,n_rep)
rmse_prod_zifa_100_50_5_0_10 <- rep(0,n_rep)
rmse_prod_zippcalpnm_100_50_5_0 <- rep(0,n_rep)


for(i in 1:n_rep){
  
  rmse_compositions_zifa_100_50_5_0_0[i] <- sqrt(mean((zifa_list_100_50_5_0_0[[i]]$Estimated_Compositions - compositions_list_100_50_5_0[[i]])^2))
  rmse_compositions_zifa_100_50_5_0_2[i] <- sqrt(mean((zifa_list_100_50_5_0_2[[i]]$Estimated_Compositions - compositions_list_100_50_5_0[[i]])^2))
  rmse_compositions_zifa_100_50_5_0_10[i] <- sqrt(mean((zifa_list_100_50_5_0_10[[i]]$Estimated_Compositions - compositions_list_100_50_5_0[[i]])^2))
  
  rmse_compositions_zippcalpnm_100_50_5_0[i] <- sqrt(mean((zippcalpnm_list_100_50_5_0[[i]]$Q - compositions_list_100_50_5_0[[i]])^2))
  

  rmse_kappa_zifa_100_50_5_0_0[i] <- sqrt(mean((kappa[-P] - (zifa_list_100_50_5_0_0[[i]]$Tau_1 / (zifa_list_100_50_5_0_0[[i]]$Tau_1+zifa_list_100_50_5_0_0[[i]]$Tau_2)))^2))
  rmse_kappa_zifa_100_50_5_0_2[i] <- sqrt(mean((kappa[-P] - (zifa_list_100_50_5_0_2[[i]]$Tau_1 / (zifa_list_100_50_5_0_2[[i]]$Tau_1+zifa_list_100_50_5_0_2[[i]]$Tau_2)))^2))
  rmse_kappa_zifa_100_50_5_0_10[i] <- sqrt(mean((kappa[-P] - (zifa_list_100_50_5_0_10[[i]]$Tau_1 / (zifa_list_100_50_5_0_10[[i]]$Tau_1+zifa_list_100_50_5_0_10[[i]]$Tau_2)))^2))
  
  rmse_kappa_zippcalpnm_100_50_5_0[i] <- sqrt(mean(((zippcalpnm_list_100_50_5_0[[i]]$lvs$gam[,1]/(zippcalpnm_list_100_50_5_0[[i]]$lvs$gam[,1]+zippcalpnm_list_100_50_5_0[[i]]$lvs$gam[,2])) - kappa)^2))
  

  rmse_prod_zifa_100_50_5_0_0[i] <- sqrt(mean((c( matrix(zifa_list_100_50_5_0_0[[i]]$A0,N, P - 1, byrow = TRUE) + tcrossprod(meanforskewnormal_100_50_5_0_0[[i]], zifa_list_100_50_5_0_0[[i]]$R)) - 
                                                 c( matrix(beta_0,N, P - 1, byrow = TRUE) + tcrossprod(f, beta)))^2)) 
  
  rmse_prod_zifa_100_50_5_0_2[i] <- sqrt(mean((c( matrix(zifa_list_100_50_5_0_2[[i]]$A0,N, P - 1, byrow = TRUE) + tcrossprod(meanforskewnormal_100_50_5_0_2[[i]], zifa_list_100_50_5_0_2[[i]]$R)) - 
                                                 c( matrix(beta_0,N, P - 1, byrow = TRUE) + tcrossprod(f, beta)))^2)) 
  rmse_prod_zifa_100_50_5_0_10[i] <- sqrt(mean((c( matrix(zifa_list_100_50_5_0_10[[i]]$A0,N, P - 1, byrow = TRUE) + tcrossprod(meanforskewnormal_100_50_5_0_10[[i]], zifa_list_100_50_5_0_10[[i]]$R)) - 
                                                  c( matrix(beta_0,N, P - 1, byrow = TRUE) + tcrossprod(f, beta)))^2)) 
  rmse_prod_zippcalpnm_100_50_5_0[i] <- sqrt(mean((c( matrix( zippcalpnm_list_100_50_5_0[[i]]$params$factor_coefs_0,N,P,byrow = TRUE) + tcrossprod(zippcalpnm_list_100_50_5_0[[i]]$lvs$factor_scores, zippcalpnm_list_100_50_5_0[[i]]$lvs$factor_coefs_j)) - 
                                                     c( matrix( c(beta_0,0),N, P , byrow = TRUE) + tcrossprod(f, rbind(beta,0))))^2))  
  
}