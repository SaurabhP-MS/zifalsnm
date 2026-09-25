# Load Output for K=6 
# Load meta data

O <- OUTPUT_ZIFA_CAVI_neg2

meanforskewnormal <- O$Xi +
  O$Omega * (O$Alpha / sqrt(1 + O$Alpha^2)) * sqrt(2/pi)
Fscore <- meanforskewnormal

grp <- metadata$Study.Group
ibd <- ifelse(grp == "Normal", "Normal", "IBD")

separation <- function(X1, X2) {
  
  n1 <- nrow(X1); n2 <- nrow(X2); p <- ncol(X1); N <- n1 + n2
  
  
  S_pooled <- ((n1 - 1) * cov(X1) + (n2 - 1) * cov(X2)) / (N - 2)
  D2 <- mahalanobis(colMeans(X1), colMeans(X2), S_pooled)
  
  
  T2 <- (n1 * n2 / N) * D2
  
  
  Fstat <- (N - p - 1) / (p * (N - 2)) * T2
  df1 <- p; df2 <- N - p - 1
  
  
  pval <- pf(Fstat, df1, df2, lower.tail = FALSE)
  
  
   D2u <- ((N - p - 3) / (N - 2)) * D2 - p * N / (n1 * n2)
  
  c(n1 = n1, n2 = n2, D2 = D2, D = sqrt(max(D2u, 0)),
    F = Fstat, p = pval)
}

results <- rbind(
  "CD vs Normal"          = separation(Fscore[grp == "Normal", ], Fscore[grp == "CD", ]),
  "IBD (CD+UC) vs Normal" = separation(Fscore[ibd == "Normal", ], Fscore[ibd == "IBD", ]),
  "UC vs Normal"          = separation(Fscore[grp == "Normal", ], Fscore[grp == "UC", ])
)

print(round(results, 4))

signif(results[, "p"], 3)
