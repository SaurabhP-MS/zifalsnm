library(sn)
library(moments)
library(ggplot2)

n      <- 100
p      <- 200
k      <- 5
R      <- 100
alphas <- c(-10, -5, 0, 5, 10)

beta   <- matrix(1, nrow = p - 1, ncol = k, byrow = TRUE)
kappa  <- runif(p, 0, 0.85)
beta_0 <- rep(0, p-1)

B0      <- matrix(beta_0, n, p - 1, byrow = TRUE)   
kap_vec <- rep(kappa, each = n)                    
taxa    <- seq_len(p - 1)


N_rows  <- length(alphas) * R * (p - 1)
res_alpha <- numeric(N_rows)
res_rep   <- integer(N_rows)
res_taxon <- integer(N_rows)
res_skew  <- numeric(N_rows)
pos <- 0L


cs <- c(1, 0.5, 1e-2, 1e-4)

alr_skew <- function(C, cs) {
  ref <- ncol(C)
  js  <- seq_len(ncol(C) - 1L)
  
  out <- list()
  for (cc in cs) {
    s <- sapply(js, function(j) skewness(log((C[, j] + cc) / (C[, ref] + cc))))
    out[[length(out) + 1L]] <- data.frame(c = cc, taxon = js, skewness = s)
  }
  do.call(rbind, out)
}

res <- list(); m <- 0L

for (alpha in alphas) {
  
  delta <- alpha / sqrt(1 + alpha^2)
  xi    <- -1 * delta * sqrt(2 / pi)
  
  for (r in seq_len(R)) {
    
    f <- if (alpha == 0) {
      matrix(rnorm(n * k), nrow = n, ncol = k)
    } else {
      matrix(rsn(n * k, xi = xi, omega = 1, alpha = alpha),
             nrow = n, ncol = k)
    }
    
    z       <- matrix(rbinom(n * p, size = 1, prob = kap_vec), nrow = n, ncol = p)
    num_mat <- (1 - z[, taxa, drop = FALSE]) * exp(B0 + tcrossprod(f, beta))
    denom   <- 1 + rowSums(num_mat)
    rho     <- cbind(num_mat / denom, 1 / denom)
    COUNT_MATRIX <- matrix(0, n, p)
    for (i in 1:n) {
       COUNT_MATRIX[i, ] <- rmultinom(1, as.integer(runif(1, 8e4, 1e5)),
                                      prob = rho[i, ])
      
    }
    
    m <- m + 1L
    res[[m]] <- cbind(alpha = alpha, rep = r, alr_skew(COUNT_MATRIX, cs))
  }
}

results <- do.call(rbind, res)
results <- results[is.finite(results$skewness), ]

results$alpha_f <- factor(paste0("alpha=", results$alpha),
                          levels = paste0("alpha=", alphas))
results$c_f <- factor(paste0("c = ", results$c),
                      levels = paste0("c = ", cs))


env <- do.call(rbind, lapply(cs, function(cc) {
  g <- results$skewness[results$c == cc & results$alpha == 0]
  data.frame(c    = cc,
             c_f  = factor(paste0("c = ", cc), levels = paste0("c = ", cs)),
             lo95 = quantile(g, 0.025),
             hi95 = quantile(g, 0.975))
}))

fill_cols <- c("alpha=-10" = "#2B4C6F", "alpha=-5" = "#8CB4D8",
               "alpha=0"   = "#BEBEBE",
               "alpha=5"   = "#E07060", "alpha=10" = "#8B2020")

border_cols <- c("alpha=-10" = "#1A3050", "alpha=-5" = "#5A8AB0",
                 "alpha=0"   = "#888888",
                 "alpha=5"   = "#C04030", "alpha=10" = "#601010")

make_plot <- function(cc) {
  d <- results[results$c == cc, ]
  e <- env[env$c == cc, ]
  
  ggplot(d, aes(alpha_f, skewness, fill = alpha_f, colour = alpha_f)) +
    geom_rect(data = e, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = lo95, ymax = hi95),
              fill = "grey60", alpha = 0.20) +
    geom_violin(trim = FALSE, alpha = 0.45, linewidth = 0.4) +
    geom_boxplot(width = 0.18, outlier.size = 0.8, outlier.alpha = 0.4,
                 alpha = 0.7, linewidth = 0.4) +
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = border_cols) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(title    = paste0("Per-taxon ALR skewness, pseudo-count c = ", cc),
         subtitle = NULL,
         x = NULL, y = "Per-taxon ALR skewness") +
    theme_minimal(base_size = 12) +
    theme(legend.position     = "none",
          plot.title          = element_text(hjust = 0.5, size = 16,          
                                             face = "bold", colour = "black"),
          plot.title.position = "plot",
          panel.grid.major.x  = element_blank(),
          panel.grid.minor    = element_blank(),
          axis.title.x        = element_text(size = 14, face = "bold", colour = "black"),
          axis.title.y        = element_text(size = 14, face = "bold", colour = "black",
                                             margin = margin(r = 8)),
          axis.text.x         = element_text(angle = 45, hjust = 1, size = 12,  
                                             face = "bold", colour = "black"),
          axis.text.y         = element_text(size = 12, face = "bold", colour = "black"),
          axis.ticks          = element_line(colour = "black"),
          panel.border        = element_rect(colour = "black", fill = NA,
                                             linewidth = 0.8))
}

plots <- lapply(cs, make_plot)
names(plots) <- paste0("c=", cs)

for (p in plots) print(p) 

# ggsave("plotc=1.png", plots$`c=1`, width = 9, height = 6, dpi = 600)
# ggsave("plotc=0.5.png", plots$`c=0.5`, width = 9, height = 6, dpi = 600)
# ggsave("plotc=1e-2.png", plots$`c=0.01`, width = 9, height = 6, dpi = 600)
# ggsave("plotc=1e-4.png", plots$`c=1e-04`, width = 9, height = 6, dpi = 600)


