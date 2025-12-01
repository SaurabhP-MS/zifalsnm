library(dplyr)
library(tidyr)
library(ggplot2)

mk <- function(x, n, method, k, p) {
  tibble(
    n = n,
    method = method,
    k = paste0("k=", k),
    p = paste0("p=", p),
    RMSE = list(x)      
  )
}

# FIRSTLY, LOAD ALL RMSE'S VALUES


#############################################
#------------ For Factor Loadings-----------#
#############################################


df_FL_2 <-
  bind_rows(
    mk(rmse_factor_loadings_zifa_50_50_2, 50, "ZIFA-LSNM", 2, 50),
    mk(rmse_factor_loadings_zifa_100_50_2, 100, "ZIFA-LSNM", 2, 50),
    mk(rmse_factor_loadings_zifa_1000_50_2, 1000, "ZIFA-LSNM", 2, 50),
    mk(rmse_factor_loadings_zippca_50_50_2, 50, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_factor_loadings_zippca_100_50_2, 100, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_factor_loadings_zippca_1000_50_2, 1000, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_factor_loadings_zifa_50_100_2,   50,   "ZIFA-LSNM", 2, 100),
    mk(rmse_factor_loadings_zifa_100_100_2,  100,  "ZIFA-LSNM", 2, 100),
    mk(rmse_factor_loadings_zifa_1000_100_2, 1000, "ZIFA-LSNM", 2, 100),
    mk(rmse_factor_loadings_zippca_50_100_2,   50,   "ZIPPCA-LPNM", 2, 100),
    mk(rmse_factor_loadings_zippca_100_100_2,  100,  "ZIPPCA-LPNM", 2, 100),
    mk(rmse_factor_loadings_zippca_1000_100_2, 1000, "ZIPPCA-LPNM", 2, 100)
  ) %>%
  unnest(RMSE) %>%                                  
  mutate(
    n = factor(n, levels = c(50, 100, 1000)),
    method = factor(method, levels = c("ZIFA-LSNM", "ZIPPCA-LPNM")),
    p = factor(p, levels = c("p=50", "p=100"))
  )

df_FL_no_outliers_2 <- df_FL_2 %>%
  group_by(n, method, p, k) %>%
  mutate(
    Q1 = quantile(RMSE, 0.25, na.rm = TRUE),
    Q3 = quantile(RMSE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1
  ) %>%
  filter(RMSE >= (Q1 - 1.5 * IQR) & RMSE <= (Q3 + 1.5 * IQR)) %>%
  ungroup()

p1 <- ggplot(df_FL_no_outliers_2, aes(x = n, y = RMSE, fill = method)) +
  
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, alpha = 0.6, scale = "width") +
  
  geom_jitter(shape = 21, aes(color = method), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), size = 1.2, alpha = 0.9) +
  
  facet_grid(rows = vars(p), cols = vars(k)) +
  scale_fill_manual(values = c("ZIFA-LSNM" = "#2C7BE5", "ZIPPCA-LPNM" = "#D43F3A")) +
  scale_color_manual(values = c("ZIFA-LSNM" = "#003366", "ZIPPCA-LPNM" = "#800000")) +
  labs(
    title = "RMSE Comparison for Factor Loadings",
    x = "Sample Size", y = "RMSE", fill = "Model"
  ) +
  coord_cartesian(ylim = c(0, 1.5)) + 
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(
      face = "bold",      
      colour = "#000000", 
      size = 15           
    ),
    legend.position = "bottom",
    
    axis.text = element_text(colour = "#000000", face = "bold", size = 15),
    axis.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.text  = element_text(colour = "#000000", face = "bold", size = 13),
    
  ) +
  guides(color = "none") 

p1


###########################################
#------------ For Factor Scores-----------#
###########################################


df_FS_2 <-
  bind_rows(
    mk(rmse_factor_scores_zifa_50_50_2, 50, "ZIFA-LSNM", 2, 50),
    mk(rmse_factor_scores_zifa_100_50_2, 100, "ZIFA-LSNM", 2, 50),
    mk(rmse_factor_scores_zifa_1000_50_2, 1000, "ZIFA-LSNM", 2, 50),
    mk(rmse_factor_scores_zippca_50_50_2, 50, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_factor_scores_zippca_100_50_2, 100, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_factor_scores_zippca_1000_50_2, 1000, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_factor_scores_zifa_50_100_2,   50,   "ZIFA-LSNM", 2, 100),
    mk(rmse_factor_scores_zifa_100_100_2,  100,  "ZIFA-LSNM", 2, 100),
    mk(rmse_factor_scores_zifa_1000_100_2, 1000, "ZIFA-LSNM", 2, 100),
    mk(rmse_factor_scores_zippca_50_100_2,   50,   "ZIPPCA-LPNM", 2, 100),
    mk(rmse_factor_scores_zippca_100_100_2,  100,  "ZIPPCA-LPNM", 2, 100),
    mk(rmse_factor_scores_zippca_1000_100_2, 1000, "ZIPPCA-LPNM", 2, 100)
  ) %>%
  unnest(RMSE) %>%                                  
  mutate(
    n = factor(n, levels = c(50, 100, 1000)),
    method = factor(method, levels = c("ZIFA-LSNM", "ZIPPCA-LPNM")),
    p = factor(p, levels = c("p=50", "p=100"))
  )


df_FS_no_outliers_2 <- df_FS_2 %>%
  group_by(n, method, p, k) %>%
  mutate(
    Q1 = quantile(RMSE, 0.25, na.rm = TRUE),
    Q3 = quantile(RMSE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1
  ) %>%
  filter(RMSE >= (Q1 - 1.5 * IQR) & RMSE <= (Q3 + 1.5 * IQR)) %>%
  ungroup()

p2 <- ggplot(df_FS_no_outliers_2, aes(x = n, y = RMSE, fill = method)) +
  
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, alpha = 0.6, scale = "width") +
  
  geom_jitter(shape = 21, aes(color = method), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), size = 1.2, alpha = 0.9) +
  
  facet_grid(rows = vars(p), cols = vars(k)) +
  scale_fill_manual(values = c("ZIFA-LSNM" = "#2C7BE5", "ZIPPCA-LPNM" = "#D43F3A")) +
  scale_color_manual(values = c("ZIFA-LSNM" = "#003366", "ZIPPCA-LPNM" = "#800000")) + 
  labs(
    title = "RMSE Comparison for Factor Scores",
    x = "Sample size", y = "RMSE", fill = "Model"
  ) +
  coord_cartesian(ylim = c(0, 1.5)) + 
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(
      face = "bold",      
      colour = "#000000", 
      size = 15           
    ),
    legend.position = "bottom",
    
    axis.text = element_text(colour = "#000000", face = "bold", size = 15),
    axis.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.text  = element_text(colour = "#000000", face = "bold", size = 13),
  ) +
  guides(color = "none") 

p2


##########################################
#------------ For Compositions-----------#
##########################################


df_C_2 <-
  bind_rows(
    mk(rmse_compositions_zifa_50_50_2, 50, "ZIFA-LSNM", 2, 50),
    mk(rmse_compositions_zifa_100_50_2, 100, "ZIFA-LSNM", 2, 50),
    mk(rmse_compositions_zifa_1000_50_2, 1000, "ZIFA-LSNM", 2, 50),
    mk(rmse_compositions_zippca_50_50_2, 50, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_compositions_zippca_100_50_2, 100, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_compositions_zippca_1000_50_2, 1000, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_compositions_zifa_50_100_2,   50,   "ZIFA-LSNM", 2, 100),
    mk(rmse_compositions_zifa_100_100_2,  100,  "ZIFA-LSNM", 2, 100),
    mk(rmse_compositions_zifa_1000_100_2, 1000, "ZIFA-LSNM", 2, 100),
    mk(rmse_compositions_zippca_50_100_2,   50,   "ZIPPCA-LPNM", 2, 100),
    mk(rmse_compositions_zippca_100_100_2,  100,  "ZIPPCA-LPNM", 2, 100),
    mk(rmse_compositions_zippca_1000_100_2, 1000, "ZIPPCA-LPNM", 2, 100)
  ) %>%
  unnest(RMSE) %>%                                  
  mutate(
    n = factor(n, levels = c(50, 100, 1000)),
    method = factor(method, levels = c("ZIFA-LSNM", "ZIPPCA-LPNM")),
    p = factor(p, levels = c("p=50", "p=100"))
  )

df_C_no_outliers_2 <- df_C_2 %>%
  group_by(n, method, p, k) %>%
  mutate(
    Q1 = quantile(RMSE, 0.25, na.rm = TRUE),
    Q3 = quantile(RMSE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1
  ) %>%
  filter(RMSE >= (Q1 - 1.5 * IQR) & RMSE <= (Q3 + 1.5 * IQR)) %>%
  ungroup()


p3 <- ggplot(df_C_no_outliers_2, aes(x = n, y = RMSE, fill = method)) +
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, alpha = 0.6, scale = "width") +
  
  geom_jitter(shape = 21, aes(color = method), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), size = 1.2, alpha = 0.9) +
  
  facet_grid(rows = vars(p), cols = vars(k)) +
  scale_fill_manual(values = c("ZIFA-LSNM" = "#2C7BE5", "ZIPPCA-LPNM" = "#D43F3A")) +
  scale_color_manual(values = c("ZIFA-LSNM" = "#003366", "ZIPPCA-LPNM" = "#800000")) + 
  labs(
    title = "RMSE Comparison for Compositions",
    x = "Sample size", y = "RMSE", fill = "Model"
  ) +
  coord_cartesian(ylim = c(0, 0.04)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(
      face = "bold",      
      colour = "#000000", 
      size = 15         
    ),
    legend.position = "bottom",
    
    axis.text = element_text(colour = "#000000", face = "bold", size = 15),
    axis.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.text  = element_text(colour = "#000000", face = "bold", size = 13),
  ) +
  guides(color = "none") 

p3


####################################
#------------ For Kappas-----------#
####################################


df_Kappa_2 <-
  bind_rows(
    mk(rmse_kappa_zifa_50_50_2, 50, "ZIFA-LSNM", 2, 50),
    mk(rmse_kappa_zifa_100_50_2, 100, "ZIFA-LSNM", 2, 50),
    mk(rmse_kappa_zifa_1000_50_2, 1000, "ZIFA-LSNM", 2, 50),
    mk(rmse_kappa_zippca_50_50_2, 50, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_kappa_zippca_100_50_2, 100, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_kappa_zippca_1000_50_2, 1000, "ZIPPCA-LPNM", 2, 50),
    mk(rmse_kappa_zifa_50_100_2,   50,   "ZIFA-LSNM", 2, 100),
    mk(rmse_kappa_zifa_100_100_2,  100,  "ZIFA-LSNM", 2, 100),
    mk(rmse_kappa_zifa_1000_100_2, 1000, "ZIFA-LSNM", 2, 100),
    mk(rmse_kappa_zippca_50_100_2,   50,   "ZIPPCA-LPNM", 2, 100),
    mk(rmse_kappa_zippca_100_100_2,  100,  "ZIPPCA-LPNM", 2, 100),
    mk(rmse_kappa_zippca_1000_100_2, 1000, "ZIPPCA-LPNM", 2, 100)
  ) %>%
  unnest(RMSE) %>%                                  # now each value is a row
  mutate(
    n = factor(n, levels = c(50, 100, 1000)),
    method = factor(method, levels = c("ZIFA-LSNM", "ZIPPCA-LPNM")),
    p = factor(p, levels = c("p=50", "p=100"))
  )


df_Kappa_no_outliers_2 <- df_Kappa_2 %>%
  group_by(n, method, p, k) %>%
  mutate(
    Q1 = quantile(RMSE, 0.25, na.rm = TRUE),
    Q3 = quantile(RMSE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1
  ) %>%
  filter(RMSE >= (Q1 - 1.5 * IQR) & RMSE <= (Q3 + 1.5 * IQR)) %>%
  ungroup()

p4 <- ggplot(df_Kappa_no_outliers_2, aes(x = n, y = RMSE, fill = method)) +
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, alpha = 0.6, scale = "width") +
  
  geom_jitter(shape = 21, aes(color = method), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), size = 1.2, alpha = 0.9) +
  
  facet_grid(rows = vars(p), cols = vars(k)) +
  scale_fill_manual(values = c("ZIFA-LSNM" = "#2C7BE5", "ZIPPCA-LPNM" = "#D43F3A")) +
  scale_color_manual(values = c("ZIFA-LSNM" = "#003366", "ZIPPCA-LPNM" = "#800000")) +
  labs(
    title = "RMSE Comparison for Kappa",
    x = "Sample size", y = "RMSE", fill = "Model"
  ) +
  coord_cartesian(ylim = c(0, 0.09)) + 
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(
      face = "bold",      
      colour = "#000000", 
      size = 15         
    ),
    legend.position = "bottom",
    
    axis.text = element_text(colour = "#000000", face = "bold", size = 15),
    axis.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.title = element_text(colour = "#000000", face = "bold", size = 13),
    legend.text  = element_text(colour = "#000000", face = "bold", size = 13),
  ) +
  guides(color = "none") 

p4
