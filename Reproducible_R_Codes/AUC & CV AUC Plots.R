library(ggplot2)

k      <- 2:10
models <- c("ZIPPCA-LPNM",
            "ZIFA-LSNM, a0 = 2",  "ZIFA-LSNM, a0 = -2",
            "ZIFA-LSNM, a0 = 5",  "ZIFA-LSNM, a0 = -5")

dat <- data.frame(
  k     = rep(k, times = 5),
  model = factor(rep(models, each = length(k)), levels = models),
  
  AUC = c(
    0.6167, 0.5910, 0.6476, 0.6270, 0.7762, 0.8292, 0.8132, 0.8055, 0.7793,
    0.7782, 0.7752, 0.7926, 0.8158, 0.8281, 0.8492, 0.8405, 0.8266, 0.8564,
    0.7299, 0.7839, 0.7716, 0.8220, 0.8323, 0.8364, 0.8590, 0.8533, 0.8605,
    0.7731, 0.7782, 0.8009, 0.8215, 0.8256, 0.8806, 0.8317, 0.8230, 0.8251,
    0.7289, 0.7736, 0.7890, 0.8194, 0.8307, 0.8425, 0.8585, 0.8580, 0.8580),
  
  CV_AUC = c(
    0.5363, 0.4606, 0.5214, 0.4726, 0.6783, 0.7176, 0.7015, 0.6734, 0.6168,
    0.7335, 0.7306, 0.7377, 0.7466, 0.7492, 0.7601, 0.7464, 0.7071, 0.7227,
    0.6719, 0.7367, 0.7075, 0.7534, 0.7606, 0.7435, 0.7524, 0.7340, 0.7168,
    0.7343, 0.7307, 0.7509, 0.7530, 0.7558, 0.7982, 0.7377, 0.6968, 0.6963,
    0.6714, 0.7278, 0.7303, 0.7469, 0.7603, 0.7516, 0.7707, 0.7316, 0.7304)
)

cols   <- c("#525252", "#0072B2", "#D55E00", "#009E73", "#7B3294")
ltys   <- c("dotted",  "solid",   "dashed",  "solid",   "dashed")
shapes <- c(15, 16, 17, 18, 8)

legend_labels <- c("ZIPPCA-LPNM",
                   expression("ZIFA-LSNM," ~ alpha[0] ==  2),
                   expression("ZIFA-LSNM," ~ alpha[0] == -2),
                   expression("ZIFA-LSNM," ~ alpha[0] ==  5),
                   expression("ZIFA-LSNM," ~ alpha[0] == -5))


make_plot <- function(yvar, ylab) {
  ggplot(dat, aes(x = k, y = .data[[yvar]],
                  colour = model, linetype = model, shape = model)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.1) +
    scale_colour_manual(values = cols,   labels = legend_labels) +
    scale_linetype_manual(values = ltys, labels = legend_labels) +
    scale_shape_manual(values = shapes,  labels = legend_labels) +
    scale_x_continuous(breaks = k) +
    scale_y_continuous(limits = c(0.45, 0.90), breaks = seq(0.45, 0.90, 0.05)) +
    labs(x = "Number of latent factors, k", y = ylab) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title     = element_blank(),
      legend.position  = "bottom",
      legend.text      = element_text(size = 18),            
      legend.key.width = unit(2.2, "lines"),                 
      legend.key.height = unit(1.2, "lines"),
      axis.title.x = element_text(face = "bold", colour = "black", size = 18,
                                  margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", colour = "black", size = 18,
                                  margin = margin(r = 8)),
      axis.text.x      = element_text(face = "bold", colour = "black", size = 11),
      axis.text.y      = element_text(face = "bold", colour = "black", size = 11),
      axis.ticks       = element_line(colour = "black")      
    ) +
    guides(colour   = guide_legend(nrow = 2, byrow = TRUE,
                                   override.aes = list(size = 3, linewidth = 1)),
           linetype = guide_legend(nrow = 2, byrow = TRUE),
           shape    = guide_legend(nrow = 2, byrow = TRUE))
}

p_auc <- make_plot("AUC",    "AUC")
p_cv  <- make_plot("CV_AUC", "Cross-validated AUC")

p_auc
p_cv

 
# ggsave("fig_auc.png",    p_auc, width = 10, height = 6, dpi = 600)
# ggsave("fig_cv_auc.png", p_cv,  width = 10, height = 6, dpi = 600)
