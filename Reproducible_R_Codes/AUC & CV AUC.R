setwd("C:/Work Files/Simulation_Study_For_Proj_1/Jacob_Data_Analysis/Jacob_Data_Analysis/ZIFA LSNM")
load("Jacobs_ibd_2016")

varimax_R <- varimax(OUTPUT_ZIFA_CAVI_2$R)

meanforskewnormal <- OUTPUT_ZIFA_CAVI_2$Xi +
  OUTPUT_ZIFA_CAVI_2$Omega*(OUTPUT_ZIFA_CAVI_2$Alpha/sqrt(1+OUTPUT_ZIFA_CAVI_2$Alpha^2))*sqrt(2/pi)
# 
# ######
# 
# varimax_R <- varimax(OUTPUT_ZIFA_CAVI_5$R)
# 
# meanforskewnormal <- OUTPUT_ZIFA_CAVI_5$Xi +
#   OUTPUT_ZIFA_CAVI_5$Omega*(OUTPUT_ZIFA_CAVI_5$Alpha/sqrt(1+OUTPUT_ZIFA_CAVI_5$Alpha^2))*sqrt(2/pi)
# 
# 
# ######
# 
# varimax_R <- varimax(OUTPUT_ZIFA_CAVI_neg2$R)
# 
# meanforskewnormal <- OUTPUT_ZIFA_CAVI_neg2$Xi +
#   OUTPUT_ZIFA_CAVI_neg2$Omega*(OUTPUT_ZIFA_CAVI_neg2$Alpha/sqrt(1+OUTPUT_ZIFA_CAVI_neg2$Alpha^2))*sqrt(2/pi)
# 
# 
# ######
# 
# varimax_R <- varimax(OUTPUT_ZIFA_CAVI_neg5$R)
# 
# meanforskewnormal <- OUTPUT_ZIFA_CAVI_neg5$Xi +
#   OUTPUT_ZIFA_CAVI_neg5$Omega*(OUTPUT_ZIFA_CAVI_neg5$Alpha/sqrt(1+OUTPUT_ZIFA_CAVI_neg5$Alpha^2))*sqrt(2/pi)


aftervarimax_mean_sn <- meanforskewnormal %*% varimax_R$rotmat

### AUC ###

y <- ifelse(metadata$Study.Group %in% c("CD", "UC"), 1, 0)

dat <- data.frame(y = y, aftervarimax_mean_sn)              

fit  <- glm(y ~ ., data = dat, family = binomial)
prob <- predict(fit, type = "response")

auc <- function(y, prob) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r  <- rank(prob)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc(dat$y, prob)

### CV-AUC ###

cv_auc <- function(dat, seed) {
  set.seed(seed)
  folds <- sample(rep(1:5, length.out = nrow(dat))) 
  
  prob_cv <- numeric(nrow(dat))
  for (f in 1:5) {
    m <- glm(y ~ ., data = dat[folds != f, ], family = binomial)
    prob_cv[folds == f] <- predict(m, newdata = dat[folds == f, ],
                                   type = "response")
  }
  auc(dat$y, prob_cv)
}

many <- sapply(1:100, function(s) cv_auc(dat, seed = s))
mean(many)              
