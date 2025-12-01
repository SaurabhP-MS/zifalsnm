


###############################################################################################
########## Multi-Step Preprocessing and Filtering Process Of Publicly Available Data ##########
###############################################################################################



load("Jacobs_ibd_2016")

df <- genera.counts
mat <- as.matrix(df)
mat <- mat[ , -1]
mat <- apply(mat, 2, as.numeric)

# Remove columns with unknown genera
filtered_mat <- mat[, !grepl("g__$", colnames(mat))] 
filtered_mat <- filtered_mat[, colnames(filtered_mat) != "Unclassified"]

taxa_names <- colnames(filtered_mat)

# Extract genus names (everything after "g__")
genus_names <- sub(".*g__", "", taxa_names)

# Define function to identify valid genus names
is_valid_genus <- function(name) {
  grepl("^[A-Z][a-z]+$", name)
}

# Logical vector of valid genus names
keep_cols <- sapply(genus_names, is_valid_genus)

# Subset the matrix and rename columns
filtered_mat <- filtered_mat[, keep_cols]

colnames(filtered_mat) <- genus_names[keep_cols]

# Output dimensions and sample column names
dim(filtered_mat)
head(colnames(filtered_mat))

# Remove zeros with more than 97% of zeros

# Calculate proportion of zeros in each column
zero_proportion <- colMeans(filtered_mat == 0)

# Keep only columns where the proportion of zeros is <= 97%
filtered_data <- filtered_mat[, zero_proportion <= 0.97]


# Prevalence Filtering

# Keep taxa present in at least 5% of samples
prevalence_thresh <- 0.05 * nrow(filtered_data)
prevalent_taxa <- colSums(filtered_data > 0) >= prevalence_thresh
mat_filt <- filtered_data[, prevalent_taxa]

cat("After prevalence filtering:", dim(mat_filt), "\n")


# Abundance Filtering

# Keep taxa with total count across all samples ≥ 20
abundance_thresh <- 20
abundant_taxa <- colSums(mat_filt) >= abundance_thresh
mat_filt <- mat_filt[, abundant_taxa]

cat("After abundance filtering:", dim(mat_filt), "\n")



##############################################
########## Reproduce Skewness Graph ##########
##############################################



Y <- mat_filt + 1e-8
install.packages("compositions")
library(compositions)
d <- alr(Y)

install.packages("moments")
library(moments)
skews <- apply(d, 2, skewness)

skew_vals <- apply(d, 2, skewness)
df_skew  <- data.frame(
  Taxon    = names(skew_vals),
  Skewness = skew_vals
) %>%
  # sort by absolute skew so the most skewed appear at the top
  arrange(desc(abs(Skewness))) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon))


p <- ggplot(df_skew, aes(x = Taxon, y = Skewness, fill = Skewness)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(option = "D", guide = guide_colorbar(title="Skew")) +
  coord_flip() +
  
  labs(
    title = "Distribution of Skewness Across Taxa",
    x     = "Taxa (Sorted by Skewness)",
    y     = "Skewness"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face="bold", hjust=0.5, size=16),
    axis.text.y     = element_blank(),
    axis.text.x     = element_text(colour = "#000000", face = "bold", size = 11),
    axis.title.x    = element_text(colour = "#000000", face = "bold", size = 13),
    axis.title.y    = element_text(colour = "#000000", face = "bold", size = 13),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    
    panel.grid.major.y = element_line(colour = "grey98"), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.position = "right"
  )
p



#########################################################
########## Reproduce Real Data Analysis Values ##########
#########################################################



# Install my package from GitHub
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("SaurabhP-MS/zifalsnm", dependencies=TRUE)
library(zifalsnm)

OUTPUT_ZIFA_LSNM <- ZIFA_LSNM(X = mat_filt, num_fac = 3)

rownames(OUTPUT_ZIFA_LSNM$R) <- colnames(mat_filt)

a2 <- varimax(OUTPUT_ZIFA_LSNM$R)

rownames(a2$loadings) <- colnames(mat_filt)

meanforskewnormal <- OUTPUT_ZIFA_LSNM$Xi + OUTPUT_ZIFA_LSNM$Omega*(OUTPUT_ZIFA_LSNM$Alpha/sqrt(1+OUTPUT_ZIFA_LSNM$Alpha^2))*sqrt(2/pi)

aftervarimax_mean_sn <- meanforskewnormal %*% a2$rotmat


install.packages("tidyverse")
library(tidyverse)



######################################################
########## Reproduce Factor Loadings Graphs ##########
######################################################



load1 <- tibble(
  Taxon   = if (!is.null(rownames(a2$loadings))) rownames(a2$loadings) else paste0("Taxon_", seq_len(nrow(a2$loadings))),
  Loading = as.numeric(a2$loadings[,1])
)

load2 <- tibble(
  Taxon   = if (!is.null(rownames(a2$loadings))) rownames(a2$loadings) else paste0("Taxon_", seq_len(nrow(a2$loadings))),
  Loading = as.numeric(a2$loadings[, 2])
)

load3 <- tibble(
  Taxon   = if (!is.null(rownames(a2$loadings))) rownames(a2$loadings) else paste0("Taxon_", seq_len(nrow(a2$loadings))),
  Loading = as.numeric(a2$loadings[, 3])
)

# Top loading taxa 
topN <- 25  
top_abs_1 <- load1 %>%
  mutate(abs_loading = abs(Loading), sign = ifelse(Loading >= 0, "Positive", "Negative")) %>%
  slice_max(abs_loading, n = topN) %>%
  arrange(Loading) %>%
  mutate(Taxon = fct_inorder(Taxon))

top_abs_2 <- load2 %>%
  mutate(abs_loading = abs(Loading), sign = ifelse(Loading >= 0, "Positive", "Negative")) %>%
  slice_max(abs_loading, n = topN) %>%
  arrange(Loading) %>%
  mutate(Taxon = fct_inorder(Taxon))

top_abs_3 <- load3 %>%
  mutate(abs_loading = abs(Loading), sign = ifelse(Loading >= 0, "Positive", "Negative")) %>%
  slice_max(abs_loading, n = topN) %>%
  arrange(Loading) %>%
  mutate(Taxon = fct_inorder(Taxon))

dark_colors <- c(
  "Negative" = "#B2182B",  
  "Positive" = "#01665E"
)

p_top_1 <- ggplot(top_abs_1, aes(x = Taxon, y = Loading, color = sign)) +
  geom_segment(aes(xend = Taxon, y = 0, yend = Loading), linewidth = 0.6) +
  geom_point(size = 2) +
  
  scale_color_manual(values = dark_colors) + 
  
  coord_flip() +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour = "black", face = "bold", size=14), 
    axis.text.x = element_text(colour = "black", face = "bold", size=14)
  )
p_top_1

p_top_2 <- ggplot(top_abs_2, aes(x = Taxon, y = Loading, color = sign)) +
  geom_segment(aes(xend = Taxon, y = 0, yend = Loading), linewidth = 0.6) +
  geom_point(size = 2) +
  
  scale_color_manual(values = dark_colors) + 
  
  coord_flip() +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour = "black", face = "bold", size=14), 
    axis.text.x = element_text(colour = "black", face = "bold", size=14)
  )
p_top_2


p_top_3 <- ggplot(top_abs_3, aes(x = Taxon, y = Loading, color = sign)) +
  geom_segment(aes(xend = Taxon, y = 0, yend = Loading), linewidth = 0.6) +
  geom_point(size = 2) +
  
  scale_color_manual(values = dark_colors) + 
  
  coord_flip() +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour = "black", face = "bold", size=14), 
    axis.text.x = element_text(colour = "black", face = "bold", size=14)
  )
p_top_3



######################################################
########## Reproduce Factor Score Graphs ##########
######################################################



disease_colors <- c(
  "Normal" = "#E64B35",  
  "UC"     = "#4DBBD5",  
  "CD"     = "#00A087"  
)

df <- as.data.frame(aftervarimax_mean_sn_2)
if (is.null(colnames(df))) colnames(df) <- paste0("F", seq_len(ncol(df)))
df$Disease <- factor(metadata$Study.Group, levels = c("Normal","UC","CD"))

# All unique column pairs
pairs_idx <- combn(seq_len(ncol(df)-1), 2, simplify = FALSE)  

plot_df <- map_dfr(
  pairs_idx,
  ~{
    xj <- .x[1]; yj <- .x[2]
    tibble(
      x       = df[[xj]],
      y       = df[[yj]],
      x_name  = names(df)[xj],
      y_name  = names(df)[yj],
      facet   = paste0(names(df)[yj], " vs ", names(df)[xj]),
      Disease = df$Disease
    )
  }
)

p_facets <- ggplot(plot_df, aes(x = x, y = y, colour = Disease)) +
  geom_point(size = 1.6) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.4, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.35, linetype = 2, show.legend = FALSE) +
  scale_color_manual(values = disease_colors) +
  facet_wrap(~ facet, scales = "free", ncol = 5) +
  labs(x = NULL, y = NULL, colour = "Disease") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.grid.minor = element_blank()
  )

p_facets
