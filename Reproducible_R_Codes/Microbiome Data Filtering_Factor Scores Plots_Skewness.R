setwd("C:/Work Files/Real Data Analysis/Microbiome data sets")
load("Jacobs_ibd_2016")

df <- genera.counts
mat <- as.matrix(df)
mat <- mat[ , -1]
mat <- apply(mat, 2, as.numeric)

# Remove columns with unknown genera
filtered_mat <- mat[, !grepl("g__$", colnames(mat))] 
filtered_mat <- filtered_mat[, colnames(filtered_mat) != "Unclassified"]     

# Assume your matrix is named 'filtered_mat'
taxa_names <- colnames(filtered_mat)

# Step 1: Extract genus names (everything after "g__")
genus_names <- sub(".*g__", "", taxa_names)

# Step 2: Define function to identify valid genus names
# Criteria:
# - Must start with a capital letter
# - Must NOT contain digits, underscores, or hyphens
# - Must only contain letters (i.e., [A-Za-z]+)

is_valid_genus <- function(name) {
  grepl("^[A-Z][a-z]+$", name)
}

# Logical vector of valid genus names
keep_cols <- sapply(genus_names, is_valid_genus)

# Step 3: Subset the matrix and rename columns
filtered_mat <- filtered_mat[, keep_cols]

colnames(filtered_mat) <- genus_names[keep_cols]

# Output dimensions and sample column names
dim(filtered_mat)
head(colnames(filtered_mat))


# Remove zeros with more than 85% of zeros

# Calculate proportion of zeros in each column
zero_proportion <- colMeans(filtered_mat == 0)

# Keep only columns where the proportion of zeros is <= 85%
filtered_data <- filtered_mat[, zero_proportion <= 0.85]

# Keep taxa with total count across all samples ≥ 20
abundance_thresh <- 20
abundant_taxa <- colSums(filtered_data) >= abundance_thresh
mat_filt <- filtered_data[, abundant_taxa]

cat("After abundance filtering:", dim(mat_filt), "\n")

# Have the last column of the matrix least abundant
mat_filt[, c(122, 132)] <- mat_filt[, c(132, 122)]
colnames(mat_filt)[c(122, 132)] <- colnames(mat_filt)[c(132, 122)]

########### Real Data analysis figures ################

rownames(OUTPUT$R) <- colnames(mat_filt[,-132])

a1 <- varimax(OUTPUT$R)

library(tidyverse)

# --- Inputs ---
# F_scores : n x K matrix/data.frame of factor scores (no rownames)
# metadata  : data.frame with columns at least: Sample, Study.Group
#            Study.Group levels should be CD, UC, Normal (any order; we set below)

# 0) Attach Sample IDs to F_scores if missing
attach_sample_ids <- function(F_scores, meta_df,
                              candidates = c("Sample","sample","sample_id","SampleID","ID")){
  FS <- as.data.frame(F_scores)
  # If F_scores already has a sample ID column, normalize its name to 'Sample'
  have <- intersect(names(FS), candidates)
  if (length(have) >= 1) {
    FS <- FS %>% rename(Sample = !!have[1]) %>% mutate(Sample = as.character(Sample))
  } else {
    # Otherwise assume row-order matches meta_df; attach meta_df$Sample
    stopifnot(nrow(FS) == nrow(meta_df))
    FS <- FS %>% mutate(Sample = as.character(meta_df$Sample))
  }
  FS
}

meanforskewnormal1 <- OUTPUT$Xi + OUTPUT$Omega*(OUTPUT$Alpha/sqrt(1+OUTPUT$Alpha^2))*sqrt(2/pi)

aftervarimax_mean_sn_1 <- meanforskewnormal1 %*% a1$rotmat

FS <- attach_sample_ids(aftervarimax_mean_sn_1, metadata)

rownames(a1$loadings) <- colnames(mat_filt[,-132])


disease_colors <- c(
  "Normal" = "#E64B35",  # Orange
  "UC"     = "#4DBBD5",  # Sky Blue
  "CD"     = "#00A087"   # Bluish Green
)

df <- as.data.frame(aftervarimax_mean_sn_1)
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

# Faceted scatterplots with 68% normal ellipses
p_facets <- ggplot(plot_df, aes(x = x, y = y, colour = Disease)) +
  geom_point(size = 1.6) +
  # Ellipses help visualize cluster separation per disease
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.4, show.legend = FALSE) +
  # Optional trend line per facet (remove if not needed)
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
ggsave("fig_factor_scores.png", p_facets, width = 10, height = 6, dpi = 600)

###### Skewness Figures in Section 8 Supplementary Material ########

Y <- mat_filt + 0.5
library(compositions)
library(tidyverse)
d <- alr(Y)

library(moments)
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
    
    # --- MODIFIED GRID LINES ---
    panel.grid.major.y = element_line(colour = "grey98"), # Light grey horizontal grid
    panel.grid.major.x = element_blank(), # No vertical grid
    panel.grid.minor = element_blank(), # No minor grid
    # --- END MODIFIED SECTION ---
    
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.position = "right"
  )
p

sum(abs(skew_vals) > 0.5) / length(skew_vals)

