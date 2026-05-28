############################################################################################################################
#### CODE USED TO GENERATE FIGURE 3
############################################################################################################################

#A

# BETADIV 

# BRAY CURTIS

library(vegan)
data1<-read.delim("data/bacteriome.txt", row.names = 1, check.names = F)
dataTransposed1<-t(data1)
dis <- vegdist(dataTransposed1, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"WeightedBrayCurtis.txt", sep = '\t')


# Reformat Betadiv superimposed Output table (BRAY CURTIS)

library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)

data <- read.delim("WeightedBrayCurtis.txt", stringsAsFactors = FALSE, row.names = 1, header = TRUE, check.names = FALSE)
meta <- read.delim("metadata.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

superimposed_matrix <- as.matrix(data)
upper_logical <- upper.tri(superimposed_matrix)
upper_triangle_subset <- matrix(NA, nrow = nrow(superimposed_matrix), ncol = ncol(superimposed_matrix))
upper_triangle_subset[upper_logical] <- superimposed_matrix[upper_logical]
row_names <- row.names(data)
col_names <- colnames(data)
rownames(upper_triangle_subset) <- row_names
colnames(upper_triangle_subset) <- col_names
print(upper_triangle_subset)
upper_triangle_df <- as.data.frame(upper_triangle_subset) %>%
  rownames_to_column(var = "sample1")

melted_data <- upper_triangle_df %>%
  pivot_longer(
    cols = -sample1,
    names_to = "sample2",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

cleaned_data <- na.omit(melted_data)
meta_sample1 <- meta
names(meta_sample1)[names(meta_sample1) == "Code"] <- "sample1"
merged_data <- merge(cleaned_data, meta_sample1[, c("sample1", "cluster", "Day")],
                     by = "sample1", all.x = TRUE)
names(merged_data)[names(merged_data) == "cluster"] <- "cluster"
names(merged_data)[names(merged_data) == "Day"] <- "Day_sample1"
meta_sample2 <- meta
names(meta_sample2)[names(meta_sample2) == "Code"] <- "sample2"
merged_data <- merge(merged_data, meta_sample2[, c("sample2", "cluster", "Day")],
                     by = "sample2", all.x = TRUE)
names(merged_data)[names(merged_data) == "cluster"] <- "cluster2"
names(merged_data)[names(merged_data) == "Day"] <- "Day_sample2"
head(merged_data)
merged_data <- merged_data[, c("sample1", "Day_sample1", "cluster1", "sample2","Day_sample2","cluster2","value")]
write.table(merged_data, file="wbc_formatted.txt", sep="\t", quote=FALSE, row.names=FALSE)


### boxplot 

#### boxplot WILCOX NORMAL

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(dplyr)

dat <- read.delim("wbc_formatted.txt", stringsAsFactors = FALSE)

head(dat)

custom_colors <- c("C1_high_risk" = "indianred1", "C2_low_risk" = "limegreen")


dat <- dat %>%
  mutate(across(
    .cols = where(is.numeric) & !all_of(c("cluster")),
    .fns = ~log10(.x + 1)
  ))

dat %>% sample_n_by(cluster, size = 2)
dat %>%
  group_by(cluster) %>%
  get_summary_stats(Betadiv, type = "median_iqr")

stat.test <- dat %>% 
  wilcox_test(Betadiv ~ cluster) %>%
  add_significance()
stat.test
dat %>% wilcox_effsize(Betadiv ~ cluster)
stat.test <- stat.test %>% add_xy_position(x = "cluster")

p <-ggplot(dat, aes(cluster, Betadiv)) +  
  # Boxplot layer: will be drawn first, in the background
  geom_boxplot(aes(fill = cluster), width = 0.2, color = "black", outlier.shape = NA, alpha = 0.2) +  
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) +  
  
  # Add p-value annotation with stat_pvalue_manual()
  stat_pvalue_manual(stat.test, tip.length = 0, size = 7, bracket.size = 2) +
  
  # Customize labels and theme
  labs(subtitle = get_test_label(stat.test, detailed = TRUE)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "none"
  ) +
  
  # Jittered points layer: will be drawn last, in the foreground
  geom_jitter(aes(color = cluster), width = 0.03, size = 10, alpha = 0.1) 
p


# === 1. Plot boxplots with same format as Avg_RPKM script ===
p <- ggplot(dat, aes(x = cluster, y = Betadiv, fill = cluster)) + 
  
  # Boxplot layer
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", size = 2) +
  
  # Jittered points
  geom_jitter(shape = 21, color = "black", aes(fill = cluster),
              width = 0.1, size = 10, stroke = 1, alpha = 0.5) +
  
  # Statistical comparison (Wilcoxon)
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("C1_high_risk","C2_low_risk")),
                     label.y = max(dat$Betadiv) * 1.05,
                     size = 10, bracket.size = 2) +
  
  # Fill colors
  scale_fill_manual(values = custom_colors) +
  
  # Theme customization
  theme_pubclean() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 15, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold", color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    legend.position = "none"
  ) +
  
  # Axis labels
  labs(y = "Betadiv", x = "")

# === 2. Display plot ===
print(p)

# === 3. Wilcoxon test ===
wilcox_result <- wilcox.test(Betadiv ~ cluster, data = dat)
print(wilcox_result)

pdf("Betadiv_boxplot.pdf",width=5,height=5);
p
dev.off()


#B

#############################PCOA#####################################

library(ggplot2)
library(vegan)
library(dplyr)

set.seed(1)
data1 <- read.delim("data/bacteriome.txt", row.names = 1)
dataTransposed1 <- t(data1)

metadata <- read.delim("data/metadata.txt", row.names = 1)

# ---- CORRECT ALIGNMENT (this is the fix) ----
# Keep only samples present in BOTH objects
common_samples <- intersect(rownames(dataTransposed1), rownames(metadata))

# Subset both objects to the same samples
dataTransposed1 <- dataTransposed1[common_samples, ]
metadata <- metadata[common_samples, ]

# Check dimensions (should match)
dim(dataTransposed1)
dim(metadata)

# Compute Bray-Curtis distance
dist.1 <- vegdist(dataTransposed1, method = "bray")

# Run adonis
adonis.test1.cont.subjID <- adonis2(dist.1 ~ GROUP, data = metadata, permutations = 9999)
adonis.test1.cont.subjID

pval <- adonis.test1.cont.subjID$`Pr(>F)`[1]
R2   <- adonis.test1.cont.subjID$R2[1]

pval_text <- paste0("p = ", signif(pval, 3))
R2_text   <- paste0("R² = ", signif(R2, 3))

# Betadisper
dispersion <- betadisper(dist.1, group = metadata$GROUP)
scores_d <- as.data.frame(scores(dispersion, display = "sites"))

scores_d$cluster <- metadata$GROUP

colors <- c("C1" = "indianred1", "C2" = "limegreen")
#colors <- c("G" = "gold", "S" = "purple")

# % variance explained
eig <- dispersion$eig / sum(dispersion$eig)
xlab1 <- paste0("PCoA1 (", round(eig[1] * 100, 1), "%)")
ylab1 <- paste0("PCoA2 (", round(eig[2] * 100, 1), "%)")

# Compute centroids
centroids <- scores_d %>%
  group_by(cluster) %>%
  summarise(
    centroid_x = mean(PCoA1),
    centroid_y = mean(PCoA2)
  )

# Join centroid coordinates back to each point
scores_seg <- scores_d %>%
  left_join(centroids, by = "cluster")
# PLOT
pdf("pcoa_bac.pdf", width = 6, height = 4)

p <- ggplot(scores_seg, aes(x = PCoA1, y = PCoA2, color = cluster)) +
  
  # Lines from points to centroid
  geom_segment(aes(
    xend = centroid_x,
    yend = centroid_y
  ),
  alpha = 0.5,
  linewidth = 2
  ) +
  
  # Points
  geom_point(size = 4, alpha=0.5) +
  
  # Centroids
  geom_point(
    data = centroids,
    aes(x = centroid_x, y = centroid_y, color = cluster),
    size = 6,
    shape = 4,
    stroke = 2
  ) +
  
  scale_color_manual(values = colors) +
  
  theme_classic(base_size = 18) +
  
  labs(
    title = "Bray-Curtis PCoA",
    x = xlab1,
    y = ylab1
  ) +
  
  annotate(
    "text",
    x = min(scores_d$PCoA1),
    y = max(scores_d$PCoA2),
    label = pval_text,
    hjust = 0.3,
    vjust = 1,
    size = 6
  ) +
  
  annotate(
    "text",
    x = min(scores_d$PCoA1),
    y = max(scores_d$PCoA2) * 0.92,
    label = R2_text,
    hjust = 0.3,
    vjust = 2,
    size = 6
  )

print(p)

dev.off()

# C

### Donut plot 

library(dplyr)
library(ggplot2)
library(reshape2)

# Load your abundance matrix
matrix <- read.delim("datat/bacteriome.txt", header = TRUE, sep = "\t", check.names = F)

# Load metadata (Sample -> Cluster)
metadata <- read.delim("data/metadata.txt", header = TRUE, sep = "\t", check.names = F)

# Melt the abundance matrix
melted_data <- reshape2::melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Genus", "SAMPLEID", "Value")

# Merge with metadata to get cluster info
melted_data <- melted_data %>%
  left_join(metadata[, c("SAMPLEID", "cluster")], by = "SAMPLEID")

# Calculate median relative abundance per Genus per cluster
median_cluster <- melted_data %>%
  group_by(cluster, Genus) %>%
  summarise(MeanAbundance  = median(Value, na.rm = TRUE), .groups = "drop")

# Normalize per cluster to get percentages
median_cluster <- median_cluster %>%
  group_by(cluster) %>%
  mutate(RelAbundance = MeanAbundance  / sum(MeanAbundance ) * 100) %>%
  ungroup()

pdf("plot_velo_bac.pdf", width = 10, height = 5)

ggplot(median_cluster, aes(x = 2, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", colour = "white") +
  geom_text(aes(label = sprintf("%.2f", RelAbundance)), 
            position = position_stack(vjust = 0.5), size = 6, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # creates space in the middle for the donut
  facet_wrap(~cluster) +
  theme_void() +
  labs(title = "Median relative abundance per cluster (donut plot)") +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom")
dev.off()


# D

############################### ANCOMBC AND Volcano plot

library("phyloseq")
library("ANCOMBC")
library("SIAMCAT")
library("ggpubr")
library("reshape2")
library("vegan")
library("dplyr")
library("tidyverse")
library(microbiome)

# Load data
otu = read.delim(file="data/bacteriome.txt", header=T, row.names=1, check.names=F)

meta = read.delim(file="meta.txt", header=T, row.names=1, check.names=F)

# Create placeholder taxonomy table
tax <- data.frame(
  Phylum = rep("_", nrow(otu)),
  Class = rep("_", nrow(otu)),
  Order = rep("_", nrow(otu)),
  Family = rep("_", nrow(otu)),
  Genus = rep("_", nrow(otu)),
  Species = rep("_", nrow(otu))
)
row.names(tax) <- row.names(otu)

# Create phyloseq object
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))
physeq = phyloseq(OTU, TAX)
physeq = merge_phyloseq(physeq, sample_data(meta))

# Normalize by relative abundance
physeq <- transform_sample_counts(physeq, function(x) x / sum(x))

# Run ANCOMBC
out = ancombc(data = physeq,
              formula = "GROUP + IFN_dose + DAY + IFN_dose*DAY+treatment_antibio+treatment_cortico",
              p_adj_method = "fdr", lib_cut = 0,
              group = "GROUP", neg_lb = TRUE, tol = 1e-5,
              max_iter = 1000, conserve = TRUE, alpha = 0.2, global = TRUE)


res = out$res

# Convert results to a data frame
df_res <- data.frame(res)
sig_otus <- subset(df_res, p_val.GROUPC2 < 0.05)
write.table(df_res, file = "ancombc_sig.txt", sep = "\t", col.names = NA)

# Step 1: Extract the significant taxons based on the condition

significant_taxons <- res$p_val$taxon[res$p_val$GROUPC2 <= 0.05]

# Step 2: Create a logical vector that identifies which rows in the OTU table correspond to significant taxons
sig_otus <- rownames(otu) %in% significant_taxons

# Step 3: Subset the OTU table to keep only the significant OTUs
otu_filtered <- otu[sig_otus, , drop = FALSE]

# Create new phyloseq object with filtered OTUs
OTU_filtered = otu_table(as.matrix(otu_filtered), taxa_are_rows = TRUE)
TAX_filtered = tax_table(as.matrix(tax[sig_otus, , drop=FALSE]))
# Create new phyloseq object with filtered OTUs and corresponding taxonomy
physeq_filtered = phyloseq(OTU_filtered, TAX_filtered)
physeq_filtered = merge_phyloseq(physeq_filtered, sample_data(meta))

# Normalize by relative abundance
physeq_filtered <- transform_sample_counts(physeq_filtered, function(x) x / sum(x))

# Create SIAMCAT object
label <- create.label(meta = meta, label = "GROUP", case = "C1", control = "C2")
sc.obj <- siamcat(phyloseq = physeq_filtered, label = label)

# Filter features
sc.obj <- filter.features(sc.obj, cutoff=0.01, filter.method = 'abundance')
sc.obj <- filter.features(sc.obj, cutoff=0.05, filter.method='prevalence', feature.type = 'filtered')

# Check associations
sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.5)

# Generate association plot
body_cols = c("C2" = "limegreen", "C1" = "indianred1")
par(cex = 0.5)

association.plot(sc.obj, 
                 color.scheme = body_cols, sort.by = 'fc', prompt = F,
                 panels = c('fc', 'prevalence', 'auroc'), max.show = 30) 

association.plot(sc.obj, fn.plot = './association_plot_vir_gen_sig.pdf', 
                 color.scheme = body_cols, sort.by = 'fc', max.show = 100,
                 panels = c('fc', 'prevalence', 'auroc'))

associations(sc.obj)

write.table(associations(sc.obj), "./table_list_signature.txt", sep = "\t", col.names = NA)

#### volcano plot bac sig

library(ggplot2)
library(ggrepel)

df <- read.delim("ancombc_sig.txt", header = TRUE)
head(df)
# Add columns for plotting
df$logP <- -log10(df$q_val)
df$direction <- ifelse(df$fc > 0, "P2_high_risk", "P1_low_risk")

df$signif <- df$q_val < 0.05
volcano <- ggplot(df, aes(x = fc, y = -10 * log10(q_val))) +
  
  geom_point(
    alpha = 0.4, size = 5,
    aes(color = ifelse(q_val < 0.05 & fc > 0.5,
                       "P2_high_risk",
                       ifelse(q_val < 0.05 & fc < -0.5,
                              "P1_low_risk",
                              "NS")))
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = ifelse(q_val < 0.05 & (fc > 0.5 | fc < -0.5), Bac, "")),
    size = 5,
    fontface = "bold.italic",
    max.overlaps = 20
  ) +
  
  geom_hline(
    yintercept = -10 * log10(0.05),
    color = "black",
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = c(-0.5, 0.5),
    color = "black",
    linetype = "dotted"
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "",
    y = "-10 * log10(Pvalue)"
  ) +
  
  scale_color_manual(values = c(
    "NS" = "grey70",
    "P1_low_risk" = "green3",
    "P2_high_risk" = "red3"
  )) +
  
  ggtitle("") +
  
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top",
    legend.title = element_blank()
  )
  

volcano

pdf("Volcano_bacsig.pdf",width=8,height=8)
volcano
dev.off()



# E

############### heatmap REORDER BY CLUSTER 

library(pheatmap)

matrix <- read.delim("bacteriome_fucntions_genes.txt", header = TRUE, sep = "\t", row.names = 1, check.names = F)
matrix<-matrix[apply(matrix,1,function(x)length(x[x>0]))>=0.5*ncol(matrix),]

my_sample_col <- read.table("meta.txt", header = TRUE, sep = "\t", row.names = 1)
my_sample_row <- read.table("metagenes.txt", header = TRUE, sep = "\t", row.names = 1)
metadata <- read.delim("meta.txt", check.names = FALSE)
stopifnot(all(colnames(matrix) %in% metadata$SAMPLEID))
metadata_ordered <- metadata[match(colnames(matrix), metadata$SAMPLEID), ]
metadata_sorted <- metadata_ordered[order(metadata_ordered$cluster), ]
matrix_sorted <- matrix[, metadata_sorted$SAMPLEID]
library(RColorBrewer)
# Extract unique functions from your annotation data
function_levels <- unique(my_sample_row$Function)
function_levels <- sort(unique(my_sample_row$Function))
# Assign Set1 palette (recycled if >9 functions)
function_colors <- colorRampPalette(brewer.pal(25, "Set1"))(length(function_levels))
names(function_colors) <- function_levels

ann_colors = list(
  cluster = c("P2_high_risk" = "indianred1", "P1_low_risk" = "limegreen")
)

heatmap <- pheatmap(
  log10(matrix_sorted),
  color=colorRampPalette(c("white", "grey","grey4"))(50),
  annotation_row = my_sample_row,
  annotation_col = my_sample_col,
  cluster_cols = F,
  #cellheight = 50,
  #cellwidth = 20,
  cluster_rows =T,
  annotation_colors = ann_colors,
  fontsize = 12, 
  legend = T,
  annotation_legend = T,
  show_colnames = F, show_rownames = T
)
heatmap

pdf("heatmap.pdf", width = 10, height = 7)
heatmap
dev.off()

# F

######### glm bac funcions cog

########
library(ggplot2)
library(dplyr)
library(broom)
library(purrr)

# Load data
expr <- read.delim("bacteriome_fucntions.txt", header = TRUE, row.names = 1, check.names = FALSE)
group <- read.delim("metadata.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Match samples
common_samples <- intersect(colnames(expr), rownames(group))
expr <- expr[, common_samples]
group <- group[common_samples, , drop = FALSE]

# Binary outcome (C1 vs C2)
y <- ifelse(group$cluster == "C1", 1, 0)
X <- t(expr)

# Covariates as factors
group$IFN_dose <- as.factor(group$IFN_dose)

group$DAY <- as.factor(group$DAY)

group$treatment_antibio <- as.factor(group$treatment_antibio)

group$treatment_cortico <- as.factor(group$treatment_cortico)

# ---- GLM per Function ----
glm_results <- lapply(colnames(X), function(g) {
  
  df_model <- data.frame(
    y = y,
    gene = X[, g],
    
    IFN_dose = group$IFN_dose,
    
    DAY = group$DAY,
    
    treatment_antibio = group$treatment_antibio,
    
    treatment_cortico = group$treatment_cortico
  )
  
  fit <- glm(
    y ~ gene +
      IFN_dose +
      treatment_antibio + DAY +
      treatment_cortico,
    
    data = df_model,
    family = binomial
  )
  
  coef_summary <- summary(fit)$coefficients
  
  data.frame(
    Function = g,
    
    estimate = coef_summary["gene", "Estimate"],
    
    pval = coef_summary["gene", "Pr(>|z|)"]
  )
})

# Combine results
glm_results <- do.call(rbind, glm_results)

# Adjust p-values
glm_results$padj <- p.adjust(glm_results$pval, method = "fdr")

# Sort by significance
glm_results <- glm_results[order(glm_results$padj), ]

head(glm_results)

glm_df <- glm_results

# Use raw p-values for significance
glm_df <- glm_df %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    ),
    direction = ifelse(estimate > 0, "C1", "C2")
  )

# Keep only significant
glm_sig <- glm_df %>%
  filter(sig != "")

# Order by effect size
glm_sig <- glm_sig %>%
  arrange(estimate)

# ---- LEfSe-like plot ----

pdf("glm_bac_func.pdf",width=13,height=5)

ggplot(glm_sig, aes(x = reorder(Function, estimate), y = estimate, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 20) +
  labs(
    title = "Significant Features (GLM, raw p<0.05)",
    x = NULL, 
    y = "Effect size (log-odds)"
  ) +
  geom_text(aes(label = sig),
            hjust = ifelse(glm_sig$estimate > 0, -0.2, 1.2),
            size = 10) +
  scale_fill_manual(values = c("C1" = "indianred1", "C2" = "limegreen")) +
  theme(
    axis.text.y = element_text(size = 11, face = "bold")
  )

dev.off()

library(vegan)
library(ggplot2)
library(dplyr)

# Distance matrix (Bray-Curtis, common in expression/omics)
dist_mat <- vegdist(t(expr), method = "bray")

# Run PERMANOVA (multivariate GLM-like test)
permanova_res <- adonis2(dist_mat ~ cluster, data = group)

print(permanova_res)


#G 

############FFT tree

library(FFTrees)
library(caret)
library(pROC)
library(ROCR)

# Load your data
Bac <- read.delim("Bacteriome.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
VAL <- read.delim("Bacteriome.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Ensure outcome variable exists in Bac and VAL (VAL can be NA for prediction)
# Assuming your outcome column is named 'C1'
VAL$C1 <- NA  

# Set seed for reproducibility
set.seed(123)

# Train FFTree on the training set
tree <- FFTrees(formula = C1 ~ ., 
                data = Bac,
                decision.labels = c("pneumotype1", "pneumotype2"),
                main = "HAP pneumotype prediction")

# Plot the trained FFTree
plot(tree)

# Predict on external validation dataset
# type = "class" gives the predicted class
# type = "path" gives the decision path (which variables were used)
hap_predicted <- predict(
  object = tree,
  newdata = VAL,
  tree = 1,
  type = "class"
)

VAL$cluster <- hap_predicted

table(hap_predicted)

metaibis <- read.delim("metadata.txt", header = TRUE, stringsAsFactors = FALSE)

# Join lengths to ids based on 'id'
joined_data <- left_join(VAL, metaibis, by = "SAMPLEID")

# View the result
print(joined_data)

write.table(joined_data, "./VAL_with_predictions.txt", sep = "\t", row.names = F)

#H 

############################### 
###
library(phyloseq)
library(ANCOMBC)
library(SIAMCAT)
library(caret)
library(dplyr)
library(pROC)
library(ggplot2)

set.seed(123)

# Load data
otu = read.delim("bacteriome.txt", header=T, row.names=1, check.names=F)
meta = read.delim("metadata.txt", header=T, row.names=1, check.names=F)

# Ensure matching sample IDs
common_samples <- intersect(colnames(otu), rownames(meta))

otu <- otu[, common_samples]
meta <- meta[common_samples, ]

# Optional but recommended: enforce same order
meta <- meta[colnames(otu), ]

# ============================================
# Helper function for safe ROC smoothing
# ============================================
safe_smooth_roc <- function(response, predictor, min_points = 4) {
  # Create ROC object
  roc_obj <- roc(response, predictor, quiet = TRUE)
  
  # Check if we have enough points for smoothing
  n_points <- length(roc_obj$sensitivities)
  
  if (n_points >= min_points) {
    # Try binormal smoothing first
    tryCatch({
      roc_smooth <- roc(response, predictor, smooth = TRUE, 
                        smooth.method = "binormal", quiet = TRUE)
      return(roc_smooth)
    }, error = function(e) {
      # Fall back to density smoothing
      tryCatch({
        roc_smooth <- roc(response, predictor, smooth = TRUE, 
                          smooth.method = "density", quiet = TRUE)
        return(roc_smooth)
      }, error = function(e2) {
        # Return unsmoothed ROC
        message("  Using unsmoothed ROC curve for this fold")
        return(roc_obj)
      })
    })
  } else {
    # Not enough points, use unsmoothed
    message(paste("  Not enough points (", n_points, ") for smoothing, using unsmoothed ROC"))
    return(roc_obj)
  }
}

# ============================================
# PART 1: ROC with 10-fold cross-validation
# ============================================

# Outer folds
outer_folds <- createFolds(meta$GROUP, k = 10, list = TRUE)

auc_outer <- c()
feature_list <- list()
all_roc_data <- list()
cv_pred_store <- list()

for (i in 1:10) {
  cat("\nOuter fold:", i, "\n")
  
  test_idx <- outer_folds[[i]]
  train_idx <- setdiff(1:nrow(meta), test_idx)
  
  meta_train <- meta[train_idx, ]
  meta_test  <- meta[test_idx, ]
  
  otu_train <- otu[, rownames(meta_train)]
  otu_test  <- otu[, rownames(meta_test)]
  
  # ---- PHYLOSEQ TRAIN ----
  OTU = otu_table(as.matrix(otu_train), taxa_are_rows = TRUE)
  TAX = tax_table(matrix("_", nrow=nrow(otu_train), ncol=6,
                         dimnames=list(rownames(otu_train),
                                       c("Phylum","Class","Order","Family","Genus","Species"))))
  physeq_train = phyloseq(OTU, TAX, sample_data(meta_train))
  physeq_train <- transform_sample_counts(physeq_train, function(x) x / sum(x))
  
  # ---- ANCOMBC ONLY ON TRAIN ----
  # Suppress messages from ancombc
  suppressMessages({
    out = tryCatch({
      ancombc(data = physeq_train,
              formula = "GROUP", p_adj_method = "fdr", 
              lib_cut = 0, group = "GROUP", neg_lb = TRUE, 
              tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.2, global = TRUE)
    }, error = function(e) {
      cat("  ANCOMBC error:", e$message, "\n")
      return(NULL)
    })
  })
  
  if (is.null(out)) next
  
  res = out$res
  
  # Check if p_val exists and has GROUPC2 column
  if (is.null(res$p_val) || !"GROUPC2" %in% colnames(res$p_val)) {
    cat("  No significant taxa found (p_val table issue)\n")
    next
  }
  
  sig_taxa <- res$p_val$taxon[res$p_val$GROUPC2 <= 0.05]
  
  if (length(sig_taxa) < 2) {
    cat("  Only", length(sig_taxa), "significant taxa found, skipping fold\n")
    next
  }
  
  cat("  Found", length(sig_taxa), "significant taxa\n")
  feature_list[[i]] <- sig_taxa
  
  # ---- FILTER DATA ----
  otu_train_f <- otu_train[rownames(otu_train) %in% sig_taxa, ]
  otu_test_f  <- otu_test[rownames(otu_test) %in% sig_taxa, ]
  
  # Check if we have any features after filtering
  if (nrow(otu_train_f) == 0 || nrow(otu_test_f) == 0) {
    cat("  No features remaining after filtering, skipping fold\n")
    next
  }
  
  # transpose for caret
  train_df <- t(otu_train_f)
  test_df  <- t(otu_test_f)
  
  train_df <- as.data.frame(train_df)
  test_df  <- as.data.frame(test_df)
  
  train_df$GROUP <- meta_train$GROUP
  test_df$GROUP  <- meta_test$GROUP
  
  # Check for near-zero variance predictors
  nzv <- nearZeroVar(train_df[, -ncol(train_df)], saveMetrics = TRUE)
  if (sum(!nzv$nzv) < 2) {
    cat("  Too few non-zero variance predictors, skipping fold\n")
    next
  }
  
  # ---- INNER CV (model tuning) ----
  ctrl <- trainControl(
    method = "cv",
    number = min(5, nrow(train_df) - 1),  # Adjust CV folds based on sample size
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  model <- tryCatch({
    train(
      GROUP ~ .,
      data = train_df,
      method = "glmnet",
      metric = "ROC",
      trControl = ctrl
    )
  }, error = function(e) {
    cat("  Model training error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(model)) next
  
  # ---- TEST PREDICTION ----
  pred <- predict(model, test_df, type = "prob")
  
  cv_pred_store[[i]] <- data.frame(
    label = test_df$GROUP,
    score = pred$C1,
    fold = i
  )
  
  # Safe ROC calculation with smoothing
  roc_obj <- safe_smooth_roc(test_df$GROUP, pred$C1)
  current_auc <- auc(roc_obj)
  auc_outer[i] <- current_auc
  
  cat("  Fold AUC:", round(current_auc, 3), "\n")
  
  # Store ROC data for plotting
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    fold = paste0("Fold ", i),
    auc = round(current_auc, 3)
  )
  all_roc_data[[i]] <- roc_df
}

cv_all_preds <- do.call(rbind, cv_pred_store)

# Remove NULL entries from all_roc_data
all_roc_data <- all_roc_data[!sapply(all_roc_data, is.null)]

# Check if we have any successful folds
if (length(all_roc_data) == 0) {
  stop("No successful folds completed. Check your data and ANCOMBC settings.")
}

# Combine all ROC data
roc_all_folds <- do.call(rbind, all_roc_data)

# Calculate mean and SD of AUC
auc_outer_valid <- auc_outer[!is.na(auc_outer) & auc_outer > 0]
mean_auc <- mean(auc_outer_valid, na.rm = TRUE)
sd_auc <- sd(auc_outer_valid, na.rm = TRUE)

cat("\n========================================\n")
cat("Nested CV AUC:", round(mean_auc, 3), "±", round(sd_auc, 3), "\n")
cat("Successful folds:", length(auc_outer_valid), "/", 10, "\n")
cat("========================================\n")

# Calculate mean ROC curve for CV folds
if (length(all_roc_data) > 0) {
  common_fpr <- seq(0, 1, length.out = 200)
  all_tpr_interp <- matrix(NA, nrow = length(common_fpr), ncol = length(all_roc_data))
  
  for (i in 1:length(all_roc_data)) {
    if (!is.null(all_roc_data[[i]]) && nrow(all_roc_data[[i]]) > 0) {
      # Interpolate TPR at common FPR points
      all_tpr_interp[, i] <- approx(all_roc_data[[i]]$fpr, 
                                    all_roc_data[[i]]$tpr, 
                                    xout = common_fpr, 
                                    rule = 2, 
                                    ties = mean)$y
    }
  }
  
  # Calculate mean and 95% confidence interval
  mean_tpr <- rowMeans(all_tpr_interp, na.rm = TRUE)
  ci_lower <- apply(all_tpr_interp, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
  ci_upper <- apply(all_tpr_interp, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
  
  # Create mean ROC data frame
  mean_roc_df <- data.frame(
    fpr = common_fpr,
    tpr = mean_tpr,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  
  # Apply loess smoothing for cleaner lines
  mean_roc_df$tpr_smooth <- predict(loess(tpr ~ fpr, data = mean_roc_df, span = 0.3))
  
  # SIAMCAT-style ROC plot with 10-fold CV
  plot_cv <- ggplot() +
    geom_line(data = roc_all_folds, 
              aes(x = fpr, y = tpr, group = fold, color = "Individual Folds"), 
              size = 0.5, alpha = 0.5) +
    geom_ribbon(data = mean_roc_df, 
                aes(x = fpr, ymin = ci_lower, ymax = ci_upper, fill = "95% CI"), 
                alpha = 0.2) +
    geom_line(data = mean_roc_df, 
              aes(x = fpr, y = tpr_smooth, color = "Mean ROC"), 
              size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray50", size = 0.8, alpha = 0.7) +
    scale_color_manual(name = "ROC Curves",
                       values = c("Individual Folds" = "gray70",
                                  "Mean ROC" = "red")) +
    scale_fill_manual(name = "", values = c("95% CI" = "red3")) +
    labs(title = "ROC Curve - 10-Fold Cross-Validation",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)",
         caption = paste0("Mean AUC = ", round(mean_auc, 3), " ± ", round(sd_auc, 3), 
                          " (", length(auc_outer_valid), " successful folds)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = c(0.7, 0.2),
          legend.background = element_rect(fill = "white", color = "gray80"),
          legend.key.size = unit(0.8, "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90"),
          plot.caption = element_text(hjust = 0, size = 9))
  
  print(plot_cv)
} else {
  cat("No ROC data available for plotting\n")
}


plot_cv <- ggplot() +
  geom_line(data = roc_all_folds, 
            aes(x = fpr, y = tpr, group = fold, color = "Individual Folds"), 
            size = 1, alpha = 0.5) +
  geom_ribbon(data = mean_roc_df, 
              aes(x = fpr, ymin = ci_lower, ymax = ci_upper, fill = "95% CI"), 
              alpha = 0.2) +
  geom_line(data = mean_roc_df, 
            aes(x = fpr, y = tpr_smooth, color = "Mean ROC"), 
            size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "gray50", size = 0.8, alpha = 0.7) +
  scale_color_manual(name = "ROC Curves",
                     values = c("Individual Folds" = "gray70",
                                "Mean ROC" = "red")) +
  scale_fill_manual(name = "", values = c("95% CI" = "red")) +
  labs(title = "ROC Curve - 10-Fold Cross-Validation",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       caption = paste0("Mean AUC = ", round(mean_auc, 3), " ± ", round(sd_auc, 3), 
                        " (", length(auc_outer_valid), " successful folds)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = c(0.7, 0.2),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.key.size = unit(0.8, "cm"),
        plot.caption = element_text(hjust = 0, size = 9))

print(plot_cv)

pdf(file = "smooth_roc_bac.pdf")
plot_cv
dev.off()

#I  

library(survival)
library(survminer)

datat <- VAL_with_predictions

# Create a survival object
surv_object <- Surv(time = datat$icu, event = datat$Mortality)

# Fit the Kaplan-Meier model
km_fit <- survfit(surv_object ~ cluster, data = datat)

summary(km_fit)

# Plot Kaplan-Meier curve
km_plot <- ggsurvplot(
  km_fit,                       
  data = datat,  
  risk.table = TRUE,         
  pval = TRUE,              
  conf.int = F,           
  break.time.by = 10,     
  ggtheme = theme_bw(), 
  risk.table.y.text.col = TRUE, 
  risk.table.y.text = TRUE,
  ylab = "Probability of Survival",
  xlab ="Time (days)",
  alette = cluster_colors,
  linetype = "cluster"
)

cox_model <- coxph(Surv(time = icu, event = Mortality) ~ cluster, data = datat)

# Extract HR and 95% CI
cox_summary <- summary(cox_model)
hr <- round(cox_summary$coefficients[1, "exp(coef)"], 1)    
ci_lower <- round(cox_summary$conf.int[1, "lower .95"], 2)     
ci_upper <- round(cox_summary$conf.int[1, "upper .95"], 2)   

km_plot$plot <- km_plot$plot +
  annotate(
    "text", x = 2, y = 0.1, 
    label = paste0("HR: ", hr, " (95% CI: ", ci_lower, "-", ci_upper, ")"),
    size = 5, hjust = 0
  ) + coord_cartesian(xlim = c(0, 90)) 

km_plot

























