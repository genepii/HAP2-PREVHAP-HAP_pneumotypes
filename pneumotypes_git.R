
### Donut plot 

library(dplyr)
library(ggplot2)
library(reshape2)

# Load your abundance matrix
matrix <- read.delim("Relab_ibis_signew.txt", header = TRUE, sep = "\t", check.names = F)

# Load metadata (Sample -> Cluster)
metadata <- read.delim("filt_bacterial_HAP.txt", header = TRUE, sep = "\t", check.names = F)

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

# Plot bar plot RPKM

ggplot(median_cluster, aes(x = cluster, y = MeanAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "white") +
  labs(title = "Median relative abundance per cluster",
       x = "Cluster",
       y = "Relative Abundance") +
  theme_classic() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14)
  )


# Plot stacked 100% bar plot
ggplot(median_cluster, aes(x = cluster, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stak", colour = "white") +
  #geom_text(aes(label = MeanAbundance ), position = position_stack(vjust = 0.5), size = 6, color = "black") +
  #scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(title = "Median relative abundance per cluster (100%)",
       x = "Cluster",
       y = "Relative Abundance (%)") +
  theme_classic() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 12),
        legend.position = "bottom",
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14))

# Circular (stacked) plot per cluster
ggplot(median_cluster, aes(x = "", y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", colour = "white") +
  #scale_fill_viridis(discrete = TRUE, option = "D") +
  coord_polar(theta = "y") +
  facet_wrap(~cluster) +
  theme_void() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 12),
        legend.position = "bottom")

pdf("plot_velo_bac_ibis.pdf", width = 10, height = 5)

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


############################### ANCOMBC AND SIAMCAT

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
otu = read.delim(file="dRPKM_genus_ibis2.txt", header=T, row.names=1, check.names=F)

#otu<-otu[apply(otu,1,function(x)length(x[x>0]))>=0.5*ncol(otu),]

meta = read.delim(file="metaibis.txt", header=T, row.names=1, check.names=F)

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
              formula = "GROUP",
              p_adj_method = "fdr", lib_cut = 0,
              group = "GROUP", neg_lb = TRUE, tol = 1e-5,
              max_iter = 100, conserve = TRUE, alpha = 0.2, global = TRUE)

res = out$res

# Convert results to a data frame
#df_res <- data.frame(res)
#sig_otus <- subset(df_res, p_val.GROUPC2 < 0.05 | p_val.GROUPC3 < 0.05)
#write.table(sig_otus, file = "significant_otus_bac_d7.txt", sep = "\t", col.names = NA)


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

#otu_data_clean <- otu_data[, colSums(is.na(otu_data)) == 0]
#OTU_clean <- otu_table(as.matrix(otu_data_clean), taxa_are_rows = TRUE)
#physeq_clean <- phyloseq(OTU_clean, tax_table(physeq_filtered), sample_data(physeq_filtered))


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
                 panels = c('fc', 'prevalence', 'auroc'), max.show = 100) 

association.plot(sc.obj, fn.plot = './association_plot_vir_gen_sig.pdf', 
                 color.scheme = body_cols, sort.by = 'fc', max.show = 100,
                 panels = c('fc', 'prevalence', 'auroc'))

associations(sc.obj)

write.table(associations(sc.obj), "./table_list_signature_BAC.txt", sep = "\t", col.names = NA)

# Test the model – model building – plot the ROC curve – interpretation

sc.obj <- normalize.features(
  sc.obj,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p =2, norm.margin=1))

sc.obj <- create.data.split(
  sc.obj,
  num.folds = 10,
  num.resample = 10)

sc.obj <- train.model(
  sc.obj,
  method = "lasso")

sc.obj <- train.model(sc.obj, method = "randomForest",param.set = list(ntree = 1000, mtry = 5))

model_type(sc.obj)
models <- models(sc.obj)
models[[1]]
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)
sc.obj <- evaluate.predictions(sc.obj)

model.evaluation.plot(sc.obj, show.all = T, colours = "black")

model.evaluation.plot(sc.obj, show.all = T, colours = "black",fn.plot = './roc_plot_VIR.pdf')


######################
#valid

# Load data
otu = read.delim(file="Bac_species_ibis.txt", header=T, row.names=1, check.names=F)

meta = read.delim(file="metaibis.txt", header=T, row.names=1, check.names=F)

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

# Normalize by relative abundance
physeq_filtered <- transform_sample_counts(physeq, function(x) x / sum(x))

# Create SIAMCAT object
label <- create.label(meta = meta, label = "GROUP", case = "C1", control = "C2")
vc.obj <- siamcat(phyloseq = physeq_filtered, label = label)

# Filter features
vc.obj <- filter.features(vc.obj, cutoff=0.01, filter.method = 'abundance')
vc.obj <- filter.features(vc.obj, cutoff=0.05, filter.method='prevalence', feature.type = 'filtered')

# Check associations
vc.obj <- check.associations(vc.obj, log.n0 = 1e-06, alpha = 0.5)

# Generate association plot
body_cols = c("C2" = "BLUE", "C1" = "red2")
association.plot(vc.obj, 
                 color.scheme = body_cols, sort.by = 'fc', prompt = F, 
                 panels = c('fc', 'prevalence', 'auroc'), max.show = 50)
associations(vc.obj)

# Test the model – model building – plot the ROC curve – interpretation
# Normalize training object

sc.obj <- normalize.features(sc.obj, norm.method = "log.unit", norm.param = list(log.n0 = 1e-06, n.p = 2, norm.margin = 1))

# Normalize validation object using parameters from training normalization to keep data comparable
vc.obj <- normalize.features(vc.obj, norm.method = "log.unit", norm.param = list(log.n0 = 1e-06, n.p = 2, norm.margin = 1))


# Create data split in training object for cross-validation during training
sc.obj <- create.data.split(sc.obj, num.folds = 1, num.resample = 1)

# Train model on training object
sc.obj <- train.model(sc.obj, method = "lasso")

sc.obj <- train.model(sc.obj, method = "randomForest",param.set = list(ntree = 1000, mtry = 5))

# Predict on validation object using trained model from training object
vc.obj <- make.predictions(sc.obj, siamcat.holdout = vc.obj, normalize.holdout = TRUE)

# Evaluate predictions on validation object
vc.obj <- evaluate.predictions(vc.obj)

# Plot evaluation results on validation object
model.evaluation.plot(vc.obj, show.all = T, colours = "BLACK")

model.evaluation.plot(vc.obj, show.all = F, colours = "black",fn.plot = './roc_plot_VIR_valid2.pdf')

###################### trials
# Step 1: Repeated CV on training data
sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 10)
sc.obj <- train.model(sc.obj, method = "lasso")
sc.obj <- make.predictions(sc.obj)
sc.obj <- evaluate.predictions(sc.obj)

# Step 2: Train final model on full training data
sc.obj <- create.data.split(sc.obj, num.folds = 1, num.resample = 1)
sc.obj <- train.model(sc.obj, method = "lasso")

# Step 3: External validation dataset prep
# (subset validation features to match training etc., normalize vc.obj with sc.obj parameters)
vc.obj <- make.predictions(sc.obj, siamcat.holdout = vc.obj, normalize.holdout = TRUE)
vc.obj <- evaluate.predictions(vc.obj)
model.evaluation.plot(vc.obj)
###############################


############FFT tree

library(FFTrees)
library(caret)
library(pROC)
library(ROCR)

# Load your data
Bac <- read.delim("Relab_bac_species_prev_sig_FFT.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
VAL <- read.delim("Relab_bac_species_ibis_sig_FFT.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

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

metaibis <- read.delim("metaibis.txt", header = TRUE, stringsAsFactors = FALSE)

# Join lengths to ids based on 'id'
joined_data <- left_join(VAL, metaibis, by = "SAMPLEID")

# View the result
print(joined_data)

#write.table(joined_data, "./VAL_with_predictions.txt", sep = "\t", row.names = F)

###########################

###pcoa

library(vegan)
set.seed(1)
data1 <- read.delim("dRPKM.txt", row.names = 1)
dataTransposed1 <- t(data1)

metadata <- read.delim("meta.txt", row.names = 1)

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
adonis.test1.cont.subjID <- adonis2(dist.1 ~ cluster, data = metadata, permutations = 9999)
adonis.test1.cont.subjID

pval <- adonis.test1.cont.subjID$`Pr(>F)`[1]
R2   <- adonis.test1.cont.subjID$R2[1]

pval_text <- paste0("p = ", signif(pval, 3))
R2_text   <- paste0("R² = ", signif(R2, 3))

# Betadisper
dispersion <- betadisper(dist.1, group = metadata$cluster)
scores_d <- as.data.frame(scores(dispersion, display = "sites"))

scores_d$cluster <- metadata$cluster

colors <- c("C1" = "indianred1", "C2" = "limegreen")

# % variance explained
eig <- dispersion$eig / sum(dispersion$eig)
xlab1 <- paste0("PCoA1 (", round(eig[1] * 100, 1), "%)")
ylab1 <- paste0("PCoA2 (", round(eig[2] * 100, 1), "%)")

# PLOT
#pdf("pcoa.pdf", width = 7, height = 7)

p <- ggplot(scores_d, aes(x = PCoA1, y = PCoA2, color = cluster)) +
  geom_point(size = 4) +
  stat_ellipse(type = "norm", linewidth = 1.3) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 18) +
  labs(
    title = "Bray-Curtis PCoA",
    x = xlab1,
    y = ylab1
  ) +
  annotate("text", x = min(scores_d$PCoA1), y = max(scores_d$PCoA2),
           label = pval_text, hjust = 0.3, vjust = 1, size = 6) +
  annotate("text", x = min(scores_d$PCoA1), y = max(scores_d$PCoA2)*0.92,
           label = R2_text, hjust = 0.3, vjust = 2, size = 6)

print(p)

#dev.off()


#####################

# Alphadiv metrics :

library(vegan)
library(ggplot2)

## refomrlat meta as otu matrix
data<-read.delim("dRPKM.txt", row.names = 1,check.names = F)

#data<-data[apply(data,1,function(x)length(x[x>0]))>=0.5*ncol(data),]

data<-t(data)
meta<-read.delim("meta.txt", row.names = 1, check.names = F)
meta <- meta[rownames(meta) %in% rownames(data), , drop = FALSE]
meta <- meta[match(rownames(data), rownames(meta)), , drop = FALSE]

data_evenness <- diversity(data) / log(specnumber(data))            

data_shannon <- diversity(data, index = "shannon")                          

data=round(data, digits = 0)

data_richness <- estimateR(data)

data_alphadiv <- cbind(meta, t(data_richness), data_shannon, data_evenness) 

rm(data_richness, data_evenness, data_shannon)                                

head(data_alphadiv)
data <-data_alphadiv 

#write.table(data_alphadiv,"data_alphadiv_VIR.txt", sep = '\t', col.names = NA)

### wilcox test + plot shannon+richness at once 

### Mann–Whitney U test for two variables (data_shannon & S.obs)

library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)  

# -----------------------------
# 🔧 Variables
group_var <- "cluster"
num_vars  <- c("data_shannon", "S.obs") 
# -----------------------------

# Custom colors (C1 and C2)
custom_colors <- c("C2" = "limegreen", "C1" = "indianred1")

# Keep only C1 and C2 samples
data <- data %>% filter(.data[[group_var]] %in% c("C1", "C2"))

# Reorder factor levels
data[[group_var]] <- factor(data[[group_var]], levels = c("C1", "C2"))

# Function to run MWU + plot for one variable
make_plot <- function(df, var, group_var, custom_colors) {
  
  # Mann–Whitney U test
  res.mannwhitney <- wilcox.test(
    as.formula(paste(var, "~", group_var)),
    data = df,
    exact = FALSE
  )
  
  # Tidy result for annotation
  res.wilcox <- df %>%
    wilcox_test(as.formula(paste(var, "~", group_var))) %>%
    add_significance() %>%
    add_xy_position(x = group_var)
  
  # Boxplot
  p <- ggboxplot(
    df, x = group_var, y = var, fill = group_var,
    alpha = 0.5, width = 0.2
  ) +
    geom_jitter(aes_string(color = group_var), width = 0.03, size = 3, alpha = 0.4) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    labs(
      title = paste("Mann–Whitney U Test"),
      subtitle = paste0("p-value: ", signif(res.mannwhitney$p.value, digits = 3)),
      x = "", y = var
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 14, colour = "black"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  
  # Add significance stars if any p < 0.05
  if (any(res.wilcox$p < 0.05)) {
    p <- p + stat_pvalue_manual(
      res.wilcox, label = "p.signif", hide.ns = TRUE, size = 5,
      tip.length = 0.01
    )
  }
  
  return(p)
}

# Generate plots for both variables
p1 <- make_plot(data, "data_shannon", group_var, custom_colors)
p2 <- make_plot(data, "S.obs", group_var, custom_colors)

# Combine side by side
final_plot <- p1 + p2
final_plot

pdf("Alphadiv_vir_clusters.pdf",width=10,height=5)
print(final_plot)
dev.off()



#########################

############## nbglm volcano gene symbol
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(cowplot)

# Read expression data and metadata
expression_df <- read.delim("Raw_coding_PREVHAP_nodup.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
metadata <- read.delim("filtered_data.txt", stringsAsFactors = FALSE)

sapply(metadata$cluster,function(x)ifelse(x=="N1",1,0))->metadata$C1
sapply(metadata$cluster,function(x)ifelse(x=="N2",1,0))->metadata$C2

# Aggregate expression data by gene
aggregate(. ~ Gene, expression_df, FUN = sum) -> expression_df

# Set gene ID as rownames
expression_df <- expression_df %>%
  tibble::column_to_rownames("Gene")

# Reorder expression data to match metadata
expression_df <- as.data.frame(expression_df) %>%
  dplyr::select(metadata$SAMPLEID)

# Check if the column order matches metadata
all.equal(colnames(expression_df), metadata$SAMPLEID)

# Filter genes based on non-zero expression
filtered_expression_df <- expression_df[apply(expression_df, 1, function(x) length(x[x > 0])) > 0* ncol(expression_df), ]

# Round values to integers (for edgeR)
filtered_expression_df <- round(filtered_expression_df)

# Normalize counts (edgeR)
y <- DGEList(filtered_expression_df)
y <- calcNormFactors(y)

# Export normalized matrix (CPM values)
norm_counts <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
#write.table(norm_counts, file = "normalized_counts_matrix_nohap.txt", sep = "\t", quote = FALSE, col.names = NA)

### NEW STEP: Subset to inflammation genes
infl_genes <- read.delim("infla_list.txt", header = FALSE, stringsAsFactors = FALSE)
infl_genes <- infl_genes[[1]]   # assuming first column has the gene symbols

# Keep only inflammation genes present in the dataset
keep_genes <- rownames(y) %in% infl_genes
y <- y[keep_genes, , keep.lib.sizes = FALSE]

y <- y[rowSums(y$counts) > 0, , keep.lib.sizes = FALSE]
y <- y[, colSums(y$counts) > 0, keep.lib.sizes = FALSE]

metadata <- metadata[metadata$SAMPLEID %in% colnames(y), ]
metadata <- metadata[match(colnames(y), metadata$SAMPLEID), ]

# Design matrix for clusters (assuming 'cluster' is a column in 'metadata')
design <- model.matrix(~ C1, data = metadata)
y <- estimateDisp(y, design, robust = TRUE)

# Fit model and perform quasi-likelihood F-test
fit <- glmQLFit(y, design, robust = TRUE)
qlt <- glmQLFTest(fit)

# Get top genes (restricted to inflammation list)
topgenes <- topTags(qlt, n = nrow(y))

# Identify enriched and depleted genes
tol10b.enriched <- topgenes$table$logFC > 0 & topgenes$table$PValue < 0.05
tol10b.depleted <- topgenes$table$logFC < 0 & topgenes$table$PValue < 0.05

# Extract results as a plain data.frame
tg_df <- topgenes$table

# Add Ensembl IDs as a column (they are stored separately)
tg_df$Ensembl <- rownames(tg_df)

# Read mapping file
ens_map <- read.delim("ens_to_symbol.txt", header = TRUE, stringsAsFactors = FALSE)

# Merge results with symbol mapping
tg_df <- merge(tg_df, ens_map, by = "Ensembl", all.x = TRUE)

# Create label column (prefer symbol, fallback = Ensembl)
tg_df$Label <- ifelse(is.na(tg_df$Symbol), tg_df$Ensembl, tg_df$Symbol)

# (Optional) reset rownames for convenience
rownames(tg_df) <- tg_df$Ensembl

write.table(topgenes, file = "inputList_prevhap_nodup_before_hap.txt", sep = "\t", col.names = NA)

##################################################################

volcano <- ggplot(tg_df, aes(x = logFC, y = -10 * log10(PValue))) +
  geom_point(alpha = 0.6, size = 4,
             aes(color = ifelse(PValue < 0.05,
                                ifelse(logFC > 0, "red", "green"),
                                "black"))) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(PValue < 0.05 & abs(logFC) > 1, Label, "")),
    size = 5,
    max.overlaps = 20
  ) +
  geom_hline(yintercept = -10 * log10(0.05), color = "black",linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "black",linetype = "dashed") +
  theme_bw() +
  labs(x = "Log2 Fold Change", y = "-10 * log10(P-Value)") +
  scale_color_manual(values = c("grey", "LIGHTBLUE", "BLUE")) +
  ggtitle("Volcano Plot") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top",
    legend.title = element_blank())
volcano

pdf("Volcano_clusters.pdf",width=20,height=15)
volcano
dev.off()

########## enrichissemnt go
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)


genes_df <- read.delim("inputList_ibis.txt", header = TRUE, check.names = FALSE)

head(genes_df)

# Split into C1 and C2 gene lists
genes_C1 <- genes_df$ID[genes_df$Group == "C1"]
genes_C2 <- genes_df$ID[genes_df$Group == "C2"]

# Run enrichment separately
ego_C1 <- enrichGO(gene         = genes_C1,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = 'ENSEMBL',
                   ont          = "BP",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

ego_C2 <- enrichGO(gene         = genes_C2,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = 'ENSEMBL',
                   ont          = "BP",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

# Add group label for plotting
df_C1 <- as.data.frame(ego_C1) %>% mutate(Group = "predicted high risk")
df_C2 <- as.data.frame(ego_C2) %>% mutate(Group = "predicted low risk")

# Combine results
ego_df <- rbind(df_C1, df_C2)

# Keep top 10 terms per group for clarity
ego_df <- ego_df %>%
  group_by(Group) %>%
  top_n(-50, p.adjust) %>%
  ungroup()

# Barplot with ggplot2
ggplot(ego_df, aes(x = reorder(Description, -Count), 
                   y = Count, 
                   fill = p.adjust)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "GO Terms", y = "Gene Count") +
  theme_bw(base_size = 20) +
  theme(axis.text.y = element_text(size = 18)) + facet_wrap(~Group)



# Add group label
ego_C1@result$Group <- "predicted high risk"
ego_C2@result$Group <- "C2"

# Combine results into one object
ego_combined <- ego_C1
ego_combined@result <- rbind(ego_C1@result, ego_C2@result)

pdf("GO_clusters.pdf",width=20,height=10)
# Dotplot comparison
dotplot(ego_combined, split = "Group") + 
  facet_wrap(~Group) +
  theme_linedraw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),    
    strip.text = element_text(size = 20, face = "bold") 
  )+
  scale_size(range = c(3, 20))
dev.off()



#########PLSDA
# SPLS

library(mixOmics)


Clinical <- read.delim("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
Y <- factor(Clinical$cluster)
Y_numeric <- as.numeric(factor(Y))
summary(Y)

X1 <- read.delim("data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = F)
X1<-t(X1)


result.plsda.srbct <- plsda(X1, Y)


plotIndiv(result.plsda.srbct, legend = T, ind.names = FALSE, star = T, 
          ellipse = T, centroid = T , comp = 1:2, col = c("limegreen", "indianred1"))

plsda.res <- plsda(X1, Y, ncomp = 10)
#perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = 10)
#plot(perf.plsda, sd = TRUE, legend.position = 'horizontal')

background.max <- background.predict(plsda.res, 
                                     comp.predicted = 2,
                                     dist = 'max.dist') 
pdf("plsda_clusters_ibis.pdf",width=10,height=10)
plotIndiv(plsda.res, comp = 1:2,
          group = Clinical$cluster,
          star = TRUE, style = "lattice",
          col = c("indianred1", "limegreen"),
          ind.names = FALSE, title = 'Maximum distance',
          legend = T, background = background.max)
dev.off()


######### glm bac funcions cog

########
library(ggplot2)
library(dplyr)
library(broom)
library(purrr)

# Load data
expr <- read.delim("data.txt", header = TRUE, row.names = 1, check.names = FALSE)
group <- read.delim("meta.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Match samples
common_samples <- intersect(colnames(expr), rownames(group))
expr <- expr[, common_samples]
group <- group[common_samples, , drop = FALSE]

# Binary outcome (C1 vs C2)
y <- ifelse(group$cluster == "C1", 1, 0)
X <- t(expr)

# ---- GLM per gene ----
glm_results <- lapply(colnames(X), function(g) {
  fit <- glm(y ~ X[, g], family = binomial)
  coef_summary <- summary(fit)$coefficients
  data.frame(
    gene = g,
    estimate = coef_summary[2, 1],
    pval = coef_summary[2, 4]
  )
})

glm_df <- do.call(rbind, glm_results)

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
glm_sig <- glm_df %>% filter(sig != "")

# Order by effect size
glm_sig <- glm_sig %>% arrange(estimate)

# ---- LEfSe-like plot ----

pdf("glm_avir_func.pdf",width=13,height=5)

ggplot(glm_sig, aes(x = reorder(gene, estimate), y = estimate, fill = direction)) +
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


############## clinical univariate analysis


library(survival)
library(dplyr)
library(purrr)

data <- read.delim("prevhap_db_selected_article_forest.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Define variables to test separately
vars_to_test <- c(
  "cluster",
  "com_age",
  "com_admission_type",
  "com_immunosuppression",
  "com_igs2_admission",
  "com_temp",
  "com_dysleuco",
  "com_pf_ratio",
  "com_haemophilus",
  "com_mssa",
  "com_pneumoco",
  "Corticosteroid",
  "Interferon"
)

# Run univariate Cox models
cox_results <- map_dfr(vars_to_test, function(var) {
  # Build model formula dynamically
  fml <- as.formula(paste0("Surv(icu, Mortality) ~ ", var))
  
  # Fit Cox model
  model <- coxph(fml, data = data)
  
  # Extract summary info
  model_summary <- summary(model)
  
  # Extract HR, CI, and p-value (first row if multiple levels)
  tibble(
    variable = var,
    HR = model_summary$coef[1, "exp(coef)"],
    lower_CI = model_summary$conf.int[1, "lower .95"],
    upper_CI = model_summary$conf.int[1, "upper .95"],
    p_value = model_summary$coef[1, "Pr(>|z|)"]
  )
})

# Save results as text file
write.table(cox_results, "cox_pvalues.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Optional: print results to console
print(cox_results)

###################### HR forestt plot

library(survival)
library(survminer)


data <- read.delim("prevhap_db_selected.txt", header = TRUE, row.names = 1, check.names = FALSE)

data$cluster <- as.factor(data$cluster)
data$cluster <- relevel(data$cluster, ref = "C2")

# Fit Cox model
cox_model <- coxph(Surv(icu, Mortality) ~  Corticosteroid + com_age +   
                     com_pf_ratio + com_admission_type+ com_haemophilus+com_mssa+
                     com_immunosuppression+Interferon, data = data)
summary(cox_model)

# Forest plot
ggforest(cox_model, data = data, main = "Hazard Ratios by Cluster and Clinical Variables",
         cpositions = c(0.02, 0.22, 0.4),
         fontsize = 2,
         refLabel = "Reference", noDigits = 2)


################ Virus-host ratio VHR


# Load necessary library
library(dplyr)

# === 1. Load data ===
# Replace with your actual filenames
viral <- read.table("dRPKM_final_gen.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
bacterial <- read.table("bac_genera.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# === 2. Ensure same genera and sample order ===
common_genera <- intersect(rownames(viral), rownames(bacterial))
common_samples <- intersect(colnames(viral), colnames(bacterial))

viral <- viral[common_genera, common_samples, drop = FALSE]
bacterial <- bacterial[common_genera, common_samples, drop = FALSE]

# === 3. Log10(x + 1) transformation ===
viral_log <- log10(viral + 1)
bacterial_log <- log10(bacterial + 1)

# === 4. Calculate Virus-Host Ratio (VHR) ===
VHR <- viral_log / bacterial_log

# Optional: replace NaN or Inf (if bacterial abundance is zero)
VHR[!is.finite(VHR)] <- NA

# === 5. Write output to text file ===
write.table(VHR, file = "Virus_Host_Ratio.txt", sep = "\t", quote = FALSE, col.names = NA)

# Print confirmation
cat("✅ Virus-Host Ratio matrix saved as 'Virus_Host_Ratio.txt'\n")



################### Virome bubbles 

# Load libraries
library(tidyverse)
library(igraph)
library(ggraph)
library(ggsci)
library(ggrepel)
library(tidygraph)

# Read your dataset
data <- read.table("Relab_allvirome_C1.txt", sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors = FALSE)

# ---- We already have Category, no need to overwrite ----
# Just ensure it's a factor for ordering/colors
data$Category <- factor(data$Category)

# ---- Create hierarchical edges ----
# Category → Class
edges_cat_class <- data %>%
  select(Category, Class) %>%
  distinct() %>%
  mutate(from = paste0("CAT_", Category),
         to   = paste0("CLS_", Class)) %>%
  select(from, to)

# Class → Family
edges_class_fam <- data %>%
  select(Class, family, score) %>%
  mutate(from = paste0("CLS_", Class),
         to   = paste0("FAM_", family)) %>%
  select(from, to, score)

# Combine edges
edges <- bind_rows(edges_cat_class, edges_class_fam)

# ---- Create vertex list ----
vertices <- data.frame(name = unique(c(edges$from, edges$to)))

# ---- Assign hierarchy levels ----
vertices <- vertices %>%
  mutate(level = case_when(
    grepl("^CAT_", name) ~ "Category",
    grepl("^CLS_", name) ~ "Class",
    grepl("^FAM_", name) ~ "Family"
  ))

# ---- Compute node sizes ----
# Family size: their score
fam_size <- data %>%
  mutate(name = paste0("FAM_", family)) %>%
  select(name, size = score)

# Class size: sum of scores of its families
class_size <- data %>%
  group_by(Class) %>%
  summarise(size = sum(score, na.rm = TRUE)) %>%
  mutate(name = paste0("CLS_", Class))

# Category size: sum of scores of all its classes
cat_size <- data %>%
  group_by(Category) %>%
  summarise(size = sum(score, na.rm = TRUE)) %>%
  mutate(name = paste0("CAT_", Category))

# Merge all sizes
vertex_sizes <- bind_rows(cat_size, class_size, fam_size)
vertices <- vertices %>%
  left_join(vertex_sizes, by = "name") %>%
  mutate(size = ifelse(is.na(size), 0.5, size))

# ---- Color mapping (NPG palette by Category) ----
unique_cats <- unique(data$Category)
npg_colors <- rep(ggsci::pal_npg("nrc")(10), length.out = length(unique_cats))
cat_colors <- setNames(npg_colors, unique_cats)

custom_colors <- c(
  phage = "#e31a1c",   # blue
  euka = "purple",    # green
  Giant = "yellow3",   # red
  other =  "grey"    # purple
)

vertices <- vertices %>%
  mutate(
    base_name = sub("^(CAT_|CLS_|FAM_)", "", name),
    parent_cat = case_when(
      level == "Category" ~ base_name,
      level == "Class" ~ data$Category[match(base_name, data$Class)],
      level == "Family" ~ data$Category[match(base_name, data$family)]
    )
  )

vertices <- vertices %>% distinct(name, .keep_all = TRUE)
# ---- Build the graph ----
graph <- graph_from_data_frame(edges, vertices = vertices)

pdf(file="C2_bubblesV2.pdf", height = 10, width = 10)
# ---- Plot ----
ggraph(graph, layout = 'circlepack', weight = size) +
  geom_node_circle(aes(fill = parent_cat), color = "grey30", alpha = 0.9) +
  
  # Labels
  geom_node_text(aes(label = ifelse(level == "Category", sub("^CAT_", "", name), "")),
                 size = 6, fontface = "bold", color = "black", repel = TRUE) +
  geom_node_text(aes(label = ifelse(level == "Class", sub("^CLS_", "", name), "")),
                 size = 4, fontface = "bold", color = "black", alpha = 0.9, repel = TRUE) +
  geom_node_text(aes(label = ifelse(level == "Family", sub("^FAM_", "", name), "")),
                 size = 2.5, color = "black", alpha = 0.85, repel = TRUE) +
  
  scale_fill_viridis_d() +  # or use scale_fill_brewer(palette = "Set3") or any other palette
  
  theme_void() +
  coord_equal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "right"
  )
dev.off()


############################## alluvial clusters

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyverse)
library(viridis)
library(RColorBrewer)

data <- read.delim(file="meta.txt", header=T, check.names=F)

head(data)

is_lodes_form(data, key = "cluster", value = "DAY", id = "SAMPLEID")

custom_colors <- brewer.pal(8, "Set2") 
#custom_colors <- colorRampPalette(brewer.pal(8, "Set2"))(48)
custom_colors <- c("indianred1","limegreen","gold","blue") 

p <-ggplot(data, aes(alluvium = Patient, x = DAY, stratum = cluster)) + 
  geom_alluvium(color = "white", aes(fill = cluster)) +
  geom_stratum(color = "white", aes(fill = cluster), width = 1) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6, discern = TRUE) +
  scale_y_discrete() +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(axis.text = element_text(size = 20), legend.position = "none")
p

pdf(file="alluvial_clusters.pdf", height = 5, width = 10)
p
dev.off()

#####


#### boxplot WILCOX NORMAL

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(dplyr)

dat <- read.delim("correl_clusters_nohap_without_outlier_K.txt", stringsAsFactors = FALSE)

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
  get_summary_stats(Streptococcus, type = "median_iqr")

stat.test <- dat %>% 
  wilcox_test(Streptococcus ~ cluster) %>%
  add_significance()
stat.test
dat %>% wilcox_effsize(Streptococcus ~ cluster)
stat.test <- stat.test %>% add_xy_position(x = "cluster")

p <-ggplot(dat, aes(cluster, Streptococcus)) +  
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


# === Compute IQR lines for WBC ===
iqr_data <- dat %>%
  group_by(cluster) %>%
  summarise(Q1 = quantile(Streptococcus, 0.25),
            Q3 = quantile(Streptococcus, 0.75),
            xnum = as.numeric(factor(cluster)))

# === 1. Plot boxplots with same format as Avg_RPKM script ===
p <- ggplot(dat, aes(x = cluster, y = Streptococcus, fill = cluster)) + 
  
  # Boxplot layer
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", size = 2) +
  
  # Jittered points
  geom_jitter(shape = 21, color = "black", aes(fill = cluster),
              width = 0.1, size = 10, stroke = 1, alpha = 0.5) +
  
  # Statistical comparison (Wilcoxon)
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("C1_high_risk","NOHAP")),
                     label.y = max(dat$Streptococcus) * 1.05,
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
  labs(y = "Streptococcus", x = "")

# === 2. Display plot ===
print(p)

# === 3. Wilcoxon test ===
wilcox_result <- wilcox.test(Streptococcus ~ cluster, data = dat)
print(wilcox_result)

pdf("Streptococcus_boxplot.pdf",width=5,height=5);
p
dev.off()


################## Xgboost

library(xgboost)
library(dplyr)
library(readr)
library(tibble)
library(tidyverse)
library(caret)
library(Matrix)
library(pROC)

# Load your data
Bac <- read.delim("Relab_bac_species_prev_sig_FFT.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
VAL <- read.delim("Relab_bac.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Prepare VAL with outcome placeholder
VAL$C1 <- NA  

set.seed(123)

# Convert train and test C1 to factors and numeric labels
train$C1 <- factor(train$C1)
test$C1 <- factor(test$C1, levels = levels(train$C1)) # maintain level order

train_label <- as.numeric(train$C1) - 1
test_label <- as.numeric(test$C1) - 1

# Select numeric feature columns (exclude C1 and non-numeric)
train_features <- train %>% select(-C1) %>% select(where(is.numeric))
test_features <- test %>% select(-C1) %>% select(where(is.numeric))

# Convert to numeric matrices
train_matrix <- as.matrix(train_features)
test_matrix <- as.matrix(test_features)

# Prepare VAL numeric matrix - keep only features present in train set
VAL_features <- VAL %>% select(colnames(train_features)) %>% select(where(is.numeric))
VAL_matrix <- as.matrix(VAL_features)

# Train XGBoost model
model <- xgboost(
  data = train_matrix,
  label = train_label,
  nrounds = 100,
  objective = "multi:softmax",
  num_class = length(levels(train$C1)),
  verbose = FALSE
)

# Predict on test set
test_pred <- predict(model, test_matrix)

# Calculate and print AUC
roc_obj <- roc(test_label, test_pred)
auc_val <- auc(roc_obj)
cat('AUC:', round(auc_val, 3), '\n')

# Optionally plot ROC curve
plot(roc_obj)

print(table(Predicted = test_pred, Actual = test_label))

cat("\nConfusion Matrix:\n")
conf_mat <- table(Predicted = test_pred, Actual = test_label)
print(conf_mat)

if (length(levels(train$C1)) == 2) {
  # Binary classification metrics
  confusion <- confusionMatrix(
    factor(test_pred, levels = 0:1),
    factor(test_label, levels = 0:1),
    positive = "1"
  )
  
  cat("\nAccuracy:", round(confusion$overall["Accuracy"], 3), "\n")
  cat("Sensitivity:", round(confusion$byClass["Sensitivity"], 3), "\n")
  cat("Specificity:", round(confusion$byClass["Specificity"], 3), "\n")
  
} else {
  # Multi-class accuracy
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  cat("\nMulti-class Accuracy:", round(accuracy, 3), "\n")
}




################## Upset plot xgboost

library(UpSetR)
library(RColorBrewer)
library(dplyr)

# === Read files ===
viral    <- read.delim("filtered_validation_viral_all.txt", stringsAsFactors = FALSE)
bacterial <- read.delim("filtered_validation_bacterial_all.txt", stringsAsFactors = FALSE)
combined  <- read.delim("filtered_validation_MIX_all.txt", stringsAsFactors = FALSE)

# ----------------------------------------------------
# 1️⃣ Extract SAMPLEIDs belonging to cluster C1
# ----------------------------------------------------

# Viral dataset: assuming "cluster" column exists
viral_C2 <- viral %>%
  filter(cluster == "C2") %>%
  pull(SAMPLEID)

# Bacterial dataset: assuming bacterial cluster stored in fam (change if needed)
bacterial_C2 <- bacterial %>%
  filter(cluster == "C2") %>%       # update if bacterial cluster column name is different
  pull(SAMPLEID)

# Combined dataset: same structure
combined_C2 <- combined %>%
  filter(cluster == "C2") %>%       # update if needed
  pull(SAMPLEID)

# ----------------------------------------------------
# 2️⃣ Build list for UpSet plot
# ----------------------------------------------------

listInput <- list(
  Viral_FFTree     = viral_C2,
  Bacterial_FFTree = bacterial_C2,
  Combined_FFTree   = combined_C2
)

# ----------------------------------------------------
# 3️⃣ Run UpSetR
# ----------------------------------------------------

Upsetplot <- upset(fromList(listInput),
                   mb.ratio = c(0.6, 0.4), 
                   order.by = "freq",
                   nsets = 3,
                   sets.bar.color = brewer.pal(3, "Set2"),
                   decreasing = TRUE,
                   text.scale = c(2, 2, 2, 2, 2, 2),
                   point.size = 5, 
                   line.size = 1,
                   matrix.color = "black", 
                   mainbar.y.label = "frequency\n(Intersection size)",  # Fixed multi-line text
                   sets.x.label = "signature\n(Set Size)")

Upsetplot

########### tkna importance 

# ==============================
# Importance score calculation
# ==============================

# Load libraries
library(dplyr)

# Read input table (replace with your file name)
df <- read.delim("probabilities_to_randomly_find_nodes.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Clean probability column (remove "%" and convert to numeric)
df$Probability <- as.numeric(gsub("%", "", df$Probability))

# Normalize Degree
df$Degree_norm <- (df$`Observed_degree` - min(df$`Observed_degree`)) /
  (max(df$`Observed_degree`) - min(df$`Observed_degree`))

# Normalize BiBC
df$BiBC_norm <- (df$`Observed_BiBC` - min(df$`Observed_BiBC`)) /
  (max(df$`Observed_BiBC`) - min(df$`Observed_BiBC`))

# Normalize Probability (invert, so lower prob = higher importance)
df$Prob_norm <- 1 - ((df$Probability - min(df$Probability)) /
                       (max(df$Probability) - min(df$Probability)))

# Compute importance score (average of 3 criteria)
df$Importance_score <- rowMeans(df[, c("Degree_norm", "BiBC_norm", "Prob_norm")])

# Select and reorder columns (base R way)
df_out <- df[, c("Node", "Observed_degree", "Observed_BiBC", "Probability",
                 "Degree_norm", "BiBC_norm", "Prob_norm", "Importance_score")]
# Export result to text file
write.table(df_out,
            file = "importance_scores.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Importance scores calculated and saved to 'importance_scores.txt'\n")

