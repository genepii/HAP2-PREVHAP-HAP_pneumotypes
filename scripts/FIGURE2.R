############################################################################################################################
#### CODE USED TO GENERATE FIGURE 2
############################################################################################################################

#A

#Script: Full MEFISO Multi-Omics Analysis Workflow

############################################################
# 0. Load required libraries
############################################################
library(MOFA2)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library(readr)
library(ggplot2)
library(ggpubr)
library(cluster)
library(factoextra)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(magrittr)

############################################################
# 1. Preparation of tables for MEFISO
############################################################

library(reshape2)
library(dplyr)
library(ggpubr)

# Function to reshape tables
reshape_file <- function(file, view, out) {
  data <- read.delim(file, row.names = 1, check.names = FALSE)
  data <- reshape2::melt(as.matrix(data))
  colnames(data) <- c("Sample", "Family", "Value")
  data$view <- view
  write.table(data, out, sep = "\t", row.names = FALSE)
}

# 1) reshape inputs
reshape_file("bacteriome.txt", "bacrep", "bacrep.txt")
reshape_file("bacteriome_functions.txt", "bacrep", "bacrepf.txt")
reshape_file("bacterial_metabolic_pathways.txt", "bacrep", "metab.txt")
reshape_file("virome.txt", "vir", "vir.txt")
reshape_file("active_virome.txt", "vir", "Avir.txt")
reshape_file("virome_functions.txt", "vir", "virf.txt")
reshape_file("active_virome_functions.txt", "vir", "Avirf.txt")

# 2) merge all tables
files <- list(
  "bacrep.txt","Avir.txt","vir.txt",
  "bacrepf.txt","Avirf.txt","virf.txt",
  "metab.txt"
)

dt <- bind_rows(lapply(files, function(f) read.delim(f, header = TRUE, stringsAsFactors = FALSE)))

# 3) standardize columns
colnames(dt) <- c("sample", "feature", "value", "view")

# 4) transform
dt$value <- log10(dt$value + 1)
dt$value[is.infinite(dt$value)] <- NA

# 5) plot (optional)
ggdensity(dt, x = "value", fill = "view") +
  facet_wrap(~view, nrow = 1, scales = "free")

# 6) export
write.table(dt, "dt.txt", sep = "\t", row.names = FALSE)


############################################################
# 2. MOFA Model Analysis
############################################################
##A - create MOFA in a linux environment##

####create mofa object

obj  <- create_mofa(data = dt)
samples_metadata(obj) <- meta

####add timing points to mofa object

obj <- set_covariates(obj, covariates = t(meta[, "DAY", drop= F]))
plot_data_overview(obj,show_covariate = T,show_dimensions = T) 

####prepare mofa object for training

data_opts <- get_default_data_options(obj)
head(data_opts)
model_opts <- get_default_model_options(obj)
#model_opts$num_factors <- 15 # this should be higher once you add more features in your expression data
head(model_opts)
train_opts <- get_default_training_options(obj)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42
head(train_opts)
mefisto_opts <- get_default_mefisto_options(obj)
#mefisto_opts$warping <- TRUE
mefisto_opts$new_values <- matrix(-7:0, nrow =1) # set time points to interpolate factors to

obj <- prepare_mofa(obj, model_options = model_opts,
                   mefisto_options = mefisto_opts,
                   training_options = train_opts,
                   data_options = data_opts)

obj <- run_mofa(obj, use_basilisk = F,
               outfile = file.path(getwd(), "model.hdf5"))

# saveRDS(obj , "obj.rds")

plot_variance_explained(obj)
plot_variance_explained(obj, plot_total = T)[[2]]

############################################################
# 3. Clustering
############################################################
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
library(MOFA2)
library(dplyr)

dt <- read.delim("dt.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
meta <- read.delim("meta.txt", stringsAsFactors = FALSE, row.names = 1)
# Views names
unique(dt$view)
# Number of samples
length(unique(dt$sample))
# Number of features per view
dt <- as.data.table(dt)
result <- dt[, .(unique_feature_count = length(unique(feature))), by = view]
result <- dt %>%
  group_by(view) %>%
  summarise(unique_feature_count = n_distinct(feature))
result
a <- ggdensity(dt, x="value", fill="view") +
  facet_wrap(~view, nrow=1, scales = "free")
a

obj  <- create_mofa(data = dt)
meta$sample = rownames(meta)
samples_metadata(obj) <- meta

####add timing points to mofa object

obj <- set_covariates(obj, covariates = t(meta[, "DAY", drop= F]))
#pdf("data_overview.pdf", width = 20, height = 15)
plot_data_overview(obj,show_covariate = T,show_dimensions = T) 
#dev.off()

######################################
obj <- load_model("model.hdf5")
samples_metadata(obj) <- meta

plot_variance_explained(obj, plot_total = T)[[2]]
p<-plot_variance_explained(obj)
# Custom gradient for ordered factors
p + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5)
plot_factor_cor(obj)
get_scales(obj)
variance_data <- get_variance_explained(obj)
str(variance_data)
variance_factors <- variance_data$r2_per_factor
#write.table(variance_factors, file = "variance_explained.txt", sep = "\t")
data <- get_data(obj)
lapply(data, function(x) lapply(x, dim))[[1]]

factors <- get_factors(obj, as.data.frame = T)
head(factors)

# Spread the data to wide format
data_wide <- factors %>%
  spread(key = factor, value = value)

formatted_factors=t(data_wide)
formatted_factors <- formatted_factors[!rownames(formatted_factors) %in% "group", ]
write.table(formatted_factors, file = "all_formatted_factors_all.txt", sep = "\t", col.names = NA, row.names = T, quote = F)

######################################### Partitionning clustering
library(stats)
library(cluster)
library(factoextra)
set.seed(123)
matrix <- read.delim("all_formatted_factors_all - Copie.txt", header = TRUE, sep = "\t", row.names = 1, check.names = F)
meta <- read.delim("meta.txt", stringsAsFactors = FALSE)
head(matrix)
#matrix=log10(matrix+0.001)

fviz_nbclust(t(matrix), kmeans, method = "wss") +
  geom_vline(xintercept = "", linetype = 3)+
  labs(subtitle = "Elbow method")

fviz_nbclust(t(matrix), kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette Method")

fviz_nbclust(t(matrix), kmeans, method = "gap_stat") +
  labs(subtitle = "Gap_stat Method")


gap_stat <- clusGap(t(matrix), FUN = kmeans, nstart = 25, K.max = 10, B = 500)
fviz_gap_stat(gap_stat)


best_k <- maxSE(gap_stat$Tab[,"gap"], gap_stat$Tab[,"SE.sim"], method = "firstSEmax")
cat("Best number of clusters according to gap statistic:", best_k, "\n")
set.seed(123)  # reproducibility
clust <- kmeans(t(matrix), centers = best_k, nstart = 25)
#cluster_df <- data.frame(cluster = clust$cluster)
#write.table(cluster_df, "cluster_assignments.txt", row.names = TRUE, col.names = NA, sep = "\t")

optimal_k <- 2 # Replace this with the actual number you get from the methods above

# Perform k-means clustering
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(t(matrix), centers = optimal_k, nstart = 25)
# View the results
print(kmeans_result)

dim(t(matrix))

#pdf("plot_mofa_clusters_noHAP.pdf", width = 7.5, height = 5)
famd_clust_a <- fviz_cluster(kmeans_result,data = t(matrix),
                             show.clust.cent = F,
                             ggtheme = theme_classic(),
                             main = "", geom = "point", palette = c("indianred1", "limegreen", "green","yellow4","blue","gold","grey"))
famd_clust_a 
#dev.off()

famd_clust_a <- fviz_cluster(kmeans_result, 
                             data = t(matrix),
                             show.clust.cent = FALSE,
                             geom = "point",
                             pointsize = 8,  # Increased from default (1-2)
                             stroke = 1,     # Point border thickness
                             alpha = 0.5,     # Solid circle shape
                             palette = c("indianred1", "limegreen", "green", "yellow4", "blue", "gold", "grey"),
                             main = "",
                             ggtheme = theme_classic(base_size = 20)) +  # Base font size
  
  # Customize for publication quality
  theme(
    # Axis and labels
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    
    # Title and legend
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "none",
    legend.key = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    # Plot margins
    plot.margin = margin(15, 15, 15, 15, "pt"),
    
    # Remove scientific notation on axes if needed
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(margin = margin(r = 5))
  ) +
  
  # Additional customization
  labs(
    x = "MEFISTO Factor 1",
    y = "MEFISTO Factor 2"
  )

# Display the plot
print(famd_clust_a)

#pdf("plot_mofa_clusters_HAPnew_prev50.pdf", width = 5, height = 5)
famd_clust_a
#dev.off()
################
p<-fviz_cluster(kmeans_result, data = t(matrix), geom = "point", shape = "cluster", pointsize = 3,
             ellipse.type = "t",
             palette = c("indianred1", "limegreen", "blue", "green","yellow4","blue","gold","grey") 
) + 
  theme_bw(base_size = 14) + 
  theme(
    legend.position = "right", 
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),  
    axis.title = element_text(face = "bold", size = 14), 
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.5), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )

p

kmeans_result_df <- data.frame(
  SAMPLEID = meta[, 1],                 # first column of meta = IDs
  cluster  = kmeans_result$cluster      # cluster assignments
)
#write.table(kmeans_result_df, "cluster_assignments.txt", row.names = TRUE, col.names = NA, sep = "\t")

# Perform the join explicitly on sampleid
joined_data <- dplyr::left_join(meta, kmeans_result_df, by = "SAMPLEID")

datat<- joined_data


#B

################# reformat meta_patients for KP analysis

library(dplyr)
library(stringr)

# Your data frame
datat = read.delim(file="metadata.txt",header=T,check.names=F)

# Extract numeric day
datat <- datat %>%
  mutate(DAY_numeric = as.numeric(str_extract(DAY, "\\d+")))

# Filter for max day per patient
df_max_day <- datat %>%
  group_by(Patient) %>%
  filter(DAY_numeric == max(DAY_numeric, na.rm = TRUE)) %>%
  ungroup()

# View result
print(df_max_day)

#write.table(df_max_day, file = "filtered_data.txt", sep = "\t", quote = FALSE,col.names = NA)

datat<- df_max_day

library(survival)
library(survminer)

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

#C

###################### HR forestt plot

library(survival)
library(survminer)

data <- read.delim("metadata_clinical.txt", header = TRUE, row.names = 1, check.names = FALSE)

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



#D DEG analysis and volcano plot

############## nbglm volcano gene symbol
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(cowplot)
library(dplyr)

# Read expression data and metadata
expression_df <- read.delim("data/prevhap_counts_human_coding_genes.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
metadata <- read.delim("metadata/metadata.txt", stringsAsFactors = FALSE)


sapply(metadata$cluster,function(x)ifelse(x=="C1",1,0))->metadata$C1
sapply(metadata$cluster,function(x)ifelse(x=="C2",1,0))->metadata$C2

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

### NEW STEP: Subset to inflammation genes
infl_genes <- read.delim("data/infla_list.txt", header = FALSE, stringsAsFactors = FALSE)
infl_genes <- infl_genes[[1]]  

# Keep only inflammation genes present in the dataset
keep_genes <- rownames(y) %in% infl_genes
y <- y[keep_genes, , keep.lib.sizes = FALSE]

y <- y[rowSums(y$counts) > 0, , keep.lib.sizes = FALSE]
y <- y[, colSums(y$counts) > 0, keep.lib.sizes = FALSE]

metadata <- metadata[metadata$SAMPLEID %in% colnames(y), ]
metadata <- metadata[match(colnames(y), metadata$SAMPLEID), ]


design <- model.matrix(~ C1 * IFN_dose * DAY + IFN_dose + Patient + DAY, data = metadata)

fit <- glmQLFit(y, design, robust = TRUE)

colnames(design)

qlt <- glmQLFTest(fit, coef = "C1:IFN_dose:DAY" )
#qlt <- glmQLFTest(fit, coef = "C1:IFN_dose" )


# Get top genes (restricted to inflammation list)
topgenes <- topTags(qlt, n = nrow(y))

# Identify enriched and depleted genes
tol10b.enriched <- topgenes$table$logFC > 0 & topgenes$table$PValue < 0.05
tol10b.depleted <- topgenes$table$logFC < 0 & topgenes$table$PValue < 0.05

tg_df <- topgenes$table
tg_df$Ensembl <- rownames(tg_df)
ens_map <- read.delim("data/ens_to_symbol.txt", header = TRUE, stringsAsFactors = FALSE)
tg_df <- merge(tg_df, ens_map, by = "Ensembl", all.x = TRUE)
tg_df$Label <- ifelse(is.na(tg_df$Symbol) | tg_df$Symbol == "", tg_df$Ensembl, tg_df$Symbol)
rownames(tg_df) <- tg_df$Ensembl

write.table(tg_df, file = "inputList.txt", sep = "\t", col.names = NA, quote = FALSE)

##################################################################

volcano <- ggplot(tg_df, aes(x = logFC, y = -10 * log10(PValue))) +
  geom_point(
    alpha = 0.6,
    size = 4,
    aes(color = case_when(
      PValue < 0.05 & logFC > 0  ~ "Up",
      PValue < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    ))
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = ifelse(PValue < 0.05 & abs(logFC) > 0, Label, "")),
    size = 5,
    max.overlaps = 20
  ) +
  color = "black", linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Up" = "indianred1",
    "Down" = "limegreen",
    "NS" = "grey"
  )) +
  
  theme_bw() +
  labs(x = "Log2 Fold Change", y = "-10 * log10(PValue)") +
  ggtitle("") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none",
    legend.title = element_blank()
  )
volcano

pdf("Volcano_plot.pdf",width=8,height=8)
volcano
dev.off()


#E

########## GO enrichment across clusters (C1, C2, N1, N2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Read file
genes_df <- read.delim("data/inputList.txt", header = TRUE, check.names = FALSE)

head(genes_df)

# Split gene lists by cluster
genes_C1 <- genes_df$ID[genes_df$Group == "C1"]
genes_C2 <- genes_df$ID[genes_df$Group == "C2"]


# Define a function for GO enrichment
run_enrichGO <- function(gene_list, group_name) {
  ego <- enrichGO(gene         = gene_list,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = 'ENSEMBL',
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05)
  ego@result$Group <- group_name
  return(ego)
}

# Run enrichment for each cluster
ego_C1 <- run_enrichGO(genes_C1, "C1")
ego_C2 <- run_enrichGO(genes_C2, "C2")

# Combine all results
ego_combined <- ego_C1
ego_combined@result <- rbind(ego_C1@result, ego_C2@result)

# Create combined dataframe for plotting
ego_df <- as.data.frame(ego_combined@result)

# Keep top 100 terms per group (by adjusted p-value)
ego_df <- ego_df %>%
  group_by(Group) %>%
  top_n(-100, p.adjust) %>%
  ungroup()

write.csv(
  ego_df,
  file = "GO_combined_results.csv"
)


##dotplot
library(readr)
library(tidyr)

go_data <- read_csv("data/GO_combined_results.csv")

# Convert GeneRatio from "11/20" to numeric
go_data <- go_data %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDen"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatioValue = GeneNum / GeneDen)

# Ensure 'Description' is a factor so that it appears in order
go_data$Description <- factor(go_data$Description, levels = rev(unique(go_data$Description)))

group_order <- c("C2", "C1")  # specify your desired cluster order

go_data <- go_data %>%
  # First, make Group a factor with desired order
  mutate(Group = factor(Group, levels = group_order)) %>%
  # Then, order Description within each group by Count (or other variable)
  arrange(Group, desc(Count)) %>%
  # Make Description a factor with the order preserved
  mutate(Description = factor(Description, levels = unique(Description)))

low_color <- "red"
high_color <- "blue"
mid_point <- median(go_data$p.adjust)

go_data <- go_data %>%
  arrange(desc(GeneRatioValue)) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

pdf("GO_clusters.pdf", width = 20, height = 15)
# Create the dot plot with numeric GeneRatio on x-axis
library(dplyr)

top_go <- go_data %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>%
  slice_head(n = 40)

ggplot(top_go,
       aes(x = GeneRatioValue, y = reorder(Description, GeneRatioValue))) +
  geom_point(aes(size = Count, fill = p.adjust),
             shape = 21,
             color = "black",
             stroke = 0.5,
             alpha = 0.7) +
  scale_fill_gradient(low = low_color, high = high_color) +
  scale_size(range = c(5, 15)) +
  facet_grid(. ~ Group, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 18, colour = "black")
  ) +
  labs(
    x = "Gene Ratio",
    y = NULL,
    fill = "Adjusted p-value",
    size = "Gene Count"
  )
dev.off()


