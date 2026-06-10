############################################################################################################################
#### CODE USED TO GENERATE Supplementary figures
############################################################################################################################

# Supplementary FIGURE 1

library(DiagrammeR)

flow <-grViz("
digraph flow {
 graph [rankdir =  TB, splines = ortho]
  node [shape=box, style=filled,fontname = Helvetica,
    fontsize = 14]

  A [label='184 total samples (N) / 94 patients (n) (PREVHAP)', fillcolor='white']
  B [label='184 samples / 94 patients, sequenced samples\n with vmetaG', fillcolor='#5B9BD5']
  B1 [label='133/77 vmetaG passed IC\n threshold (10 reads MS2)\n and rarefaction depth (10k RPKM)']
  B2 [label='60/35 HAP samples']
  B3 [label='73/42 no HAP samples']
  C [label='122/53 patients, sequenced samples\n with metaT', fillcolor='#FFC000']
  C1 [label='116/53 metaT passed IC\n threshold (10 reads MS2)\n and rarefaction depth (10k RPKM)']
  C2 [label='115/53 metaT passed Host\n threshold\n (5 log10 coding reads)']
  C3 [label='55/22 HAP samples']
  C4 [label='61/31 no HAP samples']
  C5 [label='51/22 HAP samples']
  C6 [label='64/31 no HAP samples']
  C7 [label='42/22 before HAP\n onset samples']
  C8 [label='13/12 After HAP\n onset samples']
  B4 [label='44/31 before HAP\n onset samples']
  B5 [label='16/15 After HAP\n onset samples']
  C9 [label='38/22 before HAP\n onset samples']
  C10 [label='13/12 After HAP\n onset samples']
  
  
  A -> B
  A -> C
  C -> C1
  C -> C2
  B -> B1
  B1 -> B2
  B1 -> B3
  C1 -> C3
  C1 -> C4
  C2 -> C5
  C2 -> C6
  C3 -> C7
  C3 -> C8
  C5 -> C9
  C5 -> C10
  B2 -> B4
  B2 -> B5

}
")


flow

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
svg <- export_svg(flow)
writeLines(svg, "flow.svg")

rsvg_pdf("flow.svg", "flow_diagram.pdf", height = 14, width = 20)


# Supplementary FIGURE 2-3

library(dplyr)
library(grid)  # needed for arrow()
data <- read.delim("metadata.txt", header = TRUE, stringsAsFactors = FALSE)

data <- data %>%
  arrange(IFN, HAP_day) %>%
  mutate(
    Patient = factor(Patient, levels = unique(Patient))
  )
ifn_points <- expand.grid(
  Patient = unique(data$Patient[data$IFN == "G"]),
  DAY = c(1, 3, 5, 7)
)

last_days <- data %>%
  filter(DAY %in% c(0, 3, 7)) %>%
  group_by(Patient) %>%
  filter(DAY == max(DAY, na.rm = TRUE)) %>%
  ungroup()
line_data <- data %>%
  select(Patient, DAY) %>%
  mutate(type = "DAY")

hap_links <- data %>%
  filter(HAP_day > 0) %>%
  transmute(
    Patient,
    DAY = HAP_day,
    type = "HAP"
  )

line_data <- bind_rows(line_data, hap_links) %>%
  arrange(Patient, DAY)

pdf('timeline_PREVHAPcluster_subset_pneumonia_final.pdf', height = 20, width = 10)
ggplot(data, aes(y = Patient)) +
  
  geom_line(
    data = line_data,
    aes(x = DAY, y = Patient, group = Patient),
    color = "grey"
  )+
  
  # DAY baseline points
  geom_point(
    aes(
      x = DAY
    ),
    size = 3,
    color = "black",
    shape = 16
  ) +
  geom_point(
    data = last_days,
    aes(
      x = DAY,
      y = Patient,
      color = icuo
    ),
    size = 4
  )+
  geom_point(
    data = subset(data, HAP_day > 0),
    aes(x = HAP_day, y = Patient, shape = "HAP day"),
    size = 4,
    color = "red"
  )+
  geom_point(
    data = ifn_points,
    aes(
      x = DAY,
      y = Patient,
      shape = "IFN dose"
    ),
    size = 4,
    color = "green4"
  )+
  # Continuous ICU color scale
  scale_color_gradient(
    low = "blue",
    high = "red",
    name = "ICU extubation day",
    limits = c(0, 80),
    oob = scales::squish
  ) +
  
  # Shapes for IFN injections
  scale_shape_manual(
    values = c(
      "0" = 16,
      "2" = 17,
      "5" = 15,
      "HAP day" = 8,
      "IFN dose" = 25
    ),
    name = "IFN injections"
  ) +
  
  # X axis
  scale_x_continuous(
    limits = c(min(data$DAY, na.rm = TRUE), 14),
    breaks = seq(
      floor(min(c(data$DAY, data$HAP_day), na.rm = TRUE)),
      14,
      by = 1
    )
  ) +
  
  theme_classic() +
  
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    
    legend.position = "top",
    
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
  ) +
  
  labs(
    x = "Study enrollment day",
    y = "Patient ID",
    shape = "IFN injections",
    color = "ICU extubation day"
  )
dev.off()

### new timeline by hap onset :

data <- read.delim("metadata.txt", header = TRUE, stringsAsFactors = FALSE)

library(grid)  # needed for arrow()

pdf('timeline_prevhap_subset_pneumonia_finalonset.pdf', height = 20, width = 10)

ggplot(data, aes(y = Patient)) +
  
  # Pre-HAP timeline
  geom_segment(aes(x = -7, xend = 0, yend = Patient,
                   color = "Pre-HAP timeline (referred to pneumotypes)"),
               size = 1) +
  
  # Post-HAP timeline
  geom_segment(aes(x = 0, xend = 6, yend = Patient,
                   color = "Post-HAP follow-up (antibiotic exposure)"),
               size = 1) +
  
  # Vertical line
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "red", size = 2) +
  
  # Label
  annotate("text", x = 0, y = Inf, label = "HAP day",
           vjust = 1, hjust = 0.3, size = 6, color = "red") +
  
  # Points (base)
  geom_point(aes(x = DAY, shape = IFN),
             size = 4, color = "black") +
  
  # Points with IFN color
  geom_point(aes(x = DAY, shape = IFN, color = as.factor(IFN_injections)),
             size = 4) +
  
  # COLOR SCALE (now includes BOTH timeline + IFN)
  scale_color_manual(values = c(
    "0" = "blue",
    "2" = "gold3",
    "5" = "green"
  )) +
  
  scale_x_continuous(
    limits = c(-7,0),
    breaks = seq(-7, 0, by = 1)
  ) +
  
  theme_classic() +
  
  labs(
    x = "Days related to HAP onset",
    y = "Patient ID",
    shape = "IFN",
    color = "Legend"
  )+ theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    
    legend.position = "top",
    
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
  )


dev.off()

# Supplementary FIGURE 4

library("phyloseq")

library("ggplot2")

otu = read.delim(file="bacteriome.txt",header=T,row.names=1,check.names=F) 
meta = read.delim(file="meta.txt",header=T,row.names=1,check.names=F) 

tax <- data.frame(
  Phylum = rep("_", nrow(otu)),
  Class = rep("_", nrow(otu)),
  Order = rep("_", nrow(otu)),
  Family = rep("_", nrow(otu)),
  Genus = rep("_", nrow(otu)),
  Species = rep("_", nrow(otu))
)
row.names(tax) <- row.names(otu)
#otu=round(otu, digits = 0)

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))

physeq = phyloseq(OTU, TAX)

physeq = merge_phyloseq(physeq, sample_data(meta))
###############################
otu=round(otu, digits = 0)

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))

physeq = phyloseq(OTU, TAX)

physeq = merge_phyloseq(physeq, sample_data(meta))

library("cluster"); packageVersion("cluster")
library(ggplot2)

theme_set(theme_bw())
exord = ordinate(physeq, method="MDS", distance="jsd")

pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
x = phyloseq:::scores.pcoa(exord, display="sites")
# gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn = clusGap(x[, 1:2], FUN=pam1, K.max = 6, B = 50)


gap_statistic_ordination = function(ord, FUNcluster, type="sites", K.max=6, axes=c(1:2), B=500, verbose=interactive(), ...){
  require("cluster")
  #   If "pam1" was chosen, use this internally defined call to pam
  if(FUNcluster == "pam1"){
    FUNcluster = function(x,k) list(cluster = pam(x, k, cluster.only=TRUE))     
  }
  # Use the scores function to get the ordination coordinates
  x = phyloseq:::scores.pcoa(ord, display=type)
  #   If axes not explicitly defined (NULL), then use all of them
  if(is.null(axes)){axes = 1:ncol(x)}
  #   Finally, perform, and return, the gap statistic calculation using cluster::clusGap  
  clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
}


plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
  require("ggplot2")
  gstab = data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab))
  p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
  p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
  p = p + ggtitle(title)
  return(p)
}


gs = gap_statistic_ordination(exord, "pam1", B=50, verbose=FALSE)
print(gs, method="Tibs2001SEmax")

plot_clusgap(gs)


plot(gs, main = "Gap statistic for the 'Enterotypes' data")
mtext("Looks like 4 clusters is best, with 3 and 5 close runners up.") 



### rarefaction

require(parallel)
options(mc.cores= 2)
require(vegan)
## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = "blue", plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


p <- ggrare(physeq, step = 50, color = "GROUP", se = FALSE) 

pdf("rarefaction_curves_bacteriome_all.pdf", width = 30, height = 20)
p +
  facet_wrap(~GROUP, scales = "free") +
  theme_bw() +
  geom_vline(
    xintercept = 10000,
    color = "blue",
    linetype = "dashed",
    linewidth = 2
  )
dev.off()



# Supplementary FIGURE 5

library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
library(MOFA2)
library(dplyr)
library(cowplot)

dt <- read.delim("data/dt.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

meta <- read.delim("metadata.txt", stringsAsFactors = FALSE, row.names = 1)


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
#pdf("density.pdf")
a
#dev.off()

colnames(meta)

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


# Supplementary FIGURE 6

library(stats)
library(cluster)
library(factoextra)
set.seed(123)
matrix <- read.delim("all_formatted_factors.txt", header = TRUE, sep = "\t", row.names = 1, check.names = F)
meta <- read.delim("meta.txt", stringsAsFactors = FALSE)
head(matrix)

fviz_nbclust(t(matrix), kmeans, method = "wss") +
  geom_vline(xintercept = "", linetype = 3)+
  labs(subtitle = "Elbow method")

fviz_nbclust(t(matrix), kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette Method")

gap_stat <- clusGap(t(matrix), FUN = kmeans, nstart = 25, K.max = 10, B = 500)
fviz_gap_stat(gap_stat)


###### cluster partionning 
library(ggplot2)
library(tidyr)

df <- data.frame(
  Prevalence = c(0, 10, 20, 30, 40, 50),
  Elbow = c(2, 2, 2, 2, 2, 2),
  Silhouette = c(10, 5, 2, 5, 7, 2),
  Gap_stat = c(2, 2, 1, 4, 2, 2)
)

df_long <- pivot_longer(df, cols = -Prevalence,
                        names_to = "Method",
                        values_to = "Clusters")

pdf("clustering_prevalence.pdf",width=5,height=3)

ggplot(df_long, aes(Prevalence, Clusters, color = Method, shape = Method)) +
  geom_point(size = 5, alpha = 0.7) +
  geom_line(aes(group = Method), linewidth = 0.8) +
  scale_y_continuous(breaks = seq(1, max(df_long$Clusters), by = 1)) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  ) +
  labs(
    x = "Prevalence (%)",
    y = "Number of Clusters",
    color = "Method",
    shape = "Method"
  )

dev.off()

# Supplementary FIGURE 7


library(tidyr)
library(dplyr)

# read your data
df <- read.delim("dt.txt", sep = "\t", stringsAsFactors = FALSE)

# convert to wide matrix
#wide_df <- df %>%
select(sample, feature, value) %>%
  pivot_wider(names_from = feature, values_from = value, values_fill = 0)

wide_df <- df %>%
  select(sample, feature, value) %>%
  group_by(sample, feature) %>%
  summarise(value = mean(value), .groups = "drop") %>%  # or sum(value)
  pivot_wider(names_from = feature, values_from = value, values_fill = 0)

# convert to matrix (samples as rownames)
mat <- as.data.frame(wide_df)
rownames(mat) <- mat$sample
mat$sample <- NULL

mat <- as.matrix(mat)

write.table(wide_df, "wide_df.txt", row.names = TRUE, col.names = NA, sep = "\t")

filecoinf <- read.delim("wide_df.txt", sep = "\t", stringsAsFactors = FALSE)

qualitative_vars <- c("death","interferon")

continuous_vars <- setdiff(colnames(filecoinf), qualitative_vars)

continuous_vars <- continuous_vars[sapply(filecoinf[continuous_vars], function(x) {
  all(grepl("^[0-9.]+$", na.omit(x)))
})]

library(effsize)
library(dplyr)
library(dplyr)
library(qwraps2)
library(ggplot2)


compute_all_features <- function(data, cluster_col, continuous_vars, qualitative_vars) {
  
  clusters <- unique(data[[cluster_col]])
  results_list <- list()
  
  ### =========================
  ### 1. VARIABLES CONTINUES
  ### =========================
  
  for (var in continuous_vars) {
    
    for (cl in clusters) {
      
      data_ref <- na.omit(data[data[[cluster_col]] == cl, var])
      data_others <- na.omit(data[data[[cluster_col]] != cl, var])
      data_global <- na.omit(data[[var]])
      
      if (length(data_ref) > 1 & length(data_others) > 1) {
        
        test_res <- wilcox.test(data_ref, data_others, exact = FALSE)
        eff <- cliff.delta(data_ref, data_others)$estimate
        
        med_ref <- median(data_ref)
        med_others <- median(data_others)
        
        direction <- ifelse(med_ref > med_others, "higher", "lower")
        median_diff <- med_ref - med_others
        
        z_score <- (mean(data_ref) - mean(data_global)) / sd(data_global)
        
      } else {
        test_res <- list(p.value = NA)
        eff <- NA
        direction <- NA
        median_diff <- NA
        z_score <- NA
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        variable = var,
        type = "continuous",
        cluster = cl,
        p_value = test_res$p.value,
        effect_size = eff,
        direction = direction,
        value_cluster = med_ref,
        value_others = med_others,
        z_score = z_score
      )
    }
  }
  
  ### =========================
  ### 2. VARIABLES QUALITATIVES
  ### =========================
  
  for (var in qualitative_vars) {
    
    for (cl in clusters) {
      
      # Table de contingence
      tab <- table(data[[var]], data[[cluster_col]] == cl)
      
      if (all(dim(tab) >= 2)) {
        
        test_res <- fisher.test(tab)
        
        # Proportions
        prop_cluster <- prop.table(tab, 2)[, "TRUE"]
        prop_others <- prop.table(tab, 2)[, "FALSE"]
        
        # On prend la catégorie dominante
        max_cat <- names(which.max(prop_cluster))
        
        direction <- paste0("over:", max_cat)
        
        effect <- test_res$estimate  # odds ratio
        
      } else {
        test_res <- list(p.value = NA)
        effect <- NA
        direction <- NA
        prop_cluster <- NA
        prop_others <- NA
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        variable = var,
        type = "categorical",
        cluster = cl,
        p_value = test_res$p.value,
        effect_size = effect,
        direction = direction,
        value_cluster = NA,
        value_others = NA,
        z_score = NA
      )
    }
  }
  
  ### =========================
  ### 3. COMBINAISON
  ### =========================
  
  results <- bind_rows(results_list)
  
  # Correction multiple par variable
  results <- results %>%
    group_by(variable) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
    ungroup()
  
  return(results)
}
results_all <- compute_all_features(
  data = filecoinf,
  cluster_col = "cluster",
  continuous_vars = continuous_vars,
  qualitative_vars = qualitative_vars
)

write.table(results_all, "results_all.txt", row.names = TRUE, col.names = NA, sep = "\t")


####plot
# Load libraries
library(tidyverse)

df <- read.delim("results_all.txt")

# Create dataframe (copy your table here)


# Create significance label (based on adjusted p-value)
df <- df %>%
  mutate(stars = case_when(
    p_value < 0.0001 ~ "****",
    p_value < 0.001  ~ "***",
    p_value < 0.01   ~ "**",
    p_value < 0.05   ~ "*",
    TRUE ~ ""
  ))

# Sort within each cluster by decreasing effect size
df <- df %>%
  group_by(cluster) %>%
  arrange(desc(effect_size), .by_group = TRUE) %>%
  mutate(variable = factor(variable, levels = unique(variable))) %>%
  ungroup()

custom_colors <- c("P2_high_risk" = "indianred1", "P1_low_risk" = "limegreen")

# Plot
p <- ggplot(df, aes(x = effect_size,
                    y = reorder(variable, effect_size),
                    fill = cluster)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = stars),
            hjust = 0,
            size = 12) +
  facet_grid(~ cluster, scales = "free") +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  labs(
    x = "Effect Size",
    y = "Variable",
    fill = "Cluster"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12, colour = "black")
  )

p

pdf("cluster_indcator2.pdf",width=10,height=7);
p
dev.off()


###### shap
library(dplyr)
library(vip)
library(randomForest)

feat  <- read.delim("wide_df.txt", stringsAsFactors = FALSE)
meta  <- read.delim("meta.txt", stringsAsFactors = FALSE)

set.seed(123)

y <- as.factor(meta$cluster)  

model <- randomForest(
  x = feat,
  y = y,
  ntree = 500,
  importance = TRUE
)

pfun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "c1"]
}

set.seed(123)

shap_imp <- vi_shap(
  model,
  train = feat,
  feature_names = colnames(feat),
  pred_wrapper = pfun,
  nsim = 30
)

print(shap_imp)

library(ggplot2)

shap_imp %>%
  arrange(desc(Importance)) %>%
  head(20) %>%
  ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top SHAP Feature Importance")

write.table(shap_imp, "shap_imp.txt", row.names = TRUE, col.names = NA, sep = "\t")



# Supplementary FIGURE 8

############################## alluvial clusters

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyverse)
library(viridis)
library(RColorBrewer)

data <- read.delim(file="meta.txt", header=T, check.names=F)

head(data)

is_lodes_form(data, key = "cluster", value = "Inclusion_ICU_Day", id = "SAMPLEID")

custom_colors <- brewer.pal(8, "Set2") 
#custom_colors <- colorRampPalette(brewer.pal(8, "Set2"))(48)
custom_colors <- c("limegreen","indianred1","gold","blue") 

p <-ggplot(data, aes(alluvium = Patient, x = Inclusion_ICU_Day, stratum = cluster)) + 
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


# Supplementary FIGURE 9



library(dplyr)
library(survival)
library(survminer)

datat <- read.delim("metadata", header = TRUE, row.names = 1, check.names = FALSE)

dataset_surv <- datat %>%
  mutate(
    event_extub = if_else(ventilation_stop == "yes", 1, 0),
    
    # temps = temps jusqu’à arrêt OU jusqu’à décès/censure
    time_extub = icuo,
    cluster = factor(cluster),
    cluster = relevel(cluster, ref = "N1")
  )

# Modèle de Kaplan-Meier : temps jusqu’à extubation avec censure des décès
km_extub <- survfit(Surv(time_extub, event_extub) ~ cluster, data = dataset_surv)

# Plot
extub_plot <- ggsurvplot(
  km_extub,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  linetype = "strata",
  ggtheme = theme_bw(),
  legend.title = "",
  xlim = c(0, 80),
  break.time.by = 5,
  xlab = "Jours depuis admission en ICU",
  ylab = "Probabilité de rester ventilé (non extubé)"
)

extub_plot
# ---- Modèle de Cox pour obtenir HR ----
cox_model <- coxph(Surv(time_extub, event_extub) ~ cluster, data = dataset_surv)

summary(cox_model)

# Récupérer HR et IC 95%
hr <- exp(coef(cox_model))
ci <- exp(confint(cox_model))

hr
ci

pdf("extub_plot_clusters.pdf", width = 7.5, height = 6)
extub_plot
dev.off()



# Supplementary FIGURE 10

library(MCPcounter)

expr <- read.table("normalized_human_gene_expression.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t")

# load signature manually
sig <- read.table("genes.txt", header = TRUE, sep = "\t")
probs <- read.table("probesets.txt", header = TRUE, sep = "\t")

# run MCPcounter using internal function workaround
res <- MCPcounter.estimate(
  as.matrix(expr),
  featuresType = "HUGO_symbols"
)


res


write.table(res, file = "MCPcounter_rpkm.txt", 
            sep = "\t", quote = FALSE, row.names = T)



# Supplementary FIGURE 11


### Donut plot 


#########################################BUBBLE PIE CHART donut automatic sample

library(dplyr)
library(tidyr)
library(ggplot2)
library(PieGlyph)
library(ggiraph)
library(reshape2)
library(ggplot2)
library(ggforce)

df <- read.delim("bacteriome.txt", sep = "\t", check.names = T, row.names = 1)

head(df)

meta <- read.delim("meta.txt", sep = "\t", check.names = F, row.names = 1)


df <-as.matrix(df)
df <- reshape2::melt(df)

colnames(df) <- c("bacteria", "SAMPLEID", "value")

head(df)

df$DAY = sapply(as.vector(df$SAMPLEID), function(x) meta[x,"DAY"])
df$cluster = sapply(as.vector(df$SAMPLEID), function(x) meta[x,"cluster"])
df$SAMPLEID = sapply(as.vector(df$SAMPLEID), function(x) meta[x,"SAMPLEID"])
df$Patient = sapply(as.vector(df$SAMPLEID), function(x) meta[x,"Patient"])

df <- df %>%
  group_by(SAMPLEID, DAY, cluster) %>%
  tidyr::complete(bacteria)

pdf("test2.pdf", height = 9, width = 10)

ggplot(df, aes(
  x = DAY,
  y = reorder(Patient, DAY)
)) +
  
  geom_pie_glyph(
    slices = "bacteria",   # ← OUTSIDE aes()
    values = "value",      # ← OUTSIDE aes()
    radius = 0.5,
    colour = "black"
  )  +
  scale_x_continuous(breaks = seq(-7, max(df$DAY), by = 1))+
  labs(
    title = "",
    x = "HAP onset day",
    y = "Patients"
  ) +
  
  theme_minimal() +
  
  scale_fill_manual(values = c(
    "red", "orange", "green3", "blue2", "purple", "grey"
  )) +
  
  theme(
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 16),
    legend.text = element_text(size = 6),
    legend.position = "none",
    legend.direction = "vertical"
  ) +
  
  facet_grid(~cluster, scales="free")


dev.off()


library(dplyr)
library(ggplot2)
library(reshape2)

# Load your abundance matrix
matrix <- read.delim("bacteriome.txt", header = TRUE, sep = "\t", check.names = F)

# Load metadata (Sample -> Cluster)
metadata <- read.delim("meta_new.txt", header = TRUE, sep = "\t", check.names = F)

# Melt the abundance matrix
melted_data <- reshape2::melt(matrix, id.vars = "fam")
colnames(melted_data) <- c("Genus", "SAMPLEID", "Value")

# Merge with metadata to get cluster info
melted_data <- melted_data %>%
  left_join(metadata[, c("SAMPLEID", "cluster")], by = "SAMPLEID")

# Calculate median relative abundance per Genus per cluster
median_cluster <- melted_data %>%
  group_by(cluster, Genus, SAMPLEID) %>%
  summarise(MeanAbundance  = median(Value, na.rm = TRUE), .groups = "drop")

# Normalize per cluster to get percentages
median_cluster <- median_cluster %>%
  group_by(cluster) %>%
  mutate(RelAbundance = MeanAbundance  / sum(MeanAbundance ) * 100) %>%
  ungroup()

# Plot bar plot RPKM

# Define genus colors
genus_colors <- c(
  "Prevotella"    = "green3",
  "Streptococcus" = "blue",
  "Haemophilus"   = "orange2",
  "Veillonella"   = "purple",
  "Fusobacterium" = "red"
)

median_line <- melted_data %>%
  group_by(cluster) %>%
  summarise(med = mean(Value, na.rm = TRUE))%>%
  mutate(label = paste0("mean = ", round(med, 0)))

g <- ggplot(melted_data, aes(x = SAMPLEID, y = Value, fill = Genus)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 1) +
  
  scale_fill_manual(values = genus_colors) +
  geom_hline(
    data = median_line,
    aes(yintercept = med),
    linetype = "dashed",
    linewidth = 1,
    color = "black"
  ) +
  geom_text(
    data = median_line,
    aes(x = 1, y = med, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -0.5,
    size = 5, fontface = "bold"
  ) +
  
  labs(
    title = "",
    x = "",
    y = "RPKM"
  ) +
  
  theme_classic() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.position = "none",
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(face = "bold", size = 8, angle = 90),
    axis.text.y = element_text(face = "bold", size = 14)
  ) +
  facet_grid(~cluster, scales = "free") + coord_cartesian(ylim = c(0, 100000))

g


pdf("RPKM_core_sample_clusters.pdf",width=10,height=5)
print(g)
dev.off()



# Supplementary FIGURE 12


#################### library size bacterial

df <- read.delim("meta.txt", stringsAsFactors = FALSE)

# Ensure variables are correctly formatted
df$pneumotype <- as.factor(df$pneumotype)

# Optional but recommended: log-transform skewed variables
df$log_expression <- log10(df$gene_expression + 1)
df$log_depth <- log10(df$sequencing_depth+1)
df$log_load <- log10(df$bacterial_load + 1)

# Fit linear model
model <- lm(gene_expression ~ pneumotype + bacterial_load + sequencing_depth, data = df)

# Summary of model
summary(model)

coef_table <- summary(model)$coefficients

pneumotype_rows <- grep("pneumotype", rownames(coef_table))

p_values <- coef_table[pneumotype_rows, "Pr(>|t|)"]
p_values


library(ggplot2)

coef_table <- summary(model)$coefficients

# keep only predictors (remove intercept)
plot_df <- data.frame(
  term = rownames(coef_table),
  estimate = coef_table[, "Estimate"],
  se = coef_table[, "Std. Error"],
  p = coef_table[, "Pr(>|t|)"]
)

plot_df <- plot_df[plot_df$term != "(Intercept)", ]

# create labels with significance
plot_df$label <- paste0(
  round(plot_df$estimate, 2),
  "\n(p=", round(plot_df$p, 3), ")"
)

library(ggplot2)

# add significance flag
plot_df$significant <- plot_df$p < 0.05

ggplot(plot_df, aes(x = reorder(term, estimate), y = estimate)) +
  
  # points colored by significance
  geom_point(aes(color = significant), size = 4) +
  
  # confidence intervals
  geom_errorbar(aes(
    ymin = estimate - 1.96 * se,
    ymax = estimate + 1.96 * se,
    color = significant
  ), width = 0.2, alpha = 0.8) +
  
  # zero reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  # optional labels (cleaner: comment out if too busy)
  geom_text(aes(label = paste0("p=", round(p, 3))),
            hjust = 1.5, size = 5) +
  
  coord_flip() +
  
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "#0072B2"),
    labels = c("Not significant", "Significant"),
    name = ""
  ) +
  
  labs(
    x = "",
    y = "Effect size (β ± 95% CI)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )


# Supplementary FIGURE 13


########
library(ggplot2)
library(dplyr)
library(broom)
library(purrr)

# Load data
expr <- read.delim("functions.txt", header = TRUE, row.names = 1, check.names = FALSE)
group <- read.delim("meta.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Match samples
common_samples <- intersect(colnames(expr), rownames(group))
expr <- expr[, common_samples]
group <- group[common_samples, , drop = FALSE]

# Binary outcome 
y <- ifelse(group$cluster == "P2_high_risk", 1, 0)
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

pdf("glm_avirvirnew_func.pdf",width=13,height=5)

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
  scale_fill_manual(values = c("P2_high_risk" = "indianred1", "P1_low_risk" = "limegreen")) +
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



### plot bacteria cog with covariates

# -----------------------------
# Load bacteria-function table
# -----------------------------
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

bac <- read.delim(
  "bacteriome.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# -----------------------------
# Function name cleaning
# -----------------------------
clean_fun <- function(x){
  
  x <- trimws(x)
  
  x <- gsub("^COG[: ]*", "", x)
  
  x <- gsub("\\s+", " ", x)
  
  x <- gsub("_", " ", x)
  
  x
}

bac$Function <- clean_fun(bac$Function)
glm_sig$Function <- clean_fun(glm_sig$Function)

# -----------------------------
# Top bacteria per function
# -----------------------------
top_bac <- bac %>%
  group_by(Function, bacteria) %>%
  summarise(
    mean = sum(mean),
    .groups = "drop"
  ) %>%
  group_by(Function) %>%
  slice_max(
    order_by = mean,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

# -----------------------------
# Function-level effect
# -----------------------------
glm_func <- glm_sig %>%
  group_by(Function) %>%
  summarise(
    effect = mean(estimate),
    pval = min(pval),
    .groups = "drop"
  )

# -----------------------------
# Merge
# -----------------------------
plot_df <- inner_join(
  top_bac,
  glm_func,
  by = "Function"
)

# Check overlap
cat(
  "Functions in GLM:",
  nrow(glm_func),
  "\nFunctions plotted:",
  length(unique(plot_df$Function)),
  "\n"
)

# -----------------------------
# Order functions by effect size
# -----------------------------
func_order <- glm_func %>%
  arrange(effect) %>%
  pull(Function)

plot_df$Function <- factor(
  plot_df$Function,
  levels = func_order
)

# -----------------------------
# Format taxa names
# -----------------------------
plot_df$bacteria <- gsub(
  "_",
  "\n",
  plot_df$bacteria
)

# -----------------------------
# Bubble plot
# -----------------------------
pdf("glm_COG_taxa_species.pdf",width = 14,height = 12)

ggplot(
  plot_df,
  aes(
    y = Function,
    x = bacteria
  )
) +
  geom_point(
    aes(
      size = log10(mean + 1),
      fill = effect
    ),
    shape = 21,
    color = "black",
    stroke = 0.4
  ) +
  coord_flip() +
  scale_size_continuous(
    name = "Mean abundance",
    range = c(3,10)
  ) +
  scale_fill_gradient2(
    low = "limegreen",
    mid = "white",
    high = "indianred1",
    midpoint = 0,
    name = "GLM effect\n(log-odds)"
  ) +
  scale_x_discrete(
    labels = function(x)
      str_wrap(x, width = 35)
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme(
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1, color="black"
    ),
    axis.text.y = element_text(
      size = 11,
      face = "bold.italic" ,color="black"
    ),
    plot.title = element_text(
      size = 16,
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

dev.off()


# Supplementary FIGURE 14


library("phyloseq")
library(ANCOMBC)
library(ggpubr)
library(reshape2)
library(vegan)
library(dplyr)

otu = read.delim(file="data/metabolic_pathways.txt",header=T,row.names=1,check.names=F)
meta = read.delim(file="meta.txt",header=T,row.names = 1, check.names=F) 
#tax = read.delim(file="tax.txt",header=T,row.names = 1, check.names=F)

tax <- data.frame(
  Phylum = rep("_", nrow(otu)),
  Class = rep("_", nrow(otu)),
  Order = rep("_", nrow(otu)),
  Family = rep("_", nrow(otu)),
  Genus = rep("_", nrow(otu)),
  Species = rep("_", nrow(otu))
)
row.names(tax) <- row.names(otu)
#otu<-otu[apply(otu,1,function(x)length(x[x>0]))>=0.5*ncol(otu),]

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))
physeq = phyloseq(OTU, TAX)
physeq = merge_phyloseq(physeq, sample_data(meta))

# Normalize by relative abundance
physeq <- transform_sample_counts(physeq, function(x) x / sum(x))

out  = ancombc(data = physeq,
               formula = "GROUP + IFN_dose + DAY + treatment_antibio + treatment_cortico + icuo + icua + Patient",
               p_adj_method = "fdr", lib_cut = 0,
               group = "GROUP", neg_lb = TRUE, tol = 1e-5,
               max_iter = 100, conserve = TRUE, alpha = 0.2, global = TRUE)

res = out$res


results_df <- data.frame(
  Feature = res[["lfc"]][["taxon"]],
  pvalue = res$p_val,
  qvalue = res$q_val,
  LFC = res$lfc
)

results_sig = results_df[results_df$pvalue.GROUPC2< 0.05,]

otu_sig = otu[rownames(otu) %in% results_sig$Feature, ]


library(ggplot2)

otu_presence <- otu > 0

# Step 4: Summarize the presence of each OTU across samples as a percentage
otu_presence_count <- rowSums(otu_presence)  # Counts of presence across samples
total_samples <- ncol(otu)  # Total number of samples
otu_percentage <- (otu_presence_count / total_samples) * 100  # Convert to percentage

# Step 5: Filter OTUs based on minimum sample presence (optional)
min_samples <- 1  # Set your threshold here
otu_matrix_filtered <- otu[otu_percentage >= (min_samples / total_samples) * 100, ]
otu_percentage_filtered <- otu_percentage[otu_percentage >= (min_samples / total_samples) * 100]

# Step 6: Create a data frame for plotting, merging with Family annotation
otu_summary <- data.frame(
  OTU = rownames(otu_matrix_filtered),
  Occurrence_Percentage = otu_percentage_filtered
)


results_sig$Occurrence_Percentage = sapply(as.vector(results_sig$Feature), function(x) otu_summary[x,"Occurrence_Percentage"])

write.table(results_sig, file = "results_sig_metab_clusters.txt", sep = "\t", col.names = NA)



ggplot(
  results_sig,
  aes(
    x = LFC.GROUPC2,
    y = factor(
      Feature,
      levels = results_sig %>%
        arrange(LFC.GROUPC2) %>%
        pull(Feature)
    )
  )
) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_point(
    aes(
      color = LFC.GROUPC2 > 0,
      size = Occurrence_Percentage
    )
  ) +
  scale_color_manual(
    values = c("TRUE" = "limegreen", "FALSE" = "indianred1"),
    guide = "none"
  ) +
  labs(
    x = "",
    y = "Predicted bacterial metabolic pathways",
    size = "Occurrence Percentage"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(
      size = 14,
      face = "bold.italic",
      color = ifelse(results_sig$Feature == "", "red", "black")
    )
  )



# Supplementary FIGURE 15 - see Figure4A.R

# Supplementary FIGURE 16-17


library(UpSetR)
library(RColorBrewer)
library(dplyr)

#C1 = P2_high_risk
#C2 = P1_low_risk

# === Read files ===
viral    <- read.delim("data/filtered_validation_viral_all.txt", stringsAsFactors = FALSE)
bacterial <- read.delim("data/filtered_validation_bacterial_all.txt", stringsAsFactors = FALSE)
combined  <- read.delim("data/filtered_validation_MIX_all.txt", stringsAsFactors = FALSE)

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


# Supplementary FIGURE 18


#########PLSDA
# SPLS

library(mixOmics)


Clinical <- read.delim("metadata", header = TRUE, sep = "\t", row.names = 1)
Y <- factor(Clinical$cluster)
Y_numeric <- as.numeric(factor(Y))
summary(Y)

X1 <- read.delim("data/RPKM_human.txt", header = TRUE, sep = "\t", row.names = 1, check.names = F)
X1<-t(X1)


########################

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

# Supplementary FIGURE 19 - see Figure2D.R

# Supplementary FIGURE 20

############################
# VIRUS–BACTERIA INTERACTIONS
# Dotplot using Pearson correlation coefficients
############################

######## LOAD LIBRARIES ########
library(dplyr)
library(purrr)
library(ggplot2)

######## LOAD DATA ########
data <- read.delim(
  "VH.txt",
  check.names = FALSE
)

######## DEFINE BACTERIA (GENERA) ########
bac_vars <- names(data)[
  !grepl("^a[0-9]*_", names(data)) &
    !(names(data) %in% c("SAMPLEID", "DAY", "cluster"))
]

######## KEEP VALID GENERA ########
valid_genus <- bac_vars[
  sapply(bac_vars, function(g) {
    any(grepl(
      paste0("^a[0-9]*_", g, "$"),
      names(data)
    ))
  })
]

######## CLUSTERS ########
clusters <- unique(data$cluster)

######## ANALYSIS ########
results <- map_dfr(clusters, function(cl) {
  
  data_cl <- data %>%
    filter(cluster == cl)
  
  map_dfr(valid_genus, function(genus) {
    
    ## find all viral copies for this genus
    virus_vars <- grep(
      paste0("^a[0-9]*_", genus, "$"),
      colnames(data_cl),
      value = TRUE
    )
    
    if(length(virus_vars) == 0) return(NULL)
    
    map_dfr(virus_vars, function(virus_var) {
      
      tmp <- data_cl %>%
        select(all_of(c(genus, virus_var, "DAY"))) %>%
        na.omit()
      
      if(nrow(tmp) < 5) return(NULL)
      
      if(sd(tmp[[genus]]) == 0 ||
         sd(tmp[[virus_var]]) == 0) return(NULL)
      
      ## Pearson correlation
      cor_test <- cor.test(
        tmp[[genus]],
        tmp[[virus_var]],
        method = "pearson"
      )
      
      data.frame(
        Cluster = cl,
        Genus = genus,
        Virus = virus_var,
        Correlation = unname(cor_test$estimate),
        P_value = cor_test$p.value,
        N = nrow(tmp)
      )
    })
  })
})

######## FDR CORRECTION ########
results <- results %>%
  group_by(Cluster) %>%
  mutate(
    FDR = p.adjust(
      P_value,
      method = "fdr"
    )
  ) %>%
  ungroup()

######## EXPORT ########
write.table(
  results,
  "association_hostvirus_correlations.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

######## PLOT DATA ########
plot_df <- results %>%
  mutate(
    Genus = reorder(
      Genus,
      Correlation
    ),
    
    Significance = case_when(
      P_value < 0.05 &
        Cluster == "P1_low_risk" ~ "P1_low_risk",
      
      P_value < 0.05 &
        Cluster == "P2_high_risk" ~ "P2_high_risk",
      
      TRUE ~ "Not significant"
    ),
    
    P_label = paste0(
      "P=",
      signif(P_value, 2)
    )
  )

######## PLOT ########
q <- ggplot(
  plot_df,
  aes(
    x = Correlation,
    y = reorder(Genus, Correlation),
    color = Significance
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  geom_point(
    aes(size = -log10(P_value)),
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    range = c(3, 10)
  ) +
  
  geom_text(
    aes(label = P_label),
    size = 3.5,
    vjust = -1.2,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "P1_low_risk" = "limegreen",
      "P2_high_risk" = "indianred1",
      "Not significant" = "grey60"
    )
  ) +
  
  facet_grid(
    ~ Cluster,
    scales = "free"
  ) +
  
  coord_cartesian(
    xlim = c(-1, 1)
  ) +
  
  theme_linedraw() +
  
  theme(
    axis.text = element_text(
      size = 20
    ),
    
    axis.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 20,
      face = "bold.italic"
    ),
    
    legend.title = element_blank(),
    
    strip.text = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.key.size = unit(
      0.8,
      "cm"
    ),
    
    legend.text = element_text(
      size = 14
    )
  ) +
  
  labs(
    x = "Correlation coefficient (r)",
    y = "Genus",
    size = "-log10(P-value)"
  )

q

######## SAVE ########
pdf(
  "virus_bacteria_correlation_dotplot.pdf",
  width = 12,
  height = 8
)

print(q)

dev.off()



################### ecological dynamics analysis bacteria per virus per cluster

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# =========================================
# READ
# =========================================
data <- read.delim("VH.txt",
                   check.names = FALSE)

# =========================================
# LOG TRANSFORM
# =========================================
#data <- data %>%
mutate(across(
  -c(SAMPLEID, DAY, cluster),
  ~ log10(.x + 1)
))

# =========================================
# KEEP ONLY TARGET CLUSTERS
# =========================================
data <- data %>%
  filter(cluster %in% c("P1_low_risk",
                        "P2_high_risk"))

# =========================================
# COLORS
# =========================================
custom_colors <- c(
  "P2_high_risk" = "indianred1",
  "P1_low_risk" = "limegreen"
)

# =========================================
# BACTERIA VARIABLES
# =========================================
bac_vars <- names(data)[
  !grepl("^a_", names(data)) &
    !(names(data) %in%
        c("SAMPLEID", "DAY", "cluster"))
]

glm_list <- list()

for(genus in bac_vars){
  
  virus_var <- paste0("a_", genus)
  
  if(!(virus_var %in% names(data))){
    next
  }
  
  # build model dataset
  tmp <- data %>%
    select(all_of(c(genus, virus_var, "cluster", "DAY"))) %>%
    na.omit()
  
  # skip low info cases
  if(nrow(tmp) < 4){
    next
  }
  
  if(sd(tmp[[genus]]) == 0 | sd(tmp[[virus_var]]) == 0){
    next
  }
  
  # GLM (use gaussian unless you specify otherwise)
  model <- glm(
    formula = as.formula(
      paste0(
        virus_var,
        " ~ ",
        genus,
        " + cluster + DAY"
      )
    ),
    data = tmp,
    family = gaussian()
  )
  
  sm <- summary(model)$coefficients
  
  # extract effect of genus (main association of interest)
  if(genus %in% rownames(sm)){
    
    glm_list[[length(glm_list) + 1]] <- data.frame(
      Genus   = genus,
      Estimate = sm[genus, "Estimate"],
      Pvalue   = sm[genus, "Pr(>|t|)"]
    )
  }
}

corr_df <- bind_rows(glm_list)

# FDR correction
corr_df$FDR <- p.adjust(corr_df$Pvalue, method = "fdr")

# =========================================
# KEEP GENERA SIGNIFICANT IN
# AT LEAST ONE CLUSTER
# =========================================
sig_genus <- corr_df %>%
  filter(Pvalue < 1) %>%
  pull(Genus) %>%
  unique()

# =========================================
# LONG FORMAT BACTERIA
# =========================================
bac_long <- data %>%
  select(SAMPLEID,
         DAY,
         cluster,
         all_of(sig_genus)) %>%
  pivot_longer(
    cols = all_of(sig_genus),
    names_to = "Genus",
    values_to = "Abundance"
  ) %>%
  mutate(Type = "Bacteria")

# =========================================
# LONG FORMAT VIRUSES
# =========================================
vir_long <- data %>%
  select(
    SAMPLEID,
    DAY,
    cluster,
    all_of(paste0("a_", sig_genus))
  ) %>%
  
  pivot_longer(
    cols = starts_with("a_"),
    names_to = "Genus",
    values_to = "Abundance"
  ) %>%
  
  mutate(
    Genus = str_remove(Genus, "^a_"),
    Type = "Virus"
  )

# =========================================
# MERGE
# =========================================
plot_data <- bind_rows(
  bac_long,
  vir_long
)

# =========================================
# ADD CORRELATION LABELS
# =========================================
label_df <- corr_df %>%
  filter(Genus %in% sig_genus) %>%
  mutate(
    label = paste0(
      "Estimate=",
      round(Estimate, 2),
      "\nPvalue=", signif(Pvalue, 2)
    )
  )
# =========================================
# PLOT
# =========================================
p <- ggplot(
  plot_data,
  aes(
    x = DAY,
    y = Abundance,
    color = cluster,
    fill = cluster,
    linetype = Type
  )
) +
  
  stat_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    linewidth = 2,
    alpha = 0.15
  ) +
  
  facet_wrap(
    ~ Genus,
    scales = "free"
  ) +
  
  scale_color_manual(
    values = custom_colors
  ) +
  
  scale_fill_manual(
    values = custom_colors
  ) +
  
  scale_x_continuous(
    breaks = seq(
      min(plot_data$DAY),
      max(plot_data$DAY),
      by = 1
    )
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(
      size = 20,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    strip.text = element_text(
      size = 20,
      face = "bold"
    ),
    
    legend.position = "bottom",
    legend.text = element_text(size=20)
  ) +
  
  xlab("Days related to HAP onset") +
  ylab("log10 abundance") +
  
  guides(fill = "none")


p

# Supplementary FIGURE 21

########### tkna importance 

# ==============================
# Importance score calculation
# ==============================

# Load libraries
library(dplyr)

# Read input table (replace with your file name)
df <- read.delim("data/probabilities_to_randomly_find_nodes.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

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


# Supplementary FIGURE 22 - see Figure4d.R

# Supplementary FIGURE 23 - see Figure6h-i.R


