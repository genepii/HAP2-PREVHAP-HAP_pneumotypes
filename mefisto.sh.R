Sys.setenv(RETICULATE_PYTHON = "C:/Users/debajyoti/anaconda3/envs/py38/python.exe")
reticulate::use_condaenv("C:/Users/debajyoti/anaconda3/envs/py38/",
                         conda =  "C:/Users/debajyoti/anaconda3/Scripts/conda.exe", 
                         required = TRUE)


library(reshape2)
library(MOFA2)


dt = read.delim("dt.txt")
meta = read.delim("meta.txt", row.names = 1)
meta$sample = rownames(meta)
head(meta)


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
mefisto_opts$new_values <- matrix(-11:-1, nrow =1) # set time points to interpolate factors to

obj <- prepare_mofa(obj, model_options = model_opts,
                   mefisto_options = mefisto_opts,
                   training_options = train_opts,
                   data_options = data_opts)

obj <- run_mofa(obj, use_basilisk = F,
               outfile = file.path(getwd(), "model.hdf5"))

# saveRDS(obj , "obj.rds")

plot_variance_explained(obj)
plot_variance_explained(obj, plot_total = T)[[2]]

######################

###MOFA2: training a model in R
library(tidyr)
library(tidyverse)
library(reshape2)
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
library(MOFA2)
library(basilisk)
library(dplyr)
library(compositions)
library(dplyr)
library(cowplot)


####################### reformat

df1 <- read.delim("bacrep.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df2 <- read.delim("Avir.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df3 <- read.delim("vir.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df4 <- read.delim("bacrepf.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df5 <- read.delim("Avirf.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df6 <- read.delim("virf.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df7 <- read.delim("metab.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df8 <- read.delim("hum.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
df9 <- read.delim("humd.txt", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

bacrep <- df1 %>% mutate(view = 'bacrep')
Avir <- df2 %>% mutate(view = 'vir')
vir <- df3 %>% mutate(view = 'vir')
bacrepf <- df4 %>% mutate(view = 'bacrep')
Avirf <- df5 %>% mutate(view = 'vir')
virf <- df6 %>% mutate(view = 'vir')
metab <- df7 %>% mutate(view = 'bacrep')


dt <- bind_rows(bacrep, Avir,vir,bacrepf,Avirf,virf,metab)
head(dt)


# Rename the column
colnames(dt)[colnames(dt) == "Family"] <- "sample"
colnames(dt)[colnames(dt) == "Sample"] <- "feature"
colnames(dt)[colnames(dt) == "Value"] <- "value"
head(dt)

# Swap the 'feature' and 'sample' columns
dt <- dt[, c("sample", "feature", "value", "view")]

head(dt)

dt <- dt %>%
  mutate(value = log10(value+1))

dt[dt=="-Inf"] <- NA

head(dt)

a <- ggdensity(dt, x="value", fill="view") +
  facet_wrap(~view, nrow=1, scales = "free")
a

dt$group = sapply(as.vector(dt$sample), function(x) meta[x,"group"])

head(dt)


write.table(dt, file = "dt.txt", sep = "\t", row.names = F)

setDT(dt)

unique(dt$view)

length(unique(dt$Sample))

dt[,length(unique(otu)),by="view"]


##################


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
plot_variance_explained(obj)

plot_factor_cor(obj)
get_scales(obj)

variance_data <- get_variance_explained(obj)

str(variance_data)
variance_factors <- variance_data$r2_per_factor

#write.table(variance_factors, file = "variance_explained.txt", sep = "\t")


plot_variance_explained(obj)

#
var_explained_values <- obj@cache[["variance_explained"]]$r2_per_factor
rowSums(sapply(var_explained_values, function(e) rowSums(e, na.rm = TRUE)))
#

r2 <- obj@cache$variance_explained$r2_per_factor[[1]]


r2.dt <- r2 %>%
  as.data.table %>% .[,factor:=as.factor(1:obj@dimensions$K)] %>%
  melt(id.vars=c("factor"), variable.name="view", value.name = "r2") %>%
  .[,cum_r2:=cumsum(r2), by="view"]

ggline(r2.dt, x="factor", y="cum_r2", color="view") +
  labs(x="Factor number", y="Cumulative variance explained (%)") +
  theme(
    legend.title = element_blank(), 
    legend.position = "top",
    axis.text = element_text(size=rel(0.8))
  )


##
library(tidyr)
library(dplyr)

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


############### heatmap REORDER BY CLUSTER 

library(pheatmap)

matrix <- read.delim("HM.txt", header = TRUE, sep = "\t", row.names = 1, check.names = F)
matrix<-matrix[apply(matrix,1,function(x)length(x[x>0]))>=1*ncol(matrix),]
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
  bacteria = c("healthy" = "limegreen", "others" = "purple"),
  cluster = c("C1" = "indianred1", "C2" = "limegreen"),
  Function = function_colors
)


heatmap <- pheatmap(
  log10(matrix),
  color=colorRampPalette(c("red","white","green4"))(100),
  annotation_row = my_sample_row,
  #annotation_col = my_sample_col,
  cluster_cols = F,
  #cellheight = 50,
  #cellwidth = 20,
  cluster_rows =F,
  annotation_colors = ann_colors,
  #fontsize = 16, 
  legend = T,
  #clustering_distance_cols = "correlation",
  annotation_legend = T,
  show_colnames = T, show_rownames = F
)

heatmap

pdf("heatmap.pdf", width = 10, height = 7)
heatmap
dev.off()



######################################### Partitionning clustering
library(stats)
library(cluster)
library(factoextra)

matrix <- read.delim("all_formatted_factors_all.txt", header = TRUE, sep = "\t", row.names = 1, check.names = F)
meta <- read.delim("meta.txt", stringsAsFactors = FALSE)
head(matrix)


fviz_nbclust(t(matrix), kmeans, method = "wss") +
  geom_vline(xintercept = "2", linetype = 3)+
  labs(subtitle = "Elbow method")

fviz_nbclust(t(matrix), kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette Method")

#fviz_nbclust(t(matrix), kmeans, method = "gap_stat") +
  labs(subtitle = "Gap_stat Method")

gap_stat <- clusGap(t(matrix), FUN = kmeans, nstart = 25, K.max = 10, B = 500)
fviz_gap_stat(gap_stat)


#best_k <- maxSE(gap_stat$Tab[,"gap"], gap_stat$Tab[,"SE.sim"], method = "firstSEmax")
#cat("Best number of clusters according to gap statistic:", best_k, "\n")
#set.seed(123)  # reproducibility
#clust <- kmeans(t(matrix), centers = best_k, nstart = 25)
#cluster_df <- data.frame(cluster = clust$cluster)
#write.table(cluster_df, "cluster_assignments.txt", row.names = TRUE, col.names = NA, sep = "\t")

optimal_k <- $nb

# Perform k-means clustering
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(t(matrix), centers = optimal_k, nstart = 25)
# View the results
print(kmeans_result)

dim(t(matrix))

#pdf("plot_mofa_clusters.pdf", width = 7.5, height = 5)
famd_clust_a <- fviz_cluster(kmeans_result,data = t(matrix),
                             show.clust.cent = F,
                             ggtheme = theme_classic(),
                             main = "", geom = "point", palette = c("blue", "skyblue3", "green","yellow4","blue","gold","grey"))
famd_clust_a
#dev.off()

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
  SAMPLEID = meta[, 1],              
  cluster  = kmeans_result$cluster    
)
#write.table(kmeans_result_df, "cluster_assignments.txt", row.names = TRUE, col.names = NA, sep = "\t")

# Perform the join explicitly on sampleid
joined_data <- dplyr::left_join(meta, kmeans_result_df, by = "SAMPLEID")

datat<- joined_data

################# reformat meta_patients for KP analysis

library(dplyr)
library(stringr)

# Your data frame
datat = read.delim(file="meta.txt",header=T,check.names=F)

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

write.table(df_max_day, file = "filtered_data.txt", sep = "\t", quote = FALSE,col.names = NA)

datat<- df_max_day

library(survival)
library(survminer)

#write.table(kmeans_result_df, file = "kmeans.txt", sep = '\t')
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




