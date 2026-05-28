############################################################################################################################
#### CODE USED TO GENERATE FIGURE 6
############################################################################################################################

# A
#Prepare TKNA input

# TKNA Input Preparation

## Overview

This workflow describes the preparation of input files used for the TKNA analysis integrating host transcriptomic signatures, bacterial abundance profiles, and viral abundance profiles.

The final dataset combines multi-omic features into a unified matrix suitable for network construction and downstream TKNA analysis.

---

# Input Files

## 1. Human Signature Matrix

File:

```r
inputlist.txt
```

Contains normalized host molecular signatures (e.g. genes, pathways, or transcriptomic features) across all samples.

Rows:

* Features

Columns:

* Sample IDs

---

## 2. Bacterial Differential Abundance Matrix

File:

```r
ancomb_bac.txt
```

Contains bacterial abundance or ANCOMBC-derived bacterial features.

Rows:

* Bacterial taxa/features

Columns:

* Sample IDs

---

## 3. Viral Differential Abundance Matrix

File:

```r
ancombc_vir.txt
```

Contains viral abundance or ANCOMBC-derived viral features.

Rows:

* Viral taxa/features

Columns:

* Sample IDs

---

# Data Integration

The three matrices were merged into a single feature matrix using:

```r
SampleID
```

as the common identifier across datasets.

The merge operation was performed by matching shared sample IDs between:

* host signatures
* bacterial features
* viral features

Only matched samples were retained in the final matrix.

---

# Data Transformation

After merging all datasets, a log-transformation was applied:

```r
log10(x + 1)
```

to reduce skewness and stabilize variance across omic layers.

The transformed matrix was exported as:

```r
data.csv
```

---

# TKNA Input Files

## data/data.csv

Main abundance/expression matrix used as TKNA input.

Columns:

* Samples

Rows:

* Features (genes, bacterial taxa, viral taxa)

Values:

* log10(x + 1) transformed abundances

---

## data/data_graph_map.csv

Metadata file containing sample annotations.

Includes:

* Sample IDs
* Cluster assignments/groups

Used for:

* phenotype mapping
* network coloring
* group comparisons

---

## data/type_map.csv

Feature annotation table describing feature categories.

Contains mappings between features and omic types:

* Gene
* Virus
* Bacteria

Used for:

* node annotation
* network visualization
* modality-specific analyses

---

## data/config.json

Configuration file containing TKNA analysis parameters.

Includes:

* filtering thresholds
* network construction parameters
* correlation settings
* statistical thresholds
* output settings

This file controls the behavior of the TKNA pipeline.

:

{
    "comparisonTreatments": [
        "P1",
        "P2"
    ],
    "differenceMethod": "mannwhitney",
    "differencePValueThresholds": {
        "individual": {
            "virus": 1,
            "bacteria": 1,
			"genes": 1,
            "fungi": 1
        },
        "combined": {
            "virus": 1,
            "bacteria": 1,
			"genes": 1,
            "fungi": 1
        },
        "corrected": {
            "virus": 1,
            "bacteria": 1,
			"genes": 1,
            "fungi": 1
        }
    },
    "foldChangeType": "mean",
    "foldChangeFilterMethod": "allsamesign",
    "networkTreatment": "P2",
    "correlationMethod": "spearman",
    "correlationPValueThresholds": {
        "individual": 0.05,
        "combined": 0.05,
        "corrected": 0.05
    },
    "correlationFilterMethod": "allsamesign"
}

---

## data/metadata.json

{ 
"name": "study_name", 
"experiments": [ 
{ 
"name": "example1", 
"dataFile": "data.csv", 
"treatmentMapFile": "data_graph_map.csv" 
}, 
{ 
"name": "example2", 
"dataFile": "data.csv", 
"treatmentMapFile": "data_graph_map.csv" 
} 
], 
"measurableTypeMapFile": "type_map.csv" 
} 

# Workflow Summary

```r
inputlist.txt
        +
ancomb_bac.txt
        +
ancombc_vir.txt
        ↓
Merge by SampleID
        ↓
log10(x + 1)
        ↓
data.csv
```

Additional annotation files:

* `data_graph_map.csv`
* `type_map.csv`
* `config.json`

are used for TKNA metadata integration and parameter specification.

---

# Notes

* Samples were matched using identical Sample IDs across all datasets.
* Missing samples were excluded during merging.
* Data transformation was applied uniformly across all omic layers.
* Feature annotations were harmonized prior to network analysis.


#Run TKNA python script and plos using Cytoscape
input_folder must contain data.csv, type_map.csv, data_graph_map.csv, metadata.json and config.json
### TKNA script
python ../reconstruction/intake_data.py --data-dir input_folder/ --out-file output/all_data_and_metadata.zip
python ../reconstruction/run.py --data-source output/all_data_and_metadata.zip --config-file config.json --out-file output/network_output.zip
python ../reconstruction/to_csv.py --data-file output/network_output.zip --config-file config.json --out-dir output/network_output
python ../analysis/assess_network.py --file output/network_output/correlations_bw_signif_measurables.csv --out-dir output/network_output/
python ../analysis/infomap_assignment.py --network output/network_output/network_output_comp.csv --network-format csv --map type_map.csv --out-dir output/network_output/
python ../analysis/louvain_partition.py --network output/network_output/network_output_comp.csv --network-format csv --map type_map.csv --out-dir output/network_output/
python ../analysis/find_all_shortest_paths_bw_subnets.py --network output/network_output/network_output_comp.csv --network-format csv --map type_map.csv --node-groups virus bacteria --out-dir output/network_output/
python ../analysis/calc_network_properties.py --network output/network_output/network_output_comp.csv --bibc --bibc-groups node_types --bibc-calc-type rbc --map type_map.csv --node-groups virus bacteria --out-dir output/network_output/
python ../random_networks/create_random_networks.py --template-network output/network_output/network.pickle --networks-file output/network_output/all_random_nws.zip 
python ../random_networks/compute_network_stats.py --networks-file output/network_output/all_random_nws.zip --bibc-groups node_types --bibc-calc-type rbc --stats-file output/network_output/random_network_analysis.zip --node-map type_map.csv --node-groups virus bacteria
python ../random_networks/synthesize_network_stats.py --network-stats-file output/network_output/random_network_analysis.zip --synthesized-stats-file output/network_output/random_networks_synthesized.csv
python ../visualization/dot_plots.py --pickle output/network_output/network.pickle --node-props output/network_output/node_properties.txt --network-file output/network_output/network_output_comp.csv --propx BiBC_virus_bacteria --propy Node_degrees --top-num 5 --top-num-per-type 5 --plot-dir output/network_output/plots/ --file-dir output/network_output/
python ../visualization/plot_density.py --rand-net output/network_output/random_networks_synthesized.csv --pickle output/network_output/inputs_for_downstream_plots.pickle --bibc-name BiBC_virus_bacteria
python ../visualization/plot_abundance.py --pickle output/network_output/inputs_for_downstream_plots.pickle --abund-data data.csv data.csv --metadata data_graph_map.csv data_graph_map.csv --x-axis Experiment --group-names NO_HAP HAP --group-colors blue pink



# B - G --> SEE boxplots in FIG3A.R

# K - M
# lOESS

###############loess plot per cluster###############loess plot per cluster###############loess plot per cluster
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

data <- read.csv("data/data.csv")

data <- data %>%
mutate(across(-c(SAMPLEID, DAY,  cluster), ~ log10(.x+1)))

custom_colors <- c("P2_high_risk" = "indianred1", "P1_low_risk" = "limegreen")

# 🔧 Variables
x_var <- "DAY"
group_var <- "cluster"

# Get all numeric variables except SAMPLEID and DAY
num_vars <- names(data)[sapply(data, is.numeric)]
num_vars <- setdiff(num_vars, x_var)

# Reshape into long format
data_long <- data %>%
  select(all_of(c(x_var, num_vars,group_var))) %>%
  pivot_longer(cols = all_of(num_vars),
               names_to = "Metric", values_to = "Value")

model <- lm(Value ~ DAY * cluster, data = data_long)
summary(model)

n_colors <- length(unique(data_long$Metric))

# Plot all in one LOESS curve
# Plot all in one LOESS curve with SE colored by group
loessPlot <- ggplot(data_long, aes_string(x = x_var, y = "Value", color = group_var, fill = group_var, group = group_var)) +
  stat_smooth(method = "loess", formula = y ~ x, se = T, linewidth = 2, alpha = 0.1) +  # alpha controls ribbon transparency
  scale_color_manual(values = custom_colors) +
  #geom_point()+
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(min(data_long$DAY), max(data_long$DAY), by = 1))  +
  #scale_y_continuous(limits = c(-1, 6))+
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    legend.title = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  xlab("Days related to HAP onset") +
  ylab("log10(TMM-normalized)") + facet_grid(~Metric, scales = "free")
loessPlot

pdf("cinétique_40genes_preonset.pdf",width=12,height=6)

loessPlot + coord_cartesian(xlim=c(-7,0))+

  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 1)+ 
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "HAP day",
    angle = 90,
    vjust = -0.5,
    hjust = 1.1,
    color = "grey40",
    size = 4
  )
  
dev.off() 




##############################################################