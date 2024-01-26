Construction of NicheNet’s ligand-target model
================
Robin Browaeys
2018-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/model_construction.Rmd", output_format = "github_document")
-->

This vignette shows how ligand-target prior regulatory potential scores
are inferred in the NicheNet framework. You can use the procedure shown
here to develop your own model with inclusion of context-specific
networks or removal of noisy irrelevant data sources. The networks at
the basis of NicheNet can be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7074291.svg)](https://doi.org/10.5281/zenodo.7074291).

# Background information about NicheNet’s prior ligand-target model

The prior model at the basis of NicheNet denotes how strongly existing
knowledge supports that a ligand may regulate the expression of a target
gene. To calculate this ligand-target regulatory potential, we
integrated biological knowledge about ligand-to-target signaling paths
as follows.

First, we collected multiple complementary data sources covering
ligand-receptor, signal transduction (e.g., protein-protein and
kinase-substrate interactions) and gene regulatory interactions (e.g.,
inferred from ChIP-seq and motifs). For information of all collected
data sources (link to the website of the database, etc), see [Data
source information](data_sources.xlsx)

Secondly, we integrated these individual data sources into two weighted
networks: 1) a ligand-signaling network, which contains protein-protein
interactions covering the signaling paths from ligands to downstream
transcriptional regulators; and 2) a gene regulatory network, which
contains gene regulatory interactions between transcriptional regulators
and target genes. To let informative data sources contribute more to the
final model, we weighted each data source during integration. These data
source weights were automatically determined via model-based parameter
optimization to improve the accuracy of ligand-target predictions (see
the vignette [Parameter optimization via
NSGA-II](parameter_optimization.md). In this vignette, we will show how
to construct models with unoptimized data source weigths as well.

Finally, we combined the ligand-signaling and gene regulatory network to
calculate a regulatory potential score between all pairs of ligands and
target genes. A ligand-target pair receives a high regulatory potential
if the regulators of the target gene are lying downstream of the
signaling network of the ligand. To calculate this, we used network
propagation methods on the integrated networks to propagate the signal
starting from a ligand, flowing through receptors, signaling proteins,
transcriptional regulators, and ultimately ending at target genes.

A graphical summary of this procedure is visualized here below:

![](images/workflow_model_construction.png)

# Construct a ligand-target model from all collected ligand-receptor, signaling and gene regulatory network data sources

Load the required packages and networks we will use to construct the
model.

``` r
library(nichenetr)
library(tidyverse)

# in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions

# The complete networks can be downloaded from Zenodo
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))
```

## Construct NicheNet’s ligand-target model from unoptimized data source weights

Construct the weighted integrated ligand-signaling and gene regulatory
network. In this first example, we give every data source the same
weight (as given by the `source_weights_df` data frame provided by
default by the nichenetr package). See the vignette showing how to use
NSGA-II to optimize data source weights and the hyperparameters if
interested in performing parameter optimization. For the hyperparameters
of the model (hub correction factors and damping factor), we will use
the optimized values (as given by the `hyperparameter_list` data frame
provided by default by the nichenetr package).

The ligand-signaling network hub correction factor and gene regulatory
network hub correction factor were defined as hyperparameter of the
model to mitigate the potential negative influence of over-dominant hubs
on the final model. The damping factor hyperparameter is the main
parameter of the Personalized PageRank algorithm, which we used as
network propagation algorithm to link ligands to downstream regulators.

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df = source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight)) 
```

Infer ligand-target regulatory potential scores based on the weighted
integrated networks

``` r
# in this example we will calculate target gene regulatory potential scores for TNF and the ligand combination TNF+IL6
ligands = list("TNF",c("TNF","IL6"))
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
```

Show some top target genes of the ligand TNF and the ligand combination
TNF+IL6

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##    NFKBIA      EDN1      MMP9      IRF1      IER3     NFKB1     PTGS2      PTX3     ICAM1      BMP2 
## 0.7119458 0.7118929 0.7107600 0.7105364 0.7091492 0.7091226 0.7082413 0.7082322 0.7055558 0.7038156
```

``` r
extract_top_n_targets("TNF-IL6",10,ligand_target_matrix)
##     ICAM1      IRF1      JUNB      IER3     CCND1       FOS     PTGS2     SOCS3     NFKB1      EDN1 
## 0.5947953 0.5312094 0.5247782 0.5202155 0.5002600 0.4729934 0.4636054 0.4634419 0.4603225 0.4578782
```

## Construct NicheNet’s ligand-target model from optimized data source weights

Now, we will demonstrate how you can make an alternative model with the
optimized data source weights (as given by the
`optimized_source_weights_df` data frame provided by default by the
nichenetr package)

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network,
                                                source_weights_df = optimized_source_weights_df %>% select(source, avg_weight) %>% rename(weight=avg_weight))

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight)) 

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##      PTX3      EDN1      IER3      BMP2     PTGS2   TNFAIP2      IRF1      IRF7   TNFAIP3     RIPK2 
## 0.2818886 0.2815736 0.2806793 0.2803148 0.2790848 0.2789345 0.2784311 0.2775592 0.2765948 0.2749011
```

# Change the data sources at the basis of the NicheNet ligand-target model

### Keep only specific data sources of interest

Now, we will demonstrate how you can decide which data sources to use in
the model you want to create. Let’s say for this example, that you are
interested in making a model that only consists of literature-derived
ligand-receptor interactions, signaling and gene regulatory interactions
from comprehensive databases and gene regulatory interactions inferred
from ChIP-seq. An annotation of the different data sources is given by
the `annotation_data_sources` data frame provided by default by the
nichenetr package)

``` r
annotation_data_sources$type_db %>% unique()
##  [1] "comprehensive_db" "literature"       "ptm"              "text_mining"      "directional_ppi"  "PPI"              "ChIP"            
##  [8] "motif"            "prediction"       "perturbation"
```

``` r
data_sources_to_keep = annotation_data_sources %>% filter(type_db %in% c("literature","comprehensive_db","ChIP")) %>% pull(source)

new_source_weights_df = source_weights_df %>% filter(source %in% data_sources_to_keep)
new_lr_network = lr_network %>% filter(source %in% data_sources_to_keep) 
new_sig_network = sig_network %>% filter(source %in% data_sources_to_keep)
new_gr_network = gr_network %>% filter(source %in% data_sources_to_keep)
```

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = new_lr_network, sig_network = new_sig_network, gr_network = new_gr_network, source_weights_df = new_source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##      SELE      MMP9     ICAM1     CASP3     VCAM1     CCND1    CDKN1A       MYC      TP53       JUN 
## 0.3120182 0.3116462 0.3103540 0.2993875 0.2987840 0.2097975 0.2067106 0.2059286 0.2032518 0.1960899
```

To give a second example: say that you don’t trust TF-target links
inferred from motif information and want to construct a model with all
data sources except the motif ones.

``` r
data_sources_to_remove = annotation_data_sources %>% filter(type_db %in% c("motif")) %>% pull(source)
data_sources_to_keep = annotation_data_sources$source %>% setdiff(data_sources_to_remove) 

new_source_weights_df = source_weights_df %>% filter(source %in% data_sources_to_keep)
new_lr_network = lr_network %>% filter(source %in% data_sources_to_keep) 
new_sig_network = sig_network %>% filter(source %in% data_sources_to_keep)
new_gr_network = gr_network %>% filter(source %in% data_sources_to_keep)
```

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = new_lr_network, sig_network = new_sig_network, gr_network = new_gr_network, source_weights_df = new_source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##      EDN1    NFKBIA      MMP9      PTX3     NFKB1     PTGS2      IER3      IRF1   TNFAIP2     ICAM1 
## 0.7096364 0.7080196 0.7080049 0.7077747 0.7073651 0.7070222 0.7059278 0.7055725 0.7037034 0.7032143
```

### Add own data sources to the NicheNet model

In addition to removing data sources, you can also add new data sources.
This could for example help you in making context-specific models, if
you would have a network or data containing context-specific
interactions of interest.

As input, we require a data source to contain directional interactions
between genes: these interactions are protein-protein or signaling
interactions for ligand-receptor and signaling data sources and a gene
regulatory interaction for gene regulatory data sources. The data
sources should be formatted in a data frame with following columns:
from, to and source. “from” denotes the source node “gene A” of the
directional interaction from gene A to B, “to” denotes the target node
“gene B” of this directional interaction, and “source” is a user-defined
name of this data source.

Here, we will show how you can download, process and integrate an online
data source within the NichenNet framework. As example, this is the data
source “Hub Proteins Protein-Protein Interactions” from the Harmonizome
portal
(<https://amp.pharm.mssm.edu/Harmonizome/dataset/Hub+Proteins+Protein-Protein+Interactions>).

``` r
input_file = "https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/hubs/gene_attribute_edges.txt.gz"
ppi_network = read_tsv(input_file, col_names = TRUE)

ppi_network = ppi_network %>% transmute(from=target,to=source) %>% 
    filter(from %in% geneinfo_human$symbol & to %in% geneinfo_human$symbol) # keep only interactions between genes with oficial gene symbols: optional step

# give your data source a name
ppi_network = ppi_network %>% mutate(source = "harmonizome_hub_ppi", database = "harmonizome") 

head(ppi_network)
## # A tibble: 6 × 4
##   from  to    source              database   
##   <chr> <chr> <chr>               <chr>      
## 1 LRIF1 GAD1  harmonizome_hub_ppi harmonizome
## 2 LRIF1 ID2   harmonizome_hub_ppi harmonizome
## 3 LRIF1 NOC2L harmonizome_hub_ppi harmonizome
## 4 LRIF1 ESR1  harmonizome_hub_ppi harmonizome
## 5 LRIF1 NR3C1 harmonizome_hub_ppi harmonizome
## 6 LRIF1 PPARG harmonizome_hub_ppi harmonizome
```

First, we will add this new data source to all other data sources.
Because this data sources contains intracellular protein-protein
interactions, we will consider this data source as a signaling data
source. As example, we will assign to this data source a weight of 1,
because we want it to have a strong contribution to the final model.

``` r
new_sig_network = sig_network %>% bind_rows(ppi_network)

new_network_weights_df = tibble(source = "harmonizome_hub_ppi", avg_weight = 1, median_weight = 1)
new_source_weights_df = optimized_source_weights_df %>% bind_rows(new_network_weights_df)
```

Now make this model

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = new_sig_network, gr_network = gr_network,
                                                source_weights_df = new_source_weights_df %>% select(source, avg_weight) %>% rename(weight=avg_weight))

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##      PTX3      EDN1      IER3      BMP2     PTGS2   TNFAIP2      IRF1      IRF7   TNFAIP3     RIPK2 
## 0.2818874 0.2815333 0.2806313 0.2803472 0.2790419 0.2789683 0.2783588 0.2775317 0.2765463 0.2748918
```

In some cases, it’s possible that you want that your data source will be
considered as only data source in a specific layer
(i.e. ligand-receptor, signaling or gene regulatory layer). Therefore,
we will show how you use this new data source as only data source part
of the signaling network.

``` r
new_sig_network = ppi_network

new_network_weights_df = tibble(source = "harmonizome_hub_ppi", avg_weight = 1, median_weight = 1)
new_source_weights_df = optimized_source_weights_df %>% bind_rows(new_network_weights_df)
```

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = new_sig_network, gr_network = gr_network,
                                                source_weights_df = new_source_weights_df %>% select(source, avg_weight) %>% rename(weight=avg_weight))

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = list("TNF")

ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
```

Show some top target genes of the ligand TNF

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##     CXCL8      PTX3      IRF7   TNFAIP2      EDN1      BMP2     NFKB1     RIPK2    NFKBIA      IER3 
## 0.2904951 0.2819874 0.2760319 0.2754550 0.2746321 0.2723789 0.2701484 0.2695589 0.2673814 0.2664578
```

## Final note

Most optimally, you would like to optimize the parameters again when
including own data sources. Instructions to do this are given in the
following vignette: [Parameter optimization via
NSGA-II](parameter_optimization.md):
`vignette("parameter_optimization", package="nichenetr")`

However, this optimization process takes a lot of time and requires the
availability of multiple cores to perform the optimization in parallel.
Because we demonstrate in the NicheNet paper that unoptimized models
also perform considerably well, data source weight optmization is not
necessary to have decent predictive ability.
