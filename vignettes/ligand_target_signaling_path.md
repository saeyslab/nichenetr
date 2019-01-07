### Infer signaling paths beween ligand(s) and target(s) of interest

To determine signaling paths between a ligand and target of interest, we look at which transcription factors are best regulating the target genes and are most closely downstream of the ligand. Then, the shortest paths between these transcription factors and the ligand of interests are determined and genes forming part in this path are considered as important signaling mediators. Finally, we give all interactions between the ligand, signaling mediators, transcription factors and target genes.

First, we will load the necessary networks to infer signaling paths between ligand and target genes of interest.

``` r
library(nichenetr)
library(dplyr)

weighted_networks = readRDS(url("https://zenodo.org/record/1484138/files/weighted_networks.rds"))
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/1484138/files/ligand_tf_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/1484138/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/1484138/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/1484138/files/gr_network.rds"))
```

As example, we will infer signaling paths between TGFB3 and its top-predicted target genes LAMA3, LAMC2 and TNC.

``` r
ligands_all = "TGFB3"
targets_all = c("LAMA3","LAMC2","TNC")

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix, ligands_all, targets_all, top_n_regulators = 3, weighted_networks, ligands_position = "cols")

# normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = 2*((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = 2*((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(active_signaling_network_min_max, ligands_all,targets_all, sig_color = "indianred", gr_color = "steelblue")
# DiagrammeR::render_graph(graph_min_max, layout = "tree")
```

We will now look which of the collected data sources support the interactions in this network.

``` r
data_source_network = infer_supporting_datasources(active_signaling_network,lr_network, sig_network , gr_network)
head(data_source_network) 
## # A tibble: 6 x 5
##   from  to    source                         database              layer   
##   <chr> <chr> <chr>                          <chr>                 <chr>   
## 1 SMAD1 LAMA3 harmonizome_TRANSFAC_CUR       harmonizome_gr        regulat~
## 2 SMAD1 LAMA3 pathwaycommons_controls_expre~ pathwaycommons_expre~ regulat~
## 3 SMAD2 LAMC2 harmonizome_CHEA               harmonizome_gr        regulat~
## 4 SMAD2 TNC   harmonizome_CHEA               harmonizome_gr        regulat~
## 5 SMAD2 TNC   regnetwork_source              regnetwork            regulat~
## 6 SMAD3 LAMA3 harmonizome_CHEA               harmonizome_gr        regulat~
```
