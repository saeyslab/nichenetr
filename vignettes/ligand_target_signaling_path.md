Inferring ligand-to-target signaling paths
================
Robin Browaeys
2018-11-12

### Infer signaling paths beween ligand(s) and target(s) of interest

To determine signaling paths between a ligand and target of interest, we look at which transcription factors are best regulating the target genes and are most closely downstream of the ligand. Then, the shortest paths between these transcription factors and the ligand of interests are determined and genes forming part in this path are considered as important signaling mediators. Finally, we give all interactions between the ligand, signaling mediators, transcription factors and target genes.

First, we will load the necessary networks to infer signaling paths between ligand and target genes of interest.

``` r
library(nichenetr)
library(tidyverse)

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
DiagrammeR::render_graph(graph_min_max, layout = "tree")
```

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-22ef8b264c4c9aa55410">{"x":{"diagram":"digraph {\n\ngraph [layout = \"neato\",\n       outputorder = \"edgesfirst\",\n       bgcolor = \"white\"]\n\nnode [fontname = \"Helvetica\",\n      fontsize = \"10\",\n      shape = \"circle\",\n      fixedsize = \"true\",\n      width = \"0.5\",\n      style = \"filled\",\n      fillcolor = \"aliceblue\",\n      color = \"gray70\",\n      fontcolor = \"gray50\"]\n\nedge [fontname = \"Helvetica\",\n     fontsize = \"8\",\n     len = \"1.5\",\n     color = \"gray80\",\n     arrowsize = \"0.5\"]\n\n  \"1\" [label = \"SMAD1\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"13.5,5!\"] \n  \"2\" [label = \"SMAD2\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"21,4!\"] \n  \"3\" [label = \"SMAD3\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"16,3!\"] \n  \"4\" [label = \"SMAD4\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"13,2!\"] \n  \"5\" [label = \"TGFB3\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#CD5C5C\", pos = \"15,8!\"] \n  \"6\" [label = \"TGFBR1\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"15,7!\"] \n  \"7\" [label = \"TGFBR2\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"15.5,6!\"] \n  \"8\" [label = \"LAMA3\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#4682B4\", pos = \"0.5,1!\"] \n  \"9\" [label = \"LAMC2\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#4682B4\", pos = \"13,1!\"] \n  \"10\" [label = \"TNC\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#4682B4\", pos = \"19.5,1!\"] \n\"1\"->\"2\" [penwidth = \"1.2165371328981\", color = \"indianred\"] \n\"1\"->\"3\" [penwidth = \"1.62943194594475\", color = \"indianred\"] \n\"1\"->\"4\" [penwidth = \"1.75740609306738\", color = \"indianred\"] \n\"2\"->\"1\" [penwidth = \"1.22151815883893\", color = \"indianred\"] \n\"2\"->\"3\" [penwidth = \"2.11191973227599\", color = \"indianred\"] \n\"2\"->\"4\" [penwidth = \"1.81304401432004\", color = \"indianred\"] \n\"2\"->\"6\" [penwidth = \"1.27587941235455\", color = \"indianred\"] \n\"2\"->\"7\" [penwidth = \"1.39789942561386\", color = \"indianred\"] \n\"3\"->\"1\" [penwidth = \"1.63617061797546\", color = \"indianred\"] \n\"3\"->\"2\" [penwidth = \"1.77553033990597\", color = \"indianred\"] \n\"3\"->\"4\" [penwidth = \"1.82601259288889\", color = \"indianred\"] \n\"3\"->\"6\" [penwidth = \"1.35387343185849\", color = \"indianred\"] \n\"3\"->\"7\" [penwidth = \"1.27755682515373\", color = \"indianred\"] \n\"4\"->\"1\" [penwidth = \"1.9284714460309\", color = \"indianred\"] \n\"4\"->\"2\" [penwidth = \"1.9732496749038\", color = \"indianred\"] \n\"4\"->\"3\" [penwidth = \"1.97534850369824\", color = \"indianred\"] \n\"5\"->\"1\" [penwidth = \"0.872165083151721\", color = \"indianred\"] \n\"5\"->\"2\" [penwidth = \"1.44885535243043\", color = \"indianred\"] \n\"5\"->\"3\" [penwidth = \"1.45015215990848\", color = \"indianred\"] \n\"5\"->\"4\" [penwidth = \"0.75\", color = \"indianred\"] \n\"5\"->\"6\" [penwidth = \"1.90723694124846\", color = \"indianred\"] \n\"5\"->\"7\" [penwidth = \"1.81731895449705\", color = \"indianred\"] \n\"6\"->\"1\" [penwidth = \"2.14267598976886\", color = \"indianred\"] \n\"6\"->\"2\" [penwidth = \"2.75\", color = \"indianred\"] \n\"6\"->\"3\" [penwidth = \"2.20854831149134\", color = \"indianred\"] \n\"6\"->\"4\" [penwidth = \"1.46362943457838\", color = \"indianred\"] \n\"6\"->\"7\" [penwidth = \"1.65274103743568\", color = \"indianred\"] \n\"7\"->\"1\" [penwidth = \"0.944002634602939\", color = \"indianred\"] \n\"7\"->\"2\" [penwidth = \"1.91365626280769\", color = \"indianred\"] \n\"7\"->\"3\" [penwidth = \"2.13118795313037\", color = \"indianred\"] \n\"7\"->\"4\" [penwidth = \"1.38610258052106\", color = \"indianred\"] \n\"7\"->\"6\" [penwidth = \"2.33636391701077\", color = \"indianred\"] \n\"1\"->\"8\" [penwidth = \"1.80913614222845\", color = \"steelblue\"] \n\"2\"->\"9\" [penwidth = \"0.805311951308608\", color = \"steelblue\"] \n\"2\"->\"10\" [penwidth = \"0.967941587850231\", color = \"steelblue\"] \n\"3\"->\"8\" [penwidth = \"0.780306278334167\", color = \"steelblue\"] \n\"3\"->\"9\" [penwidth = \"2.40065558831022\", color = \"steelblue\"] \n\"3\"->\"10\" [penwidth = \"2.75\", color = \"steelblue\"] \n\"4\"->\"8\" [penwidth = \"1.73102851351999\", color = \"steelblue\"] \n\"4\"->\"9\" [penwidth = \"2.40065558831022\", color = \"steelblue\"] \n\"4\"->\"10\" [penwidth = \"1.59931117765936\", color = \"steelblue\"] \n\"5\"->\"10\" [penwidth = \"0.876701313520233\", color = \"steelblue\"] \n\"7\"->\"8\" [penwidth = \"0.75\", color = \"steelblue\"] \n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>
<!--/html_preserve-->
We will now look which of the collected data sources support the interactions in this network.

``` r
data_source_network = infer_supporting_datasources(active_signaling_network,lr_network, sig_network , gr_network)
data_source_network 
## # A tibble: 249 x 5
##    from  to    source                        database             layer   
##    <chr> <chr> <chr>                         <chr>                <chr>   
##  1 SMAD1 LAMA3 harmonizome_TRANSFAC_CUR      harmonizome_gr       regulat~
##  2 SMAD1 LAMA3 pathwaycommons_controls_expr~ pathwaycommons_expr~ regulat~
##  3 SMAD2 LAMC2 harmonizome_CHEA              harmonizome_gr       regulat~
##  4 SMAD2 TNC   harmonizome_CHEA              harmonizome_gr       regulat~
##  5 SMAD2 TNC   regnetwork_source             regnetwork           regulat~
##  6 SMAD3 LAMA3 harmonizome_CHEA              harmonizome_gr       regulat~
##  7 SMAD3 LAMC2 harmonizome_CHEA              harmonizome_gr       regulat~
##  8 SMAD3 LAMC2 harmonizome_TRANSFAC_CUR      harmonizome_gr       regulat~
##  9 SMAD3 LAMC2 pathwaycommons_controls_expr~ pathwaycommons_expr~ regulat~
## 10 SMAD3 TNC   harmonizome_CHEA              harmonizome_gr       regulat~
## # ... with 239 more rows
```
