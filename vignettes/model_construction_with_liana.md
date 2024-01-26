Using LIANA ligand-receptor databases to construct the ligand-target
model
================
Chananchida Sang-aram
2023-08-23

<!-- github markdown built using
rmarkdown::render("vignettes/model_construction_with_liana.Rmd", output_format = "github_document")
-->

Following the [Construction of NicheNet’s ligand-target
model](model_construction.md) vignette, we will now demonstrate how to
use ligand-receptor reactions from LIANA to build the ligand-target
model. LIANA is a framework that combines both resources and
computational tools for ligand-receptor cell-cell communication
inference (Dimitrov et al., 2022). As the NicheNet prior model is built
by integrating ligand-receptor, signaling, and gene regulatory
databases, each part can be replaced with external data sources. We will
show how the first part, the ligand-receptor database, can be replaced
with those from LIANA, and how to run the model afterward.

**Important**: Since LIANA also offers functions to calculate
ligand-receptor interactions of interest, it is also possible to use
them to select which ligands are of interest to do the ligand activity
analysis. This is explained further in the [LIANA
vignette](https://saezlab.github.io/liana/articles/liana_nichenet.html).

First, we will install LIANA:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github('saezlab/liana')
```

Load necessary packages.

``` r
library(liana)
library(nichenetr)
library(tidyverse)
library(Seurat)
```

To check which resources are present in LIANA, we can use the
`show_resources()` function. These are then accessed via
`select_resource()`.

``` r
show_resources()
##  [1] "Default"          "Consensus"        "Baccin2019"       "CellCall"         "CellChatDB"       "Cellinker"        "CellPhoneDB"     
##  [8] "CellTalkDB"       "connectomeDB2020" "EMBRACE"          "Guide2Pharma"     "HPMR"             "ICELLNET"         "iTALK"           
## [15] "Kirouac2010"      "LRdb"             "Ramilowski2015"   "OmniPath"         "MouseConsensus"
```

Next, we will calculate how much overlap there is between the ligands
and receptors in the LIANA and NicheNet databases. If the overlap
between LIANA receptors and NicheNet signaling network is too low, the
integration will probably not work very well. Furthermore, The
`decomplexify()` function of LIANA is crucial in our case, as we would
like to separate receptors into their respective subunits.

``` r
# Load signaling networks of NicheNet
lr_network_human <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network_mouse <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))

sig_network_human <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
sig_network_mouse <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_mouse_21122021.rds"))

overlap_df <- lapply(show_resources()[-1], function(resource) {
  db <- select_resource(resource)[[1]] %>% decomplexify()

  lr_network <- paste0("lr_network_", ifelse(grepl("Mouse", resource), "mouse", "human"))
  sig_network <- paste0("sig_network_", ifelse(grepl("Mouse", resource), "mouse", "human"))
  
  data.frame(row.names = resource,
             n_ligands = length(unique(db$source_genesymbol)),
             n_receptors = length(unique(db$target_genesymbol)),
             n_ligands_overlap = length(intersect(db$source_genesymbol, get(lr_network)$from)),
             n_receptors_overlap_lr = length(intersect(db$target_genesymbol, get(lr_network)$to)),
             n_receptors_overlap_sig = length(intersect(db$target_genesymbol, get(sig_network)$from))
  ) %>% mutate(frac_ligands_overlap = n_ligands_overlap / n_ligands,
               frac_receptors_overlap_lr = n_receptors_overlap_lr / n_receptors,
               frac_receptors_overlap_sig = n_receptors_overlap_sig / n_receptors)
  
}) %>% do.call(rbind, .)

overlap_df
##                  n_ligands n_receptors n_ligands_overlap n_receptors_overlap_lr n_receptors_overlap_sig frac_ligands_overlap frac_receptors_overlap_lr
## Consensus             1032         934               923                    794                     927            0.8943798                 0.8501071
## Baccin2019             650         612               550                    535                     611            0.8461538                 0.8741830
## CellCall               276         195               272                    193                     195            0.9855072                 0.9897436
## CellChatDB             531         450               515                    430                     447            0.9698682                 0.9555556
## Cellinker             1216        1010               879                    817                    1003            0.7228618                 0.8089109
## CellPhoneDB            500         461               479                    415                     458            0.9580000                 0.9002169
## CellTalkDB             812         780               741                    683                     776            0.9125616                 0.8756410
## connectomeDB2020       813         681               768                    640                     678            0.9446494                 0.9397944
## EMBRACE                462         473               435                    436                     473            0.9415584                 0.9217759
## Guide2Pharma           293         243               289                    242                     243            0.9863481                 0.9958848
## HPMR                   365         492               260                    407                     492            0.7123288                 0.8272358
## ICELLNET               318         252               313                    246                     252            0.9842767                 0.9761905
## iTALK                  714         699               666                    625                     697            0.9327731                 0.8941345
## Kirouac2010            143          81               135                     80                      81            0.9440559                 0.9876543
## LRdb                   801         745               726                    662                     742            0.9063670                 0.8885906
## Ramilowski2015         639         589               607                    564                     588            0.9499218                 0.9575552
## OmniPath              1219         907               956                    773                     901            0.7842494                 0.8522602
## MouseConsensus         874         796               763                    685                     789            0.8729977                 0.8605528
##                  frac_receptors_overlap_sig
## Consensus                         0.9925054
## Baccin2019                        0.9983660
## CellCall                          1.0000000
## CellChatDB                        0.9933333
## Cellinker                         0.9930693
## CellPhoneDB                       0.9934924
## CellTalkDB                        0.9948718
## connectomeDB2020                  0.9955947
## EMBRACE                           1.0000000
## Guide2Pharma                      1.0000000
## HPMR                              1.0000000
## ICELLNET                          1.0000000
## iTALK                             0.9971388
## Kirouac2010                       1.0000000
## LRdb                              0.9959732
## Ramilowski2015                    0.9983022
## OmniPath                          0.9933848
## MouseConsensus                    0.9912060
```

On average, ~90% of the ligands and receptors of LIANA databases are in
the NicheNet LR network (`frac_ligands_overlap`,
`frac_receptors_overlap_lr`), and almost all of the receptors in LIANA
databases are present in the NicheNet signaling network
(`frac_receptors_overlap_sig`). When using the “Consensus” database of
LIANA, there are ~100 ligands that are not present in NicheNet; in
contrast, there are 303 ligands in NicheNet that are not present in the
LIANA consensus database.

To build the ligand-target model, we can use a very similar code to the
[Construction of NicheNet’s ligand-target model](model_construction.md)
vignette. Users can choose between replacing the NicheNet LR database
entirely with LIANA’s (`replace_nichenet_lr = TRUE`), or just adding the
LIANA database as an additional data source, which may contain a lot of
redundant information.

``` r
# Load liana consensus LR network, and rename columns to be compatible with nichenet function
# Users can choose here whether or not to replace the NicheNet LR network entirely, or just add the LIANA db to it
replace_nichenet_lr <- TRUE
liana_db <- select_resource("Consensus")[[1]] %>% decomplexify()
liana_db <- liana_db %>% rename(from = source_genesymbol, to = target_genesymbol) %>% select(from, to) %>% mutate(source = "liana") %>%
  {if (replace_nichenet_lr) (.) else (bind_rows(lr_network_human, .))}

# Change source weights dataframe (but in this case all source weights are 1)
source_weights <- source_weights_df %>%
  add_row(source = "liana", weight = 1, .before = 1)

# Load the gene regulatory network
gr_network_human <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))

# Aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks <- construct_weighted_networks(lr_network = liana_db,
                                                 sig_network = sig_network_human,
                                                 gr_network = gr_network_human,
                                                 source_weights_df = source_weights)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks <- apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# in this example we will calculate target gene regulatory potential scores for TNF and the ligand combination TNF+IL6
ligands <- list("TNF",c("TNF","IL6"))
ligand_target_matrix_liana <- construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
```

### Running NicheNet on the LIANA ligand-target model

In this section compare the results between using the LIANA LR network
and NicheNet one in a typical [NicheNet analysis](seurat_steps.md). As
this is mouse scRNA-seq data, we will build the model using mouse
networks (“MouseConsensus” resource in LIANA). Furthermore, instead of
using the same source weights for all data sources as in the previous
section, we will use the optimized data source weights. Here, we will
use the same data source weight for the LIANA model as the one that was
computed for the NicheNet LR network. This may not be the optimum
weight, but it requires a lot of runtime to optimize this parameter (see
[Parameter optimization vignette](parameter_optimization.md) for more
information).

``` r
# Typical nichenet analysis
# Load seurat object
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse") # convert gene names

# Load LIANA mouse consensus ligand-receptor network
lr_network_liana <- select_resource("MouseConsensus")[[1]] %>% decomplexify()
lr_network_liana <- lr_network_liana %>% rename(from = source_genesymbol, to = target_genesymbol) %>% select(from, to) %>% mutate(source = "liana")

# Define receiver cell type
receiver = "CD8 T"
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

# Define sender cell types
sender_celltypes <- c("CD4 T","Treg", "Mono", "NK", "B", "DC")

# Get expressed genes in the sender
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
 
# Get gene set of interest through DE analysis
seurat_obj_receiver <- subset(seuratObj, idents = receiver)
seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver$aggregate)
condition_oi <- "LCMV"; condition_reference <- "SS"
DE_table_receiver <- FindMarkers(object = seurat_obj_receiver,
                                 ident.1 = condition_oi, ident.2 = condition_reference,
                                 min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)

# Get potential ligands
ligands <- lr_network_liana %>% pull(from) %>% unique()
receptors <- lr_network_liana %>% pull(to) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_sender)
expressed_receptors <- intersect(receptors,expressed_genes_receiver)
potential_ligands <- lr_network_liana %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands
##  [1] "Adam10" "Adam23" "Icosl"  "Thy1"   "H2-K1"  "App"    "Itgb2"  "Icam1"  "Lck"    "Fcer2a" "Cd48"   "Icam2"  "H2-M3"  "Vcam1"  "Cd22"   "Cd86"  
## [17] "B2m"    "Sirpa"  "Adam17" "Selplg" "Btla"   "F11r"   "Fn1"    "Ebi3"   "Hp"     "Tgfb1"  "Lgals1" "Vcan"   "Mif"    "Il15"   "Copa"   "Ybx1"  
## [33] "Fam3c"  "Ccl2"   "Crlf2"  "Ccl22"  "Cxcl10" "Cxcl16" "Cd72"   "Txlna"  "Psen1"  "F13a1"  "C1qb"   "Gnai2"  "Tgm2"   "Gnas"   "Ccl5"   "Mia"   
## [49] "H2-T23" "Lyz1"   "Lyz2"

# Constructing the ligand-target matrix
gr_network_mouse <- readRDS(url("https://zenodo.org/records/7074291/file/gr_network_mouse_21122021.rds"))

# Define optimum weight as the one calculated for the NicheNet LR network
optim_weight <- optimized_source_weights_df %>% filter(source == "nichenet_verschueren") %>% pull(avg_weight) # 0.2788104
source_weights <- optimized_source_weights_df %>%
  add_row(source = "liana", avg_weight = optim_weight, .before = 1) %>% rename(weight=avg_weight)

# Aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks_liana <- construct_weighted_networks(lr_network = lr_network_liana, sig_network = sig_network_mouse, gr_network = gr_network_mouse, source_weights_df = source_weights)

# Downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks_liana <- apply_hub_corrections(weighted_networks = weighted_networks_liana,
                                           lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                           gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# Construct ligand-target matrix using only the potential ligands to save time
ligand_target_matrix_liana <- construct_ligand_target_matrix(weighted_networks = weighted_networks_liana, ligands = as.list(potential_ligands), algorithm = "PPR",
                                                       damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                       ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))

# Filter background genes and the gene set of interest to only the ones in the ligand-target matrix
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix_liana)]
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix_liana)]

# Perform ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix_liana, potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
## # A tibble: 51 × 6
##    test_ligand auroc  aupr aupr_corrected pearson  rank
##    <chr>       <dbl> <dbl>          <dbl>   <dbl> <dbl>
##  1 Ebi3        0.663 0.389          0.244  0.300      1
##  2 H2-M3       0.608 0.292          0.147  0.179      2
##  3 H2-T23      0.611 0.278          0.133  0.153      3
##  4 Lck         0.621 0.268          0.122  0.103      4
##  5 H2-K1       0.605 0.268          0.122  0.141      5
##  6 Sirpa       0.613 0.263          0.117  0.143      6
##  7 Cd48        0.612 0.257          0.112  0.0860     7
##  8 Tgfb1       0.597 0.254          0.108  0.202      8
##  9 Ccl22       0.605 0.250          0.104  0.125      9
## 10 Ccl5        0.606 0.248          0.102  0.0343    10
## # ℹ 41 more rows
```

#### Compare results with NicheNet ligand-target model

Run NicheNet using the wrapper function.

``` r
# Use the aggregate function
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj,
  receiver = "CD8 T",
  condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS",
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network_mouse %>% distinct(from, to),
  weighted_networks = weighted_networks)
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

nichenet_output$ligand_activities
## # A tibble: 73 × 6
##    test_ligand auroc  aupr aupr_corrected pearson  rank
##    <chr>       <dbl> <dbl>          <dbl>   <dbl> <dbl>
##  1 Ebi3        0.663 0.390          0.244   0.301     1
##  2 Ptprc       0.642 0.310          0.165   0.167     2
##  3 H2-M3       0.608 0.292          0.146   0.179     3
##  4 H2-M2       0.611 0.279          0.133   0.153     5
##  5 H2-T10      0.611 0.279          0.133   0.153     5
##  6 H2-T22      0.611 0.279          0.133   0.153     5
##  7 H2-T23      0.611 0.278          0.132   0.153     7
##  8 H2-K1       0.605 0.268          0.122   0.142     8
##  9 H2-Q4       0.605 0.268          0.122   0.141    10
## 10 H2-Q6       0.605 0.268          0.122   0.141    10
## # ℹ 63 more rows
```

Compare the results between LIANA and NicheNet. Here we see that half of
the top 20 ligands between the two are the same.

``` r
intersect(ligand_activities[1:20,]$test_ligand, nichenet_output$ligand_activities[1:20,]$test_ligand)
##  [1] "Ebi3"   "H2-M3"  "H2-T23" "H2-K1"  "Sirpa"  "Cd48"   "Tgfb1"  "Ccl22"  "App"    "Selplg" "Cxcl10" "Btla"
```

Below is a comparison of the rankings. Some of the NicheNet top-ranked
ligands are not present in the LIANA LR network, such as Ptprc and some
H2- ligands. Nonetheless, the LIANA network seems to have made some new
links that results in new ligands appearing, such as Lck, Ccl5, and
Crlf2.

``` r
rankings_df <- bind_rows(ligand_activities %>% select(test_ligand, rank) %>% mutate(db = "LIANA"),
                         nichenet_output$ligand_activities %>% select(test_ligand, rank) %>% mutate(db = "NicheNet"))
rankings_df <- rankings_df %>% group_by(db) %>% mutate(new_rank = 1:n()) %>%
  group_by(db, rank) %>% mutate(ties = n() > 1, ties = factor(ties, levels = c(TRUE, FALSE))) %>%
  ungroup() %>% mutate(x_label = case_when(db == "NicheNet" ~ 0.8, db == "LIANA" ~ 2.2),
                       db = factor(db, levels = c("NicheNet", "LIANA")))

# Top 20
top_n <- 20
p1 <- ggplot(rankings_df %>% filter(rank <= top_n),
       aes(x=db, y=new_rank, label=test_ligand, group = test_ligand, color = ties)) +
  geom_point() +
  geom_line() +
  geom_text(aes(x = x_label), show.legend = FALSE) +
  scale_color_manual(values = c("tomato", "black")) +
  scale_y_reverse() + theme_classic() +
  labs(y = "Ligand rankings", color = "Tied rankings", title = "Top 20 ligands") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=10), legend.position = "none")

# All ligands
p2 <- ggplot(rankings_df, aes(x=db, y=new_rank, label=test_ligand, group = test_ligand, color = ties)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("tomato", "black")) +
  scale_y_reverse() + theme_classic() +
  labs(y = "Ligand rankings", color = "Tied rankings", title = "All ligands") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=10))


cowplot::plot_grid(p1, p2)
```

![](model_construction_with_liana_files/figure-gfm/line-plot-1.png)<!-- -->

## References

Dimitrov, D., Türei, D., Garrido-Rodriguez M., Burmedi P.L., Nagai,
J.S., Boys, C., Flores, R.O.R., Kim, H., Szalai, B., Costa, I.G.,
Valdeolivas, A., Dugourd, A. and Saez-Rodriguez, J. Comparison of
methods and resources for cell-cell communication inference from
single-cell RNA-Seq data. Nat Commun 13, 3224 (2022).
<https://doi.org/10.1038/s41467-022-30755-0>
