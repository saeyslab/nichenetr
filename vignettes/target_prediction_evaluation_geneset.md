Assess how well top-ranked ligands can predict a gene set of interest
================
Robin Browaeys
2019-02-19

<!-- github markdown built using
rmarkdown::render("vignettes/target_prediction_evaluation_geneset.Rmd", output_format = "github_document")
-->

This vignette assesses the ligands prioritized by NicheNet in their
ability to predict a gene set of interest. We will first follow the
steps of [Perform NicheNet analysis starting from a Seurat
object](seurat_wrapper.md)) to obtain ligands rankings. Make sure you
understand the steps and output of a basic NicheNet analysis (more
information in [Perform NicheNet analysis starting from a Seurat object:
step-by-step analysis](seurat_steps.md). You can also apply this
vignette to the [NicheNet’s ligand activity analysis on a gene set of
interest](ligand_activity_geneset.md) vignette.

## Run NicheNet

``` r
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
```

``` r
# Read in networks
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

lr_network <- lr_network %>% distinct(from, to)

# Read in expression data and update
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

# Run NicheNet
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), 
  receiver = "CD8 T", 
  condition_colname = "aggregate",
  condition_oi = "LCMV",
  condition_reference = "SS",
  expression_pct = 0.05,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
  )
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"

best_upstream_ligands <- nichenet_output$ligand_activities %>%
  top_n(30, aupr_corrected) %>% arrange(desc(aupr_corrected)) %>% pull(test_ligand)
```

## Assess how well top-ranked ligands can predict a gene set of interest

For the top 30 ligands, we will now build a multi-ligand model that uses
all top-ranked ligands to predict whether a gene belongs to the gene set
of interest (differentially expressed genes in CD8 T cells after LCMV
infection) or not. This classification model will be trained via
cross-validation and returns a probability for every gene.

``` r
# change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
k <- 3 # 3-fold
n <- 2 # 2 rounds

gene_predictions_top30_list <- lapply(1:n, assess_rf_class_probabilities,
                                      folds = k,
                                      geneset = nichenet_output$geneset_oi,
                                      background_expressed_genes = nichenet_output$background_expressed_genes,
                                      ligands_oi = best_upstream_ligands,
                                      ligand_target_matrix = ligand_target_matrix)
```

Evaluate now how well the target gene probabilities accord to the gene
set assignments.

``` r
# get performance: auroc-aupr-pearson
target_prediction_performances_cv <- gene_predictions_top30_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>%
  bind_rows() %>% mutate(round=seq(1:nrow(.)))
```

What is the AUROC, AUPR and PCC of this model (averaged over
cross-validation rounds)?

``` r
target_prediction_performances_cv$auroc %>% mean()
## [1] 0.8044756
target_prediction_performances_cv$aupr %>% mean()
## [1] 0.4975771
target_prediction_performances_cv$pearson %>% mean()
## [1] 0.541401
```

Evaluate now whether genes belonging to the gene set are more likely to
be top-predicted. We will look at the top 5% of predicted targets here.

``` r
# get performance: how many viral response genes and non-viral response-genes among top 5% predicted targets
target_prediction_performances_discrete_cv <- gene_predictions_top30_list %>%
  lapply(calculate_fraction_top_predicted,
         quantile_cutoff = 0.95) %>%
  bind_rows(.id = "round") 
```

What is the fraction of viral response genes that belongs to the top 5%
predicted targets?

``` r
target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean()
## [1] 0.45
```

What is the fraction of non-viral-response genes that belongs to the top
5% predicted targets?

``` r
target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean()
## [1] 0.0179845
```

We see that the viral response genes are enriched in the top-predicted
target genes. To test this, we will now apply a Fisher’s exact test for
every cross-validation round and report the average p-value.

``` r
target_prediction_performances_discrete_fisher <- gene_predictions_top30_list %>%
  lapply(calculate_fraction_top_predicted_fisher,
         quantile_cutoff = 0.95)

target_prediction_performances_discrete_fisher %>% unlist() %>% mean()
## [1] 2.332346e-96
```

Finally, we will look at which p-EMT genes are well-predicted in every
cross-validation round.

``` r
# get top predicted genes
top_predicted_genes <- lapply(1:n, get_top_predicted_genes,
                              gene_predictions_top30_list) %>%
  reduce(full_join, by = c("gene","true_target"))

top_predicted_genes %>% filter(true_target)
## # A tibble: 125 × 4
##    gene   true_target predicted_top_target_round1 predicted_top_target_round2
##    <chr>  <lgl>       <lgl>                       <lgl>                      
##  1 Gbp4   TRUE        TRUE                        TRUE                       
##  2 Gbp9   TRUE        TRUE                        TRUE                       
##  3 Ifi203 TRUE        TRUE                        TRUE                       
##  4 Ifi209 TRUE        TRUE                        TRUE                       
##  5 Ifi213 TRUE        TRUE                        TRUE                       
##  6 Ifi208 TRUE        TRUE                        TRUE                       
##  7 Mndal  TRUE        TRUE                        TRUE                       
##  8 Ifi206 TRUE        TRUE                        TRUE                       
##  9 Phf11  TRUE        TRUE                        TRUE                       
## 10 Ifit3b TRUE        TRUE                        TRUE                       
## # ℹ 115 more rows
```
