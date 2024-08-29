Evaluation of NicheNet’s ligand-target predictions
================
Robin Browaeys
2018-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/model_evaluation.Rmd", output_format = "github_document")
-->

This vignette shows how the ligand-target predictions of NicheNet were
evaluated. For validation, we collected transcriptome data of cells
before and after they were treated by one or two ligands in culture.
Using these ligand treatment datasets for validation has the advantage
that observed gene expression changes can be directly attributed to the
addition of the ligand(s). Hence, differentially expressed genes can be
considered as a gold standard of target genes of a particular ligand.

You can use the procedure shown here to evaluate your own model and
compare its performance to NicheNet. In NicheNet v2, we added more
ligand treatment validation datasets from CytoSig. We will also
demonstrate the better performance of NicheNet v2’s ligand-target
matrix. [Ligand treatment validation datasets, NicheNet’s v1
ligand-target model](https://doi.org/10.5281/zenodo.3260758), and
[NicheNet’s v2 ligand-target
model](https://doi.org/10.5281/zenodo.7074290) can be downloaded from
Zenodo.

### Load nichenetr, the model we want to evaluate, and the datasets on which we want to evaluate it.

``` r
library(nichenetr)
library(tidyverse)

# Load in the ligand-target model
ligand_target_matrix_v1 = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix_v2 = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))

# The ligand treatment expression datasets used for validation can be downloaded from Zenodo:
expression_settings_validation = readRDS(url("https://zenodo.org/record/8010790/files/expression_settings"))
cytosig_settings = readRDS(url("https://zenodo.org/record/8010790/files/cytosig_signatures_settings_loose")) # new evaluation dataset

#Ligand treatment datasets show the log fold change in expression of genes after treatment with one or more specific ligands. Here: example for the ligand NODAL:
head(expression_settings_validation$nodal_Nodal$diffexp)
## # A tibble: 6 × 3
##     lfc     qval gene  
##   <dbl>    <dbl> <chr> 
## 1 -4.44 5.33e-15 IGFBP5
## 2 -4.03 6.58e-16 NID1  
## 3 -3.97 1.96e-17 CLPS  
## 4  3.95 1.97e-15 CLDN6 
## 5  3.45 1.07e-18 GPC3  
## 6  3.20 1.16e-16 CLIC6
```

### Example: transcriptional response prediction evaluation

First, we will demonstrate how to evaluate the transcriptional response
(i.e. target gene prediction) performance for all ligand treatment
expression datasets. For this, we determine how well the model predicts
which genes are differentially expressed after treatment with a ligand.
Ideally, target genes with high regulatory potential scores for a
ligand, should be differentially expressed in response to that ligand.

For information of all collected ligand treatment datasets, see [Dataset
information](evaluation_datasets.xlsx)

For the sake of simplicity, we exclude in this vignette the
ligand-treatment datasets profiling the response to multiple ligands. To
see how to build a ligand-target model with target predictions for
multiple ligands at once: see vignette [Construction of NicheNet’s
ligand-target model](model_construction.md):
`vignette("model_construction", package="nichenetr")`.

Step 1: convert expression datasets to the required format to perform
target gene prediction

``` r
settings = expression_settings_validation %>% lapply(convert_expression_settings_evaluation)
settings = settings %>% discard(~length(.$from) > 1)
```

Step 2: calculate the target gene prediction performances

``` r
# Evaluate transcriptional response prediction on every dataset
performances = settings %>% lapply(evaluate_target_prediction, ligand_target_matrix_v2) %>% bind_rows()
```

Step 3: visualize the results: show here different classification
evaluation metrics
(Note: Mean-rank gene-set enrichment will only be calculated if limma is installed)

``` r
# Visualize some classification evaluation metrics showing the target gene prediction performance
performances = performances %>% select(-aupr, -auc_iregulon,-pearson_log_pval,-spearman_log_pval ,-sensitivity_roc, -specificity_roc) %>% gather(key = scorename, value = scorevalue, auroc:spearman)

scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)", auc_iregulon_corrected = "AUC-iRegulon (corrected)",pearson = "Pearson correlation", spearman = "Spearman's rank correlation",mean_rank_GST_log_pval = "Mean-rank gene-set enrichment")
scorerandom = c(auroc=0.5, aupr_corrected=0, auc_iregulon_corrected = 0, pearson = 0, spearman = 0,mean_rank_GST_log_pval = 0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

performances %>%
  mutate(model = "NicheNet v2") %>%
  ggplot() +
  geom_violin(aes(model, scorevalue, group=model, fill = model)) +
  geom_boxplot(aes(model, scorevalue, group = model),width = 0.05) +
  scale_y_continuous("Score target prediction") +
  facet_wrap(~scorename, scales = "free", labeller=as_labeller(scorelabels)) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, color = "red") +
  theme_bw()
```

![](model_evaluation_files/figure-gfm/target-prediction-v2-results-1.png)<!-- -->

We will now compare performances between NicheNet v1 and v2 on both
ligand treatment datasets. Note that although the performance of v2 is
much better here, the CytoSig experiments were also included during
model construction of v2. To get the results in the MultiNicheNet paper,
you will have to follow the `model_construction.Rmd` vignette and filter
out the CytoSig data sources during model construction.

``` r

performances_df <- lapply(c("nichenet_gs", "cytosig_gs"), function(gs) {
  
  # Choose 'settings' based on type of gold standard
  if (gs == "nichenet_gs") {
    settings = expression_settings_validation %>% lapply(convert_expression_settings_evaluation)
    settings = settings %>% discard(~length(.$from) > 1)
  } else {
    settings = cytosig_settings
  }
  
  lapply(c("v1", "v2"), function(ver){
    # Get the ligand_target_matrix according to the version, Evaluate transcriptional response prediction on every dataset
    performances = settings %>% lapply(evaluate_target_prediction, get(paste0("ligand_target_matrix_", ver))) %>% bind_rows()
  
    # Select some classification evaluation metrics showing the target gene prediction performance
    performances = performances %>% select(-aupr, -auc_iregulon,-pearson_log_pval,-spearman_log_pval ,-sensitivity_roc, -specificity_roc) %>%
      pivot_longer(auroc:spearman, names_to = "scorename", values_to = "scorevalue") %>% mutate(dataset = gs, version = ver)
  
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)


scorelabels = str_wrap(c("AUROC", "AUPR (corrected)", "Mean-rank gene-set enrichment", "AUC-iRegulon (corrected)","Pearson correlation", "Spearman's rank correlation"), width=10) %>%
  setNames(performances$scorename %>% unique)
scorerandom = c(auroc=0.5, aupr_corrected=0, auc_iregulon_corrected = 0, pearson = 0, spearman = 0,mean_rank_GST_log_pval = 0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

ggplot(performances_df, aes(y=scorevalue, x=version)) +
  geom_violin(aes(fill=version)) +
  geom_boxplot(width = 0.05) +
  labs(y="Target prediction score", x="NicheNet version") +
  facet_grid(scorename~dataset, scales = "free",
             labeller=labeller(scorename=scorelabels,
                               dataset=c(cytosig_gs="CytoSig datasets", nichenet_gs="NicheNet datasets"))) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, linewidth=0.2, color = "red") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        legend.position = "none")
```

![](model_evaluation_files/figure-gfm/target-prediction-comparison-1.png)<!-- -->

### Example: ligand activity prediction evaluation

Now we will show how to assess the accuracy of the model in predicting
whether cells were treated by a particular ligand or not. In other
words, we will evaluate how well NicheNet prioritizes active ligand(s),
given a set of differentially expressed genes. For this procedure, we
assume the following: the better a ligand predicts the transcriptional
response compared to other ligands, the more likely it is that this
ligand is active. Therefore, we first get ligand activity (or ligand
importance or feature importance) scores for all ligands on all
ligand-treatment expression datasets of which the true acive ligand is
known. Then we assess whether the truly active ligands get indeed higher
ligand activity scores as should be for a good ligand-target model.

A graphical summary of this procedure is visualized here below:

![](images/ligand_activity_prediction_workflow_new.png)

Step 1: convert expression datasets to the required format to perform
ligand activity prediction

``` r
# convert expression datasets to correct format for ligand activity prediction
settings = expression_settings_validation %>% lapply(convert_expression_settings_evaluation)
settings = settings %>% discard(~length(.$from) > 1)
all_ligands = settings %>% extract_ligands_from_settings(combination = FALSE) %>% unlist()
settings_ligand_prediction = settings %>% convert_settings_ligand_prediction(all_ligands = all_ligands, validation = TRUE)
```

Step 2: calculate the ligand importances (these are classification
evaluation metrics indicating how well a ligand can predict the observed
DE genes in a specific ligand treatment dataset)

``` r
# infer ligand importances: for all ligands of interest, we assess how well a ligand explains the differential expression in a specific datasets (and we do this for all datasets).
ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances,ligand_target_matrix_v2) %>% bind_rows()
```

Step 3: evaluate how separate ligand importances can predict ligand
activity

``` r
# Look at predictive performance of single/individual importance measures to predict ligand activity: of all ligands tested, the ligand that is truly active in a dataset should get the highest activity score (i.e. best target gene prediction performance)

# Replace infinite values with 10000
ligand_importances = ligand_importances %>% mutate(across(ends_with("log_pval"), ~ifelse(is.infinite(.x), 10000, .x)))

evaluation_ligand_prediction = ligand_importances$setting %>% unique() %>% lapply(function(x){x}) %>%
    lapply(wrapper_evaluate_single_importances_ligand_prediction,ligand_importances) %>%
    bind_rows() %>% inner_join(ligand_importances %>% distinct(setting,ligand))
```

Step 4: visualize the results: show here different classification
evaluation metrics

``` r
# Visualize some classification evaluation metrics showing the ligand activity prediction performance
evaluation_ligand_prediction = evaluation_ligand_prediction %>% select(-aupr, -sensitivity_roc, -specificity_roc, -pearson, -spearman, -mean_rank_GST_log_pval) %>% gather(key = scorename, value = scorevalue, auroc:aupr_corrected)
scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)")
scorerandom = c(auroc=0.5, aupr_corrected=0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

evaluation_ligand_prediction %>%
 filter(importance_measure %in% c("auroc", "aupr_corrected", "mean_rank_GST_log_pval", "auc_iregulon_corrected", "pearson", "spearman")) %>%
  ggplot() +
  geom_violin(aes(importance_measure, scorevalue, group=importance_measure, fill = importance_measure)) +
  geom_boxplot(aes(importance_measure, scorevalue, group = importance_measure),width = 0.1) +
  scale_y_continuous("Evaluation ligand activity prediction") +
  scale_x_discrete("Ligand activity measure") +
  facet_wrap(~scorename, scales = "free", labeller=as_labeller(scorelabels)) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
```

![](model_evaluation_files/figure-gfm/ligand-activity-v2-results-1.png)<!-- -->

We will again compare performances between NicheNet v1 and v2 on both
ligand treatment datasets.

``` r

performances_df <- lapply(c("nichenet_gs", "cytosig_gs"), function(gs) {
  
  # Choose 'settings' based on type of gold standard
  if (gs == "nichenet_gs") {
    settings = expression_settings_validation %>% lapply(convert_expression_settings_evaluation)
    settings = settings %>% discard(~length(.$from) > 1)
  } else {
    settings = cytosig_settings
  }
  
  all_ligands = settings %>% extract_ligands_from_settings(combination = FALSE) %>% unlist()
  settings_ligand_prediction = settings %>% convert_settings_ligand_prediction(all_ligands = all_ligands, validation = TRUE)
  
  lapply(c("v1", "v2"), function(ver){
    
    # Look at predictive performance of single/individual importance measures to predict ligand activity: of all ligands tested,
    # the ligand that is truly active in a dataset should get the highest activity score (i.e. best target gene prediction performance)
    # infer ligand importances: for all ligands of interest, we assess how well a ligand explains the differential expression in a specific datasets (and we do this for all datasets).
    ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances, get(paste0("ligand_target_matrix_", ver))) %>% bind_rows()
    
    # Filter out column with infinite values
    ligand_importances = ligand_importances %>% mutate(across(ends_with("log_pval"), ~ifelse(is.infinite(.x), 10000, .x)))

    evaluation_ligand_prediction = ligand_importances$setting %>% unique() %>% lapply(function(x){x}) %>%
      lapply(wrapper_evaluate_single_importances_ligand_prediction,ligand_importances) %>%
      bind_rows() %>% inner_join(ligand_importances %>% distinct(setting,ligand))

    # Select some classification evaluation metrics showing the ligand activity prediction performance
    evaluation_ligand_prediction %>% select(-aupr, -sensitivity_roc, -specificity_roc, -pearson, -spearman, -mean_rank_GST_log_pval) %>% pivot_longer(auroc:aupr_corrected, names_to = "scorename", values_to = "scorevalue") %>% mutate(dataset = gs, version = ver)
  
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)
```

``` r
scorelabels = str_wrap(c("AUROC", "AUPR", "AUPR     (corrected)", "Sensitivity ROC", "Specificity ROC",
                         "Mean-rank gene-set enrichment", "AUC-iRegulon", "AUC-iRegulon (corrected)",
                        "Pearson log p-val", "Spearman log p-val", "Pearson correlation", "Spearman's rank correlation"), width=15) %>% setNames(performances_df$importance_measure %>% unique)
scorerandom = c(auroc=0.5, aupr_corrected=0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

ggplot(performances_df %>% filter(importance_measure %in%
                                                 c("auroc", "aupr_corrected", "mean_rank_GST_log_pval",
                                                   "auc_iregulon_corrected", "pearson", "spearman")),
                                               aes(y=scorevalue, x=version, color=importance_measure, group=interaction(importance_measure, version))) +
  geom_violin(position=position_dodge(0.75)) +
  geom_boxplot(width = 0.05, position=position_dodge(0.75)) +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  scale_fill_discrete(labels=scorelabels) +
  guides(color="none") +
  labs(y="Evaluation ligand activity prediction", x="Ligand activity measure for two NicheNet versions",
       fill = "Importance\nmeasure") +
  facet_grid(scorename~dataset, scales = "free",
             labeller=labeller(scorename=scorelabels,
                               dataset=c(cytosig_gs="CytoSig datasets", nichenet_gs="NicheNet datasets"))) +
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, linewidth=0.2, color = "red") +
  theme_bw() +
  theme(strip.text.y = element_text(angle=0),
        strip.background.x = element_blank(),
        strip.background.y = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        legend.position = "bottom")
```

![](model_evaluation_files/figure-gfm/ligand-activity-comparison-1.png)<!-- -->
