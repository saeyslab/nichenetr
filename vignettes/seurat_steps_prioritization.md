Prioritization of ligands based on expression values
================
Robin Browaeys & Chananchida Sang-aram
2023-10-02

<!-- github markdown built using 
rmarkdown::render("vignettes/seurat_steps_prioritization.Rmd", output_format = "github_document")
-->

In this vignette, we will extend the basic NicheNet analysis analysis
from [Perform NicheNet analysis starting from a Seurat object:
step-by-step analysis](seurat_steps.md) by incorporating gene expression
as part of the prioritization This is a generalization of the
[Differential NicheNet](differential_nichenet.md) and
[MultiNicheNet](https://github.com/saeyslab/multinichenetr) approach.
While the original NicheNet only ranks ligands based on the ligand
activity analysis, it is now also possible to prioritize ligands based
on cell type and condition specificity of the ligand and receptor.

We will again make use of mouse NICHE-seq data to explore intercellular
communication in the T cell area in the inguinal lymph node before and
72 hours after lymphocytic choriomeningitis virus (LCMV) infection
(Medaglia et al. 2017). The used [ligand-target
matrix](https://doi.org/10.5281/zenodo.7074290) and the [Seurat object
of the processed NICHE-seq single-cell
data](https://doi.org/10.5281/zenodo.3531889) can be downloaded from
Zenodo.

Make sure you understand the different steps in a NicheNet analysis that
are described in the basic vignette before proceeding with this
vignette.

# Prepare NicheNet analysis

Load required packages, read in the Seurat object with processed
expression data of interacting cells and NicheNet’s ligand-target prior
model, ligand-receptor network and weighted integrated networks.

``` r
library(nichenetr) # Please update to v2.0.6
library(Seurat)
library(SeuratObject)
library(tidyverse)
```

``` r
# Read Seurat object
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

# Load in networks
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

lr_network <- lr_network %>% distinct(from, to)
```

# Perform the NicheNet analysis

We will use the sender-focused approach here.

``` r
# 1. Define set of potential ligands
receiver <- "CD8 T"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

sender_celltypes <- c("CD4 T", "Treg", "Mono", "NK", "B", "DC")
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

all_ligands <- unique(lr_network$from)
all_receptors <- unique(lr_network$to)

expressed_ligands <- intersect(all_ligands, expressed_genes_sender)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>% unique()


# 2. Define the gene set of interest
condition_oi <-  "LCMV"
condition_reference <- "SS"

seurat_obj_receiver <- subset(seuratObj, idents = receiver)
seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["aggregate", drop=TRUE]])

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# 3. Define background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# 4. Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))
```

# Perform prioritization of ligand-receptor pairs

We will prioritize ligand-receptor pairs based on the following criteria
(with their corresponding weight names):

- Upregulation of the ligand in a sender cell type compared to other
  cell types: `de_ligand`
- Upregulation of the receptor in a receiver cell type: `de_receptor`
- Average expression of the ligand in the sender cell type:
  `exprs_ligand`
- Average expression of the receptor in the receiver cell type:
  `exprs_receptor`
- Condition-specificity of the ligand across all cell types:
  `ligand_condition_specificity`
- Condition-specificity of the receptor across all cell types:
  `receptor_condition_specificity`

This means that we will have to calculate:

- Differential expression of the ligand/receptor in a sender/receiver
  cell type
- The average expression of each ligand/receptor in each sender/receiver
  cell type
- Differential expression of the ligand/receptor between the two
  conditions

We provide a wrapper function `generate_info_tables` that will calculate
all these values for you. This function returns a list with three
dataframes:

- `sender_receiver_de`: differential expression of the ligand and
  receptor in the sender-receiver cell type pair. These were first
  calculated separately (i.e., DE of ligand in sender cell type, DE of
  receptor in receiver cell type based on FindAllMarkers) and then
  combined based on possible interactions from the lr_network.
- `sender_receiver_info`: the average expression of the ligand and
  receptor in sender-receiver cell type pairs
- `lr_condition_de`: differential expression of the ligand and receptor
  between the two conditions across all cell types.

Note that cell type specificity (i.e., the first four conditions) is
calculated only in the condition of interest.

The “scenario” argument can be either “case_control” or “one_condition”.
In “case_control” scenario, condition specificity is calculated.

``` r
lr_network_filtered <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)

info_tables <- generate_info_tables(seuratObj,
                                    celltype_colname = "celltype",
                                    senders_oi = sender_celltypes,
                                    receivers_oi = receiver,
                                    lr_network = lr_network_filtered,
                                    condition_colname = "aggregate",
                                    condition_oi = condition_oi,
                                    condition_reference = condition_reference,
                                    scenario = "case_control")

names(info_tables)
## [1] "sender_receiver_de"   "sender_receiver_info" "lr_condition_de"
info_tables$sender_receiver_de %>% head()
##   sender receiver ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg  p_val_ligand  p_adj_ligand p_val_receptor p_adj_receptor
## 1     DC    CD8 T   Ccl5    Cxcr3   6.432043   0.16714791                3.299595  1.893317e-25  2.563740e-21   7.758812e-05   1.000000e+00
## 2   Mono    CD8 T   Lyz2    Itgal   5.493265  -0.01687003                2.738198 1.728697e-160 2.340828e-156   4.973381e-02   1.000000e+00
## 3     DC    CD8 T  H2-M2     Cd8a   3.416479   1.94059972                2.678539 1.017174e-272 1.377355e-268  5.250531e-206  7.109745e-202
## 4     DC    CD8 T Cxcl16    Cxcr6   4.182085   0.54826454                2.365175 1.138617e-243 1.541801e-239   5.987787e-21   8.108063e-17
## 5   Mono    CD8 T  Cxcl9    Cxcr3   4.328801   0.16714791                2.247975 3.834954e-124 5.192911e-120   7.758812e-05   1.000000e+00
## 6   Mono    CD8 T  Cxcl9     Dpp4   4.328801   0.16416445                2.246483 3.834954e-124 5.192911e-120   6.628900e-04   1.000000e+00
##   pct_expressed_sender pct_expressed_receiver
## 1                1.000                  0.042
## 2                0.933                  0.188
## 3                0.429                  0.659
## 4                0.929                  0.089
## 5                0.547                  0.042
## 6                0.547                  0.148
info_tables$sender_receiver_info %>% head()
## # A tibble: 6 × 7
##   sender receiver ligand receptor avg_ligand avg_receptor ligand_receptor_prod
##   <chr>  <chr>    <chr>  <chr>         <dbl>        <dbl>                <dbl>
## 1 DC     Mono     B2m    Tap1           216.         8.59                1856.
## 2 DC     NK       B2m    Klrd1          216.         7.43                1607.
## 3 DC     B        B2m    Tap1           216.         7.35                1588.
## 4 DC     Treg     B2m    Tap1           216.         7.18                1552.
## 5 Mono   Mono     B2m    Tap1           158.         8.59                1353.
## 6 DC     DC       B2m    Tap1           216.         5.91                1277.
info_tables$lr_condition_de %>% head()
##   ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg  p_val_ligand  p_adj_ligand p_val_receptor p_adj_receptor
## 1 H2-Ab1      Cd4  2.4021254   0.11569357               1.2589095  4.424390e-06  5.991066e-02   5.634068e-02   1.000000e+00
## 2 Cxcl10     Dpp4  1.6066163   0.35175421               0.9791853  6.700636e-29  9.073332e-25   1.170731e-06   1.585287e-02
## 3    B2m     Tap1  0.7071427   1.13931050               0.9232266 6.936359e-174 9.392524e-170   3.585450e-52   4.855057e-48
## 4 H2-T22    Klrd1  1.5223370  -0.05659737               0.7328698 1.006291e-111 1.362618e-107   6.202530e-01   1.000000e+00
## 5 H2-T23    Klrd1  1.4651999  -0.05659737               0.7043013 1.789643e-114 2.423356e-110   6.202530e-01   1.000000e+00
## 6 Cxcl10    Cxcr3  1.6066163  -0.25400642               0.6763049  6.700636e-29  9.073332e-25   1.918372e-06   2.597667e-02
```

Next, we generate the prioritization table. This table contains the
rankings of ligand-receptor pairs based on the different criteria. We
provide two scenarios: `case_control` and `one_condition`. In the
“case_control” scenario, all weights are set to 1. If “one_condition”,
the weights are set to 0 for condition specificity, 2 for ligand
activity, and 1 for the remaining criteria. Users can also provide their
own weights using the `prioritizing_weights` argument.

``` r
prior_table <- generate_prioritization_tables(info_tables$sender_receiver_info,
                                              info_tables$sender_receiver_de,
                                              ligand_activities,
                                              info_tables$lr_condition_de,
                                              scenario = "case_control")

prior_table %>% head
## # A tibble: 6 × 51
##   sender receiver ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg p_val_ligand p_adj_ligand p_val_receptor p_adj_receptor
##   <chr>  <chr>    <chr>  <chr>         <dbl>        <dbl>                   <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
## 1 NK     CD8 T    Ptprc  Dpp4          0.596        0.164                   0.380    2.18e-  7    2.96e-  3       0.000663              1
## 2 Mono   CD8 T    Ptprc  Dpp4          0.438        0.164                   0.301    3.52e-  5    4.77e-  1       0.000663              1
## 3 Mono   CD8 T    Cxcl10 Dpp4          4.27         0.164                   2.22     2.53e- 79    3.43e- 75       0.000663              1
## 4 Mono   CD8 T    Cxcl9  Dpp4          4.33         0.164                   2.25     3.83e-124    5.19e-120       0.000663              1
## 5 Treg   CD8 T    Ptprc  Dpp4          0.282        0.164                   0.223    1.44e-  2    1   e+  0       0.000663              1
## 6 Mono   CD8 T    Cxcl11 Dpp4          2.36         0.164                   1.26     9.28e-121    1.26e-116       0.000663              1
## # ℹ 40 more variables: pct_expressed_sender <dbl>, pct_expressed_receiver <dbl>, avg_ligand <dbl>, avg_receptor <dbl>, ligand_receptor_prod <dbl>,
## #   lfc_pval_ligand <dbl>, p_val_ligand_adapted <dbl>, scaled_lfc_ligand <dbl>, scaled_p_val_ligand <dbl>, scaled_lfc_pval_ligand <dbl>,
## #   scaled_p_val_ligand_adapted <dbl>, activity <dbl>, rank <dbl>, activity_zscore <dbl>, scaled_activity <dbl>, lfc_pval_receptor <dbl>,
## #   p_val_receptor_adapted <dbl>, scaled_lfc_receptor <dbl>, scaled_p_val_receptor <dbl>, scaled_lfc_pval_receptor <dbl>,
## #   scaled_p_val_receptor_adapted <dbl>, scaled_avg_exprs_ligand <dbl>, scaled_avg_exprs_receptor <dbl>, lfc_ligand_group <dbl>,
## #   p_val_ligand_group <dbl>, lfc_pval_ligand_group <dbl>, p_val_ligand_adapted_group <dbl>, scaled_lfc_ligand_group <dbl>,
## #   scaled_p_val_ligand_group <dbl>, scaled_lfc_pval_ligand_group <dbl>, scaled_p_val_ligand_adapted_group <dbl>, lfc_receptor_group <dbl>, …
```

As you can see, the resulting table now show the rankings for
*ligand-receptor interactions of a sender-receiver cell type pair*,
instead of just the prioritized ligands. Cxcl10 now went up in the
rankings due to both the high expression of its potential receptor Dpp4
and its high celltype specificity (`scaled_lfc_ligand`). You can also
see this in the visualizations further below.

We included all columns here, but if you just want relevant columns that
were used to calculate the ranking:

``` r
prior_table %>% select(c('sender', 'receiver', 'ligand', 'receptor', 'scaled_p_val_ligand_adapted', 'scaled_p_val_receptor_adapted', 'scaled_avg_exprs_ligand', 'scaled_avg_exprs_receptor', 'scaled_p_val_ligand_adapted_group', 'scaled_p_val_receptor_adapted_group', 'scaled_activity'))
## # A tibble: 1,272 × 11
##    sender receiver ligand receptor scaled_p_val_ligand_…¹ scaled_p_val_recepto…² scaled_avg_exprs_lig…³ scaled_avg_exprs_rec…⁴ scaled_p_val_ligand_…⁵
##    <chr>  <chr>    <chr>  <chr>                     <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>
##  1 NK     CD8 T    Ptprc  Dpp4                      0.869                  0.829                  1.00                   1.00                   0.850
##  2 Mono   CD8 T    Ptprc  Dpp4                      0.841                  0.829                  0.867                  1.00                   0.850
##  3 Mono   CD8 T    Cxcl10 Dpp4                      0.960                  0.829                  1.00                   1.00                   0.929
##  4 Mono   CD8 T    Cxcl9  Dpp4                      0.975                  0.829                  1.00                   1.00                   0.787
##  5 Treg   CD8 T    Ptprc  Dpp4                      0.756                  0.829                  0.741                  1.00                   0.850
##  6 Mono   CD8 T    Cxcl11 Dpp4                      0.973                  0.829                  1.00                   1.00                   0.732
##  7 B      CD8 T    Ptprc  Dpp4                      0.748                  0.829                  0.666                  1.00                   0.850
##  8 DC     CD8 T    Icam1  Il2rg                     0.876                  0.714                  1.00                   0.995                  0.717
##  9 DC     CD8 T    Ccl22  Dpp4                      0.997                  0.829                  1.00                   1.00                   0.539
## 10 NK     CD8 T    Cd320  Jaml                      0.889                  0.943                  0.905                  1.00                   0.472
## # ℹ 1,262 more rows
## # ℹ abbreviated names: ¹​scaled_p_val_ligand_adapted, ²​scaled_p_val_receptor_adapted, ³​scaled_avg_exprs_ligand, ⁴​scaled_avg_exprs_receptor,
## #   ⁵​scaled_p_val_ligand_adapted_group
## # ℹ 2 more variables: scaled_p_val_receptor_adapted_group <dbl>, scaled_activity <dbl>
```

Note that we appended the suffix ’\_group’ to columns that refer to
differential expression between conditions, e.g., `lfc_ligand_group` and
`lfc_receptor_group.`

## Step-by-step prioritization

`generate_info_tables` is a wrapper function that calculates all the
information needed for the prioritization. However, in some cases you
may need more flexibility on how these values are calculated (but note
that you can pass extra arguments to `generate_info_tables` that will
get passed on to `FindMarkers`, `FindAllMarkers`, and
`AverageExpression`). Below, we show how we use helper functions
`calculate_de` and `get_exprs_avg` to calculate the DE and get the
average expression used for cell type specificity. `process_table_to_ic`
transforms these different dataframes so they are compatible with the
`generate_prioritization_tables` function.

``` r
lr_network_renamed <- lr_network_filtered %>% rename(ligand=from, receptor=to)

# Only calculate DE for LCMV condition, with genes that are in the ligand-receptor network
DE_table <- calculate_de(seuratObj, celltype_colname = "celltype",
                         condition_colname = "aggregate", condition_oi = condition_oi,
                         features = unique(unlist(lr_network_filtered)))

# Average expression information - only for LCMV condition
expression_info <- get_exprs_avg(seuratObj, "celltype", condition_colname = "aggregate", condition_oi = condition_oi,
                                 features = unique(unlist(lr_network_filtered)))

# Calculate condition specificity - only for datasets with two conditions!
condition_markers <- FindMarkers(object = seuratObj, ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "aggregate", min.pct = 0, logfc.threshold = 0,
                                 features = unique(unlist(lr_network_filtered))) %>% rownames_to_column("gene")

# Combine DE of senders and receivers -> used for prioritization
processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_renamed,
                                         senders_oi = sender_celltypes, receivers_oi = receiver)
  
processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_renamed)

processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network_renamed)
```

And here is how you can define custom weights:

``` r
prioritizing_weights = c("de_ligand" = 1,
                          "de_receptor" = 1,
                          "activity_scaled" = 1,
                          "exprs_ligand" = 1,
                          "exprs_receptor" = 1,
                         "ligand_condition_specificity" = 1,
                         "receptor_condition_specificity" = 1)
```

``` r
prior_table <- generate_prioritization_tables(processed_expr_table,
                                              processed_DE_table,
                                              ligand_activities,
                                              processed_condition_markers,
                                              prioritizing_weights)

prior_table %>% head
## # A tibble: 6 × 51
##   sender receiver ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg p_val_ligand p_adj_ligand p_val_receptor p_adj_receptor
##   <chr>  <chr>    <chr>  <chr>         <dbl>        <dbl>                   <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
## 1 NK     CD8 T    Ptprc  Dpp4          0.596        0.164                   0.380    2.18e-  7    2.96e-  3       0.000663              1
## 2 Mono   CD8 T    Ptprc  Dpp4          0.438        0.164                   0.301    3.52e-  5    4.77e-  1       0.000663              1
## 3 Mono   CD8 T    Cxcl10 Dpp4          4.27         0.164                   2.22     2.53e- 79    3.43e- 75       0.000663              1
## 4 Mono   CD8 T    Cxcl9  Dpp4          4.33         0.164                   2.25     3.83e-124    5.19e-120       0.000663              1
## 5 Treg   CD8 T    Ptprc  Dpp4          0.282        0.164                   0.223    1.44e-  2    1   e+  0       0.000663              1
## 6 Mono   CD8 T    Cxcl11 Dpp4          2.36         0.164                   1.26     9.28e-121    1.26e-116       0.000663              1
## # ℹ 40 more variables: pct_expressed_sender <dbl>, pct_expressed_receiver <dbl>, avg_ligand <dbl>, avg_receptor <dbl>, ligand_receptor_prod <dbl>,
## #   lfc_pval_ligand <dbl>, p_val_ligand_adapted <dbl>, scaled_lfc_ligand <dbl>, scaled_p_val_ligand <dbl>, scaled_lfc_pval_ligand <dbl>,
## #   scaled_p_val_ligand_adapted <dbl>, activity <dbl>, rank <dbl>, activity_zscore <dbl>, scaled_activity <dbl>, lfc_pval_receptor <dbl>,
## #   p_val_receptor_adapted <dbl>, scaled_lfc_receptor <dbl>, scaled_p_val_receptor <dbl>, scaled_lfc_pval_receptor <dbl>,
## #   scaled_p_val_receptor_adapted <dbl>, scaled_avg_exprs_ligand <dbl>, scaled_avg_exprs_receptor <dbl>, lfc_ligand_group <dbl>,
## #   p_val_ligand_group <dbl>, lfc_pval_ligand_group <dbl>, p_val_ligand_adapted_group <dbl>, scaled_lfc_ligand_group <dbl>,
## #   scaled_p_val_ligand_group <dbl>, scaled_lfc_pval_ligand_group <dbl>, scaled_p_val_ligand_adapted_group <dbl>, lfc_receptor_group <dbl>, …
```

# Prioritizing across multiple receivers

As NicheNet is a receiver-based pipeline, to prioritize ligand-receptor
pairs across multiple receivers, we need to perform the NicheNet
analysis for each receiver separately. Let’s suppose we want to
prioritize ligand-receptor pairs across all T cells (CD4, CD8, and
Tregs). The CD8 T analysis has already been performed above. We will use
the wrapper function to perform a basic NicheNet analysis on the other
two:

``` r
nichenet_output <- lapply(c("CD4 T", "Treg"), function(receiver_ct){
  nichenet_seuratobj_aggregate(receiver = receiver_ct,
                             seurat_obj = seuratObj,
                             condition_colname = "aggregate",
                             condition_oi = condition_oi,
                             condition_reference = condition_reference,
                             sender = sender_celltypes,
                             ligand_target_matrix = ligand_target_matrix,
                             lr_network = lr_network,
                             weighted_networks = weighted_networks,
                             expression_pct = 0.05)
  
}) %>% setNames(c("CD4 T", "Treg"))
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"
```

To generate the dataframes used for prioritization, we will simply
change the `lr_network_filtered` argument to only calculate DE and
expression values for ligand-receptor pairs of interest.

``` r
info_tables2 <- lapply(names(nichenet_output), function(receiver_ct) {
  generate_info_tables(seuratObj,
                       celltype_colname = "celltype",
                       senders_oi = sender_celltypes,
                        receivers_oi = receiver_ct,
                          lr_network_filtered = lr_network %>%
                            filter(from %in% nichenet_output[[receiver_ct]]$ligand_activities$test_ligand &
                                     to %in% nichenet_output[[receiver_ct]]$background_expressed_genes),
                          condition_colname = "aggregate",
                          condition_oi = condition_oi,
                          condition_reference = condition_reference,
                          scenario = "case_control")
})
```

We can then combine the results from `generate_info_tables` using
`bind_rows`, which will concatenate the rows together. For the ligand
activities, we will also add an additional column containing the
receiver cell type. Note that for the average expression table
(`sender_receiver_info`) and condition specificity (`lr_condition_de`),
we need to remove duplicate rows.

``` r
# Add CD8 T to list
info_tables2[[3]] <- info_tables

# bind rows of each element of info_tables using pmap
info_tables_combined <- purrr::pmap(info_tables2, bind_rows)

# Combine ligand activities and add receiver information
ligand_activities_combined <- bind_rows(nichenet_output$`CD4 T`$ligand_activities %>% mutate(receiver = "CD4 T"),
                                        nichenet_output$Treg$ligand_activities %>% mutate(receiver = "Treg"),
                                        ligand_activities %>% mutate(receiver = "CD8 T"))

prior_table_combined <- generate_prioritization_tables(
  sender_receiver_info = info_tables_combined$sender_receiver_info %>% distinct,
  sender_receiver_de = info_tables_combined$sender_receiver_de,
  ligand_activities = ligand_activities_combined,
  lr_condition_de = info_tables_combined$lr_condition_de %>% distinct,
  scenario = "case_control")

head(prior_table_combined)
## # A tibble: 6 × 51
##   sender receiver ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg p_val_ligand p_adj_ligand p_val_receptor p_adj_receptor
##   <chr>  <chr>    <chr>  <chr>         <dbl>        <dbl>                   <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
## 1 NK     CD8 T    Ptprc  Dpp4          0.596        0.164                   0.380  0.000000218      0.00296       6.63e- 4       1   e+ 0
## 2 NK     CD4 T    Ptprc  Cd4           0.596        0.996                   0.796  0.000000218      0.00296       2.63e-34       3.56e-30
## 3 B      CD4 T    H2-Eb1 Cd4           4.02         0.996                   2.51   0                0             2.63e-34       3.56e-30
## 4 Mono   CD8 T    Ptprc  Dpp4          0.438        0.164                   0.301  0.0000352        0.477         6.63e- 4       1   e+ 0
## 5 Mono   CD4 T    Ptprc  Cd4           0.438        0.996                   0.717  0.0000352        0.477         2.63e-34       3.56e-30
## 6 NK     CD4 T    Ptprc  Cd247         0.596        0.457                   0.526  0.000000218      0.00296       5.61e- 4       1   e+ 0
## # ℹ 40 more variables: pct_expressed_sender <dbl>, pct_expressed_receiver <dbl>, avg_ligand <dbl>, avg_receptor <dbl>, ligand_receptor_prod <dbl>,
## #   lfc_pval_ligand <dbl>, p_val_ligand_adapted <dbl>, scaled_lfc_ligand <dbl>, scaled_p_val_ligand <dbl>, scaled_lfc_pval_ligand <dbl>,
## #   scaled_p_val_ligand_adapted <dbl>, activity <dbl>, rank <dbl>, activity_zscore <dbl>, scaled_activity <dbl>, lfc_pval_receptor <dbl>,
## #   p_val_receptor_adapted <dbl>, scaled_lfc_receptor <dbl>, scaled_p_val_receptor <dbl>, scaled_lfc_pval_receptor <dbl>,
## #   scaled_p_val_receptor_adapted <dbl>, scaled_avg_exprs_ligand <dbl>, scaled_avg_exprs_receptor <dbl>, lfc_ligand_group <dbl>,
## #   p_val_ligand_group <dbl>, lfc_pval_ligand_group <dbl>, p_val_ligand_adapted_group <dbl>, scaled_lfc_ligand_group <dbl>,
## #   scaled_p_val_ligand_group <dbl>, scaled_lfc_pval_ligand_group <dbl>, scaled_p_val_ligand_adapted_group <dbl>, lfc_receptor_group <dbl>, …
```

### Extra visualization of ligand-receptor pairs

In addition to the usual heatmap visualizations, we provide a function
`make_circos_lr` to visualize the ligand-receptor pairs in a circos
plot. This was originally written for the (now deprecated) Differential
NicheNet vignettes. The function takes in a prioritization table and a
named vector for the color of senders and receivers. We first specify
the number of top ligand-receptor pairs to show with `n`.

``` r
# Get top n ligand-receptor pairs
prior_table_oi <- prior_table_combined %>% slice_max(prioritization_score, n = 50)

# Define colors for senders and receivers
senders_receivers <- prior_table_oi %>% select(sender, receiver) %>% unlist %>% unique %>% sort
celltype_colors <- RColorBrewer::brewer.pal(length(senders_receivers), name = 'Set3') %>%
  magrittr::set_names(senders_receivers)

circos_plot <- make_circos_lr(prior_table_oi,
               colors_sender = celltype_colors, colors_receiver = celltype_colors)
```

``` r
circos_plot
```

![](seurat_steps_prioritization_files/figure-gfm/lr-circos-1.png)<!-- -->

Furthermore, we provide the function `make_mushroom_plot` which allows
you to display expression of ligand-receptor pairs in a specific
receiver. By default, the fill gradient shows the LFC between cell
types, while the size of the semicircle corresponds to the scaled mean
expression. You can also choose to show the rankings of each
ligand-receptor-sender pair with `show_rankings`, as well as show all
data points for context (`show_all_datapoints`).
`true_color_range = TRUE` will adjust the limits of the color gradient
to the min-max of the values, instead of the limit being from 0 to 1.
Note that the numbers displayed here are the rankings within the chosen
cell type and not across all receiver cell types (in case of multiple
receivers).

``` r
receiver_oi <- "CD8 T"
legend_adjust <- c(0.7, 0.7)
make_mushroom_plot(prior_table_combined %>% filter(receiver == receiver_oi),
                   top_n = 30, 
                   true_color_range = TRUE,
                   show_rankings = TRUE,
                   show_all_datapoints = TRUE) +
  theme(legend.justification = legend_adjust,
        axis.title.x = element_text(hjust = 0.25))
```

![](seurat_steps_prioritization_files/figure-gfm/mushroom-plot-1-1.png)<!-- -->

Furthermore, you can change the “size” and “fill” values to certain
columns from the prioritization table (those with the `_ligand` or
`_receptor` suffix).

``` r
print(paste0("Column names that you can use are: ", paste0(prior_table %>% select(ends_with(c("_ligand", "_receptor", "_sender", "_receiver"))) %>% colnames() %>%
  str_remove("_ligand|_receptor|_sender|_receiver") %>% unique, collapse = ", ")))
## [1] "Column names that you can use are: lfc, p_val, p_adj, avg, lfc_pval, scaled_lfc, scaled_p_val, scaled_lfc_pval, scaled_avg_exprs, pct_expressed"

# Change size and color columns
make_mushroom_plot(prior_table, top_n = 30, size = "pct_expressed", color = "scaled_avg_exprs") +
  theme(legend.justification = legend_adjust,
        axis.title.x = element_text(hjust = 0.25))
```

![](seurat_steps_prioritization_files/figure-gfm/mushroom-plot-2-1.png)<!-- -->

``` r
sessionInfo()
## R version 4.3.2 (2023-10-31)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS Stream 8
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
##  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Asia/Bangkok
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] forcats_1.0.0      stringr_1.5.0      dplyr_1.1.4        purrr_1.0.2        readr_2.1.2        tidyr_1.3.0        tibble_3.2.1      
##  [8] ggplot2_3.4.4      tidyverse_1.3.1    SeuratObject_5.0.1 Seurat_4.4.0       nichenetr_2.0.5   
## 
## loaded via a namespace (and not attached):
##   [1] fs_1.6.3               matrixStats_1.2.0      spatstat.sparse_3.0-3  bitops_1.0-7           lubridate_1.9.3        httr_1.4.7            
##   [7] RColorBrewer_1.1-3     doParallel_1.0.17      tools_4.3.2            sctransform_0.4.0      backports_1.4.1        utf8_1.2.4            
##  [13] R6_2.5.1               lazyeval_0.2.2         uwot_0.1.16            GetoptLong_1.0.5       withr_2.5.2            sp_2.1-2              
##  [19] gridExtra_2.3          fdrtool_1.2.17         progressr_0.14.0       cli_3.6.2              spatstat.explore_3.2-1 labeling_0.4.3        
##  [25] spatstat.data_3.0-3    randomForest_4.7-1.1   proxy_0.4-27           ggridges_0.5.5         pbapply_1.7-2          foreign_0.8-85        
##  [31] parallelly_1.36.0      limma_3.56.2           readxl_1.4.3           rstudioapi_0.15.0      gridGraphics_0.5-1     visNetwork_2.1.2      
##  [37] generics_0.1.3         shape_1.4.6            ica_1.0-3              spatstat.random_3.2-2  car_3.1-2              Matrix_1.6-4          
##  [43] fansi_1.0.6            S4Vectors_0.38.1       abind_1.4-5            lifecycle_1.0.4        yaml_2.3.8             carData_3.0-5         
##  [49] recipes_1.0.7          Rtsne_0.17             grid_4.3.2             promises_1.2.1         crayon_1.5.2           miniUI_0.1.1.1        
##  [55] lattice_0.21-9         haven_2.4.3            cowplot_1.1.2          pillar_1.9.0           knitr_1.45             ComplexHeatmap_2.16.0 
##  [61] rjson_0.2.21           future.apply_1.11.0    codetools_0.2-19       leiden_0.3.9           glue_1.6.2             data.table_1.14.10    
##  [67] vctrs_0.6.5            png_0.1-8              spam_2.10-0            cellranger_1.1.0       gtable_0.3.4           assertthat_0.2.1      
##  [73] gower_1.0.1            xfun_0.41              mime_0.12              prodlim_2023.08.28     survival_3.5-7         timeDate_4032.109     
##  [79] iterators_1.0.14       hardhat_1.3.0          lava_1.7.3             DiagrammeR_1.0.10      ellipsis_0.3.2         fitdistrplus_1.1-11   
##  [85] ROCR_1.0-11            ipred_0.9-14           nlme_3.1-163           RcppAnnoy_0.0.21       irlba_2.3.5.1          KernSmooth_2.23-22    
##  [91] rpart_4.1.21           colorspace_2.1-0       BiocGenerics_0.46.0    DBI_1.1.3              Hmisc_5.1-0            nnet_7.3-19           
##  [97] tidyselect_1.2.0       compiler_4.3.2         rvest_1.0.2            htmlTable_2.4.1        xml2_1.3.6             plotly_4.10.0         
## [103] shadowtext_0.1.2       checkmate_2.3.1        scales_1.3.0           caTools_1.18.2         lmtest_0.9-40          digest_0.6.33         
## [109] goftest_1.2-3          spatstat.utils_3.0-4   rmarkdown_2.11         htmltools_0.5.7        pkgconfig_2.0.3        base64enc_0.1-3       
## [115] highr_0.10             dbplyr_2.1.1           fastmap_1.1.1          rlang_1.1.2            GlobalOptions_0.1.2    htmlwidgets_1.6.2     
## [121] shiny_1.7.1            farver_2.1.1           zoo_1.8-12             jsonlite_1.8.8         ModelMetrics_1.2.2.2   magrittr_2.0.3        
## [127] Formula_1.2-5          dotCall64_1.1-1        patchwork_1.1.3        munsell_0.5.0          Rcpp_1.0.11            ggnewscale_0.4.9      
## [133] reticulate_1.34.0      stringi_1.7.6          pROC_1.18.5            MASS_7.3-60            plyr_1.8.9             parallel_4.3.2        
## [139] listenv_0.9.0          ggrepel_0.9.4          deldir_2.0-2           splines_4.3.2          tensor_1.5             hms_1.1.3             
## [145] circlize_0.4.15        igraph_1.2.11          ggpubr_0.6.0           spatstat.geom_3.2-7    ggsignif_0.6.4         reshape2_1.4.4        
## [151] stats4_4.3.2           reprex_2.0.1           evaluate_0.23          modelr_0.1.8           tzdb_0.4.0             foreach_1.5.2         
## [157] tweenr_2.0.2           httpuv_1.6.13          RANN_2.6.1             polyclip_1.10-6        future_1.33.0          clue_0.3-64           
## [163] scattermore_1.2        ggforce_0.4.1          broom_0.7.12           xtable_1.8-4           e1071_1.7-14           rstatix_0.7.2         
## [169] later_1.3.2            viridisLite_0.4.2      class_7.3-22           IRanges_2.34.1         cluster_2.1.4          timechange_0.2.0      
## [175] globals_0.16.2         caret_6.0-94
```

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-medaglia_spatial_2017" class="csl-entry">

Medaglia, Chiara, Amir Giladi, Liat Stoler-Barak, Marco De Giovanni,
Tomer Meir Salame, Adi Biram, Eyal David, et al. 2017. “Spatial
Reconstruction of Immune Niches by Combining Photoactivatable Reporters
and <span class="nocase">scRNA</span>-Seq.” *Science*, December,
eaao4277. <https://doi.org/10.1126/science.aao4277>.

</div>

</div>
