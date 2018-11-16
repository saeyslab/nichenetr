Single-cell NicheNet's ligand activity analysis
================
Robin Browaeys
2018-11-12

This vignette shows how NicheNet can be used to predict which ligands might be active in single-cells. If a ligand has a high activity in a cell, this means that target genes of that ligand are stronger expressed in that cell than in other cells. In this example, we will use data from Puram et al. to explore intercellular communication in the tumor microenvironment in head and neck squamous cell carcinoma (HNSCC) (See Puram et al. 2017). More specifically, we will assess the activity of fibroblast ligands in malignant cells. The used ligand-target matrix and example expression data of interacting cells can be downloaded from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484138.svg)](https://doi.org/10.5281/zenodo.1484138).

### Load nichenetr and tidyverse

``` r
library(nichenetr)
library(tidyverse)
```

### Read in expression data of interacting cells

First, we will read in the single-cell data from fibroblast and malignant cells from HNSCC tumors (See Puram et al. 2017).

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/1484138/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

Secondly, we will determine which genes are expressed in fibroblasts and malignant cells from high quality primary tumors. Therefore, we wil not consider cells from tumor samples of less quality or from lymph node metastases. To determine expressed genes, we use the definition used by of Puram et al.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

fibroblast_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter((tumor %in% tumors_remove == FALSE)) %>% filter(`non-cancer cell type` == "Fibroblast") %>% .$cell
malignant_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1) %>% filter((tumor %in% tumors_remove == FALSE)) %>% .$cell

expressed_genes_fibroblasts = expression[fibroblast_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
```

### Load the ligand-target model we want to use

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/1484138/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     0.0007737318 0.0006554957 6.209527e-04 0.0005819125 6.004103e-04
## A1BG-AS1 0.0003130167 0.0002617725 2.637518e-04 0.0002528345 2.670414e-04
## A1CF     0.0007869330 0.0006304809 5.715833e-04 0.0005322643 5.775608e-04
## A2M      0.0013779875 0.0011799983 1.118986e-03 0.0010957258 1.145126e-03
## A2M-AS1  0.0001141186 0.0001077393 9.434595e-05 0.0000942375 9.862858e-05
```

### Perform NicheNet's single-cell ligand activity analysis

In a first step, we will define a set of potentially active ligands. As potentially active ligands, we will use ligands that are 1) expressed by fibroblasts and 2) can bind a (putative) receptor expressed by malignant cells. Putative ligand-receptor links were gathered from NicheNet's ligand-receptor data sources.

``` r
lr_network = readRDS(url("https://zenodo.org/record/1484138/files/lr_network.rds"))
ligands = lr_network$from %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_fibroblasts)
receptors = lr_network$to %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% .$from %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"
```

In a second step, we will scale the single-cell expression data (including only expressed genes).

``` r
apply_quantile_scale = function (x, addend, multiplier)
    # same function as apply_quantile_scale from dynutils (copied here for use in vignette to avoid having dynutils as dependency)
  # credits to the great (w/z)outer and r(obrecht)cannood(t) from dynverse (https://github.com/dynverse)! 
{
    if (is.null(dim(x))) {
        sc <- apply_quantile_scale(matrix(x, ncol = 1), addend = addend, 
            multiplier = multiplier)
        out <- sc[, 1]
        names(out) <- names(x)
        attr(out, "addend") <- attr(sc, "addend")
        attr(out, "multiplier") <- attr(sc, "multiplier")
        out
    }
    else {
        y <- apply_uniform_scale(x, addend, multiplier)
        y[y > 1] <- 1
        y[y < 0] <- 0
        y
    }
}
apply_uniform_scale = function (x, addend, multiplier) 
    # same function as apply_uniform_scale from dynutils (copied here for use in vignette to avoid having dynutils as dependency)
  # credits to the great (w/z)outer and r(obrecht)cannood(t) from dynverse (https://github.com/dynverse)! 
{
    if (is.null(dim(x))) {
        sc <- apply_uniform_scale(matrix(x, ncol = 1), addend = addend,multiplier = multiplier)
        out <- sc[, 1]
        names(out) <- names(x)
        attr(out, "addend") <- attr(sc, "addend")
        attr(out, "multiplier") <- attr(sc, "multiplier")
        out
    }
    else {
        y <- x %>% sweep(2, addend, "+") %>% sweep(2, multiplier, 
            "*")
        attr(y, "addend") <- addend
        attr(y, "multiplier") <- multiplier
        y
    }
}
scale_quantile = function (x, outlier_cutoff = 0.05) 
  # same function as scale_quantile from dynutils (copied here for use in vignette to avoid having dynutils as dependency)
  # credits to the great (w/z)outer and r(obrecht)cannood(t) from dynverse (https://github.com/dynverse)! 
{
    if (is.null(dim(x))) {
        sc <- scale_quantile(matrix(x, ncol = 1), outlier_cutoff = outlier_cutoff)
        out <- sc[, 1]
        names(out) <- names(x)
        attr(out, "addend") <- attr(sc, "addend")
        attr(out, "multiplier") <- attr(sc, "multiplier")
        out
    }
    else {
        quants <- apply(x, 2, stats::quantile, c(outlier_cutoff, 
            1 - outlier_cutoff), na.rm = TRUE)
        addend <- -quants[1, ]
        divisor <- apply(quants, 2, diff)
        divisor[divisor == 0] <- 1
        apply_quantile_scale(x, addend, 1/divisor)
    }
}
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
expression_scaled = expression %>% .[malignant_ids,background_expressed_genes] %>% scale_quantile()
```

Now perform the ligand activity analysis: infer how well NicheNet's ligand-target potential scores can predict whether a gene belongs to most strongly expressed genes in a cell compared to other cells. To reduce the running time for this vignette, we will perform the analysis only on 10 example cells from the HN5 tumor. In practice, ligand activity analysis for several cells can be better run in parallel!

``` r
convert_single_cell_expression_to_settings = function(cell_id, expression_matrix, setting_name, setting_from, regression = FALSE){
  # input check
  requireNamespace("dplyr")
  
  if (regression == TRUE){
    response = expression_matrix[cell_id,]
  } else {
    response_continuous = expression_matrix[cell_id,]
    response = response_continuous >= quantile(response_continuous,0.975)
  }
  return(list(name = paste0(setting_name,"_",cell_id), from = setting_from, response = response))
}

predict_single_cell_ligand_activities = function(cell_ids, expression_scaled,ligand_target_matrix, potential_ligands, single = TRUE,...){
  settings_single_cell_ligand_pred = cell_ids %>% lapply(convert_single_cell_expression_to_settings, expression_scaled, "", potential_ligands)
  if (single == TRUE){
    settings_ligand_prediction = settings_single_cell_ligand_pred %>% convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = TRUE)
    
    ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE) %>% bind_rows() %>% mutate(setting = gsub("^_","",setting))
      
  } else {
    settings_ligand_prediction = settings_single_cell_ligand_pred %>% convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = FALSE)
        
    ligand_importances = settings_ligand_prediction %>% lapply(get_multi_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE, ...) %>% bind_rows() %>% mutate(setting = gsub("^_","",setting))
      
  }
  return(ligand_importances)
}

malignant_hn18_ids = sample_info %>% filter(tumor == "HN5") %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1)  %>% .$cell %>% head(10)

ligand_importances = predict_single_cell_ligand_activities(cell_ids = malignant_hn18_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

### Ligand prioritization by regression analysis

Furthermore, we will also show how you can perform additional analyses by linking the ligand activity in cells to other properties of cells in order to prioritize ligands. As toy example, we will score malignant cells here on the extent to which they express the core p-EMT gene "TGFBI".

``` r
# pemt_scores_tbl = tibble(cell = names(scores_pemt), pemt_score = scores_pemt)
```

Then, we will determine the correlation between these p-EMT scores and ligand activities over all cells to prioritize p-EMT-inducing ligands. We hypothesize that ligands might be potential regulators of the p-EMT program if higher ligand activities are associated with higher p-EMT scores. Based on this correlation, we obtained a ranking of potential p-EMT-inducing ligands.

``` r
# single_ligand_importances = ligand_importances_cell_analysis_single %>% select(setting, test_ligand, pearson) %>% filter(test_ligand %in% expressed_ligands)
# 
# single_ligand_importances_pearson_norm = single_ligand_importances %>% group_by(setting) %>% mutate(pearson = nichenetr::scaling_modified_zscore(pearson)) %>% ungroup()
# single_ligand_importances_pearson_norm = single_ligand_importances_pearson_norm %>% rename(cell = setting, ligand = test_ligand) 
# single_ligand_importances_pearson_norm = single_ligand_importances_pearson_norm %>% distinct(cell,ligand,pearson) %>% spread(cell, pearson,fill = min(.$pearson)) 
# single_ligand_importances_pearson_norm_matrix = single_ligand_importances_pearson_norm %>% select(-ligand) %>% t() %>% magrittr::set_colnames(single_ligand_importances_pearson_norm$ligand)
# rownames(single_ligand_importances_pearson_norm_matrix) = single_ligand_importances_pearson_norm_matrix %>% rownames() 
# single_importances_norm_pearson_df = single_ligand_importances_pearson_norm_matrix %>% data.frame() %>% rownames_to_column("cell") %>% tbl_df()
# 
# # combine ligand activities and pemt scores
# regr_pearson_pemt_df = single_ligand_activity_score_regression(single_importances_norm_pearson_df,pemt_scores_tbl %>% rename(score = pemt_score))
# regr_pearson_pemt_df %>% arrange(-pearson_regression) %>% select(r_squared, pearson_regression, ligand)
```

### References

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel, Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017. “Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24. doi:[10.1016/j.cell.2017.10.044](https://doi.org/10.1016/j.cell.2017.10.044).
