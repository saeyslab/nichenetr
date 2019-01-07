This vignette shows how NicheNet can be used to predict which ligands might regulate a given set of genes. In this example, we will use data from Puram et al. to explore intercellular communication in the tumor microenvironment in head and neck squamous cell carcinoma (HNSCC) (See Puram et al. 2017). More specifically, we will look at which ligands expressed by cancer-associated fibroblast can induce a specific gene program in neighboring malignant cells. This program, a partial epithelial-mesenschymal transition (p-EMT) program, could be linked by Puram et al. to metastasis. For this analysis, we will assess the ligand activity of each ligand, or in other words, we will assess how well each fibroblast ligand can predict the p-EMT gene set compared to the background of expressed genes. This allows us to prioritize p-EMT-regulating ligands. In a final step, we will infer target genes of these top ligands and signaling paths between these targets and ligands. The used ligand-target matrix and example expression data of interacting cells can be downloaded from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484138.svg)](https://doi.org/10.5281/zenodo.1484138).

### Load nichenetr and tidyverse

``` r
library(nichenetr)
library(dplyr)
library(ggplot2)
```

### Read in expression data of interacting cells

First, we will read in the single-cell data from fibroblast and malignant cells from HNSCC tumors.

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

### Load the gene set of interest and background of genes

Because we want to investigate how fibroblast regulate the expression of p-EMT genes in malignant cells, we will use the p-EMT gene set defined by Puram et al. as gene set of interset and use all genes expressed in malignant cells as background of genes.

``` r
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/1484138/files/pemt_signature.txt"), col_names = "gene") %>% .$gene %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(pemt_geneset)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
```

### Perform NicheNet's ligand activity analysis on the gene set of interest

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

Now perform the ligand activity analysis: infer how well NicheNet's ligand-target potential scores can predict whether a gene belongs to the p-EMT program or not.

``` r
predict_ligand_activities = function(geneset,background_expressed_genes,ligand_target_matrix, potential_ligands, single = TRUE,...){
  setting = list(geneset) %>% 
    lapply(convert_gene_list_settings_evaluation, name = "gene set", ligands_oi = potential_ligands, background = background_expressed_genes)
  if (single == TRUE){
    settings_ligand_prediction = setting %>% 
      convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = TRUE)
    ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE) %>% bind_rows()
      
  } else {
    settings_ligand_prediction = setting %>% 
      convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = FALSE)
    ligand_importances = settings_ligand_prediction %>% lapply(get_multi_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE, ...) %>% bind_rows()
      
  }
  return(ligand_importances)
}

ligand_importances = predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

Now, we want to rank the ligands based on their ligand activity. In our study, we showed that the pearson correlation between a ligand's target predictions and the observed transcriptional response was the most informative measure for ligand activity. Therefore, we will rank the ligands based on their pearson correlation coefficient.

``` r
ligand_importances %>% arrange(-pearson) %>% select(test_ligand,auroc,aupr,pearson)
## # A tibble: 134 x 4
##    test_ligand auroc   aupr pearson
##    <chr>       <dbl>  <dbl>   <dbl>
##  1 AGT         0.676 0.0624   0.122
##  2 CXCL12      0.676 0.0549   0.121
##  3 TGFB3       0.696 0.0505   0.119
##  4 IL6         0.681 0.0502   0.115
##  5 CTGF        0.690 0.0502   0.114
##  6 INHBA       0.698 0.0523   0.113
##  7 TNC         0.691 0.0469   0.108
##  8 PTHLH       0.655 0.0511   0.108
##  9 ADAM17      0.663 0.0454   0.105
## 10 BMP5        0.713 0.0423   0.104
## # ... with 124 more rows
best_upstream_ligands = ligand_importances %>% top_n(20, pearson) %>% arrange(-pearson) %>% .$test_ligand
```

We see here that the top-ranked ligands can predict the p-EMT genes reasonably, this implies that ranking of the ligands might be accurate as shown in our study. However, it is possible that for some gene sets, the target gene prediction performance of the top-ranked ligands would not be much better than random prediction. In that case, prioritization of ligands will probably not be accurate.

### Infer target genes of top-ranked ligands

Now we will show how you can look at the regulatory potential scores between ligands and target genes of interest. In this case, we will look at links between top-ranked p-EMT regulating ligands and p-EMT genes.

``` r
infer_ligand_target_links = function(ligands_oi, targets_oi, background_expressed_genes, ligand_target_matrix, k = 2){
  ligand_target_matrix_oi = ligand_target_matrix[background_expressed_genes,ligands_oi]
  
  cutoff = k*length(targets_oi)/length(background_expressed_genes) # consider the top "k times the fraction of positive outcomes" ("true target genes") as postive predictions
  ligand_target_matrix_discrete = make_discrete_ligand_target_matrix(ligand_target_matrix = ligand_target_matrix_oi, error_rate = cutoff, cutoff_method = "quantile")
  ligand_target_matrix_oi[!ligand_target_matrix_discrete] = 0
  
  ligand_target_vis = ligand_target_matrix_oi[targets_oi,ligands_oi]
  ligand_target_vis_filtered = ligand_target_vis[ligand_target_vis %>% apply(1,sum) > 0,ligand_target_vis %>% apply(2,sum) > 0]
  
  distoi = dist(1-cor(t(ligand_target_vis_filtered)))
  hclust_obj = hclust(distoi, method = "ward.D2")
  order_targets = hclust_obj$labels[hclust_obj$order]
  
  distoi_targets = dist(1-cor(ligand_target_vis_filtered))
  hclust_obj = hclust(distoi_targets, method = "ward.D2")
  order_ligands = hclust_obj$labels[hclust_obj$order]
  
  vis_ligand_target_network = ligand_target_vis_filtered[order_targets,order_ligands]
  
}

active_ligand_target_links = infer_ligand_target_links(ligands_oi = best_upstream_ligands, targets_oi = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix)
```

The putatively active ligand-target links will be visualized in a heatmap.

``` r
p_ligand_target_network = active_ligand_target_links %>% t() %>% make_heatmap_ggplot("Ligand","Target", color = "purple",legend_position = "top", x_axis_position = "top") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01))

p_ligand_target_network
```

![](C:/Users/rbrowaey/AppData/Local/Temp/RtmpsP37uZ/preview-1ca4334d53e0.dir/ligand_activity_geneset_files/figure-markdown_github/unnamed-chunk-10-1.png)

### Infer signaling paths beween ligand(s) and target(s) of interest

As follow-up analysis, you can infer possible signaling paths between ligands and targets of interest. You can read how to do this in the following vignette [Inferring ligand-to-target signaling paths](vignettes/ligand_target_signaling_path.md): `vignette("ligand_target_signaling_path", package="nichenetr")`.

### References

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel, Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017. “Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24. doi:[10.1016/j.cell.2017.10.044](https://doi.org/10.1016/j.cell.2017.10.044).
