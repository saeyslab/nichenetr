NicheNet's ligand activity analysis on a gene set of interest
================
Robin Browaeys
2019-01-17

<!-- github markdown built using 
rmarkdown::render("vignettes/ligand_activity_geneset.Rmd", output_format = "github_document")
-->
This vignette shows how NicheNet can be used to predict which ligands might regulate a given set of genes. For this analysis, you need to define:

-   a set of genes of which expression in a "receiver cell" is possibly affected by extracellular protein signals (ligands) (e.g. genes differentially expressed upon cell-cell interaction )
-   a set of potentially active ligands (e.g. ligands expressed by interacting "sender cells")

Therefore, you often first need to process expression data of interacting cells to define both.

In this example, we will use data from Puram et al. to explore intercellular communication in the tumor microenvironment in head and neck squamous cell carcinoma (HNSCC) (See Puram et al. 2017). More specifically, we will look at which ligands expressed by fibroblasts can induce a specific gene program in neighboring malignant cells. This program, a partial epithelial-mesenschymal transition (p-EMT) program, could be linked by Puram et al. to metastasis.

For this analysis, we will assess the ligand activity of each ligand, or in other words, we will assess how well each fibroblast ligand can predict the p-EMT gene set compared to the background of expressed genes. This allows us to prioritize p-EMT-regulating ligands. In a final step, we will then infer target genes of these top ligands.

The used ligand-target matrix and example expression data of interacting cells can be downloaded from Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484138.svg)](https://doi.org/10.5281/zenodo.1484138)

### Load packages required for this vignette

``` r
library(nichenetr)
library(dplyr)
library(ggplot2)
```

### Read in expression data of interacting cells

First, we will read in the publicly available single-cell data from fibroblast and malignant cells from HNSCC tumors.

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/1484138/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

Secondly, we will determine which genes are expressed in fibroblasts and malignant cells from high quality primary tumors. Therefore, we wil not consider cells from tumor samples of less quality or from lymph node metastases. To determine expressed genes, we use the definition used by of Puram et al.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

fibroblast_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "Fibroblast") %>% pull(cell)
malignant_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)

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

As gene set of interest, we consider the genes of which the expression is possibly affected due to communication with other cells.

Because we here want to investigate how fibroblast regulate the expression of p-EMT genes in malignant cells, we will use the p-EMT gene set defined by Puram et al. as gene set of interset and use all genes expressed in malignant cells as background of genes.

``` r
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/1484138/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(pemt_geneset)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"

background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"
```

### Perform NicheNet's ligand activity analysis on the gene set of interest

In a first step, we will define a set of potentially active ligands. As potentially active ligands, we will use ligands that are 1) expressed by fibroblasts and 2) can bind a (putative) receptor expressed by malignant cells. Putative ligand-receptor links were gathered from NicheNet's ligand-receptor data sources.

``` r
lr_network = readRDS(url("https://zenodo.org/record/1484138/files/lr_network.rds"))

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_fibroblasts)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"
```

Now perform the ligand activity analysis: infer how well NicheNet's ligand-target potential scores can predict whether a gene belongs to the p-EMT program or not.

``` r
ligand_activities = predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

Now, we want to rank the ligands based on their ligand activity. In our validation study, we showed that the pearson correlation between a ligand's target predictions and the observed transcriptional response was the most informative measure to define ligand activity. Therefore, we will rank the ligands based on their pearson correlation coefficient.

``` r
ligand_activities %>% arrange(-pearson) 
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
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)
## [1] "AGT"    "CXCL12" "TGFB3"  "IL6"    "CTGF"   "INHBA"
```

We see here that the top-ranked ligands can predict the p-EMT genes reasonably, this implies that ranking of the ligands might be accurate as shown in our study. However, it is possible that for some gene sets, the target gene prediction performance of the top-ranked ligands would not be much better than random prediction. In that case, prioritization of ligands will be less trustworthy.

### Infer target genes of top-ranked ligands and visualize in a heatmap

Now we will show how you can look at the regulatory potential scores between ligands and target genes of interest. In this case, we will look at links between top-ranked p-EMT regulating ligands and p-EMT genes. In the ligand-target heatmaps, we showed regulatory potential scores for interactions between the 20 top-ranked ligands and following target genes: genes that belong to the gene set of interest and to the 250 most strongly predicted targets of at least one of the 20 top-ranked ligands. For visualization purposes, we adapted the ligand-target regulatory potential matrix as follows. Regulatory potential scores were set as 0 if their score was below a predefined threshold, which was here the 0.25 quantile of scores of interactions between the 20 top-ranked ligands and each of their respective 250 top targets.

``` r
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = pemt_geneset, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
## [1] "AGT"
##  [1] "SERPINE1" "MMP2"     "F3"       "PTHLH"    "TIMP3"    "PLAU"    
##  [7] "INHBA"    "COL1A1"   "COL4A2"   "VIM"     
## [1] "CXCL12"
##  [1] "MT2A"     "SERPINE1" "MMP2"     "MMP10"    "MMP1"     "LTBP1"   
##  [7] "PLAU"     "ITGA5"    "VIM"      "CAV1"     "FHL2"     "GJA1"    
## [1] "TGFB3"
## [1] "TNC"      "SERPINE1" "LTBP1"    "VIM"      "PTHLH"    "PLAU"    
## [7] "LAMA3"    "LAMC2"   
## [1] "IL6"
##  [1] "IGFBP3"   "ITGA5"    "MMP1"     "COL1A1"   "MT2A"     "PTHLH"   
##  [7] "SLC39A14" "MMP2"     "C1S"      "TIMP3"    "ITGA6"    "MMP10"   
## [13] "PDLIM7"  
## [1] "CTGF"
##  [1] "MMP2"     "SERPINE1" "MMP1"     "TNC"      "PLAU"     "LTBP1"   
##  [7] "CAV1"     "VIM"      "PTHLH"    "ITGA5"   
## [1] "INHBA"
## [1] "SERPINE1" "LTBP1"    "TNC"      "VIM"      "LAMC2"    "PTHLH"   
## [1] "TNC"
##  [1] "MMP2"     "TIMP3"    "SERPINE1" "MMP1"     "LTBP1"    "VIM"     
##  [7] "CAV1"     "PLAU"     "PTHLH"    "GJA1"    
## [1] "PTHLH"
##  [1] "MMP1"     "PLAU"     "COL1A1"   "MMP2"     "P4HA2"    "SERPINE1"
##  [7] "LTBP1"    "ITGA5"    "VIM"      "PTHLH"    "CAV1"    
## [1] "ADAM17"
##  [1] "TIMP3"    "MMP1"     "PLAU"     "TNC"      "APP"      "FSTL3"   
##  [7] "SERPINE1" "VIM"      "PTHLH"    "LTBP1"   
## [1] "BMP5"
## [1] "SERPINE1" "LTBP1"    "LAMA3"    "VIM"     
## [1] "FN1"
## [1] "MMP10"    "MMP1"     "SERPINE1" "LTBP1"    "CAV1"     "VIM"     
## [7] "PLAU"     "PTHLH"    "MMP2"    
## [1] "BMP4"
## [1] "MMP2"     "SERPINE1" "P4HA2"    "LTBP1"    "VIM"      "PLAU"    
## [1] "COL4A1"
##  [1] "SERPINE1" "LTBP1"    "PLAU"     "VIM"      "CAV1"     "PTHLH"   
##  [7] "LAMA3"    "ITGA5"    "MMP2"     "APP"      "MMP1"    
## [1] "IBSP"
##  [1] "SERPINE1" "MMP1"     "LTBP1"    "COL1A1"   "MMP2"     "FHL2"    
##  [7] "PLAU"     "PTHLH"    "ITGA5"    "CAV1"     "VIM"      "THBS1"   
## [1] "FBN1"
##  [1] "MMP1"     "SERPINE1" "MMP2"     "LTBP1"    "VIM"      "PLAU"    
##  [7] "CAV1"     "THBS1"    "ITGA5"    "GJA1"    
## [1] "JAM2"
## [1] "SERPINE1" "LTBP1"    "VIM"      "MMP2"     "ITGA5"    "GJA1"    
## [7] "TPM1"     "PLAU"    
## [1] "FGF1"
## [1] "SERPINE1" "PLAU"     "MMP1"     "LTBP1"    "VIM"      "ITGA5"   
## [7] "CAV1"     "MMP2"     "PTHLH"   
## [1] "LAMB2"
##  [1] "SERPINE1" "LTBP1"    "VIM"      "CAV1"     "ITGA5"    "PLAU"    
##  [7] "MMP2"     "MMP1"     "PTHLH"    "FHL2"    
## [1] "IL24"
## [1] "MMP1"     "SERPINE1" "LTBP1"    "ITGA5"    "VIM"      "MMP2"    
## [7] "FHL2"     "PLAU"     "CAV1"    
## [1] "TGFB2"
## [1] "SERPINE1" "COL1A1"   "C1S"      "PTHLH"    "DKK3"

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
```

The putatively active ligand-target links will be visualized in a heatmap.

``` r
p_ligand_target_network = active_ligand_target_links %>% t() %>% make_heatmap_ggplot("Ligand","Target", color = "purple",legend_position = "top", x_axis_position = "top") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01))

p_ligand_target_network
```

![](ligand_activity_geneset_files/figure-markdown_github/unnamed-chunk-10-1.png)

### Infer signaling paths beween ligand(s) and target(s) of interest

As follow-up analysis, you can infer possible signaling paths between ligands and targets of interest. You can read how to do this in the following vignette [Inferring ligand-to-target signaling paths](ligand_target_signaling_path.md):`vignette("ligand_target_signaling_path", package="nichenetr")`.

Another follow-up analysis is getting a "tangible" measure of how well top-ranked ligands predict the gene set of interest and assess which genes of the gene set can be predicted well. You can read how to do this in the following vignette [Assess how well top-ranked ligands can predict a gene set of interest](target_prediction_evaluation_geneset.md):`vignette("target_prediction_evaluation_geneset", package="nichenetr")`.

### References

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel, Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017. “Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24. doi:[10.1016/j.cell.2017.10.044](https://doi.org/10.1016/j.cell.2017.10.044).
