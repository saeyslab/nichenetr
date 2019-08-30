Single-cell NicheNet’s ligand activity analysis
================
Robin Browaeys
2018-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/ligand_activity_single_cell.Rmd", output_format = "github_document")
-->

This vignette shows how NicheNet can be used to predict which ligands
might be active in single-cells. If a ligand has a high activity in a
cell, this means that target genes of that ligand are stronger expressed
in that cell than in other cells. In this example, we will use data from
Puram et al. to explore intercellular communication in the tumor
microenvironment in head and neck squamous cell carcinoma (HNSCC) (See
Puram et al. 2017). More specifically, we will assess the activity of
cancer-associated fibroblast (CAF) ligands in malignant cells. The used
ligand-target matrix and example expression data of interacting cells
can be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758).

In order to prioritize ligands regulating a process of interest, you can
perform a regression/correlation analysis between ligand activities in
cells, and scores of a cell corresponding to the process of interest.
For example, in this case study we were interested in finding ligands
regulating p-EMT. Therefore we correlated ligand activities to the p-EMT
scores of cells.

The purpose of this single-cell ligand activity analysis is to offer a
complementary way to prioritize ligands driving the process of interest
and to better analyze heterogeneity in ligand activity between different
cells.

### Load nichenetr and tidyverse

``` r
library(nichenetr)
library(tidyverse)
```

### Read in expression data of interacting cells

First, we will read in the single-cell data from CAF and malignant cells
from HNSCC tumors (See Puram et al.
2017).

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

Secondly, we will determine which genes are expressed in CAFs and
malignant cells from high quality primary tumors. Therefore, we wil not
consider cells from tumor samples of less quality or from lymph node
metastases. To determine expressed genes, we use the definition used by
of Puram et
al.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter((tumor %in% tumors_remove == FALSE)) %>% filter(`non-cancer cell type` == "CAF") %>% .$cell
malignant_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1) %>% filter((tumor %in% tumors_remove == FALSE)) %>% .$cell

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
```

### Load the ligand-target model we want to use

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

### Perform NicheNet’s single-cell ligand activity analysis

In a first step, we will define a set of potentially active ligands. As
potentially active ligands, we will use ligands that are 1) expressed by
CAFs and 2) can bind a (putative) receptor expressed by malignant cells.
Putative ligand-receptor links were gathered from NicheNet’s
ligand-receptor data
sources.

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network$from %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_CAFs)
receptors = lr_network$to %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% .$from %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"
```

In a second step, we will scale the single-cell expression data
(including only expressed
genes).

``` r
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
expression_scaled = expression %>% .[malignant_ids,background_expressed_genes] %>% scale_quantile()
```

Now perform the ligand activity analysis: infer how well NicheNet’s
ligand-target potential scores can predict whether a gene belongs to
most strongly expressed genes in a cell compared to other cells. To
reduce the running time for this vignette, we will perform the analysis
only on 10 example cells from the HN5 tumor. This vignette’s only
purpose is to illustrate the analysis.

In practice, ligand activity analysis for several cells can be better
run in parallel (via
e.g. parallel::mclapply)\!

``` r
malignant_hn5_ids = sample_info %>% filter(tumor == "HN5") %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1)  %>% .$cell %>% head(10)

ligand_activities = predict_single_cell_ligand_activities(cell_ids = malignant_hn5_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

### Ligand prioritization by regression analysis

Furthermore, we will also show how you can perform additional analyses
by linking the ligand activity in cells to other properties of cells in
order to prioritize ligands. As toy example, we will score malignant
cells here on the extent to which they express the core p-EMT gene
“TGFBI”.

``` r
cell_scores_tbl = tibble(cell = malignant_hn5_ids, score = expression_scaled[malignant_hn5_ids,"TGFBI"])
```

Then, we will determine the correlation between these p-EMT scores and
ligand activities over all cells to prioritize p-EMT-inducing ligands.
We hypothesize that ligands might be potential regulators of the p-EMT
program if higher ligand activities are associated with higher p-EMT
scores. Based on this correlation, we obtained a ranking of potential
p-EMT-inducing ligands.

To do so, we frist need to process and normalize the ligand activities
(i.e. pearson correlation values) to make different cells comparable.
Here we use modified z-score
normalization.

``` r
normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)
```

Then, we combine the ligand activities and cell property scores and
perform correlation and regression analysis. We can prioritize ligands
by ranking them based on the pearson correlation between activity scores
and property
scores.

``` r
output_correlation_analysis = single_ligand_activity_score_regression(normalized_ligand_activities,cell_scores_tbl)
output_correlation_analysis %>% arrange(-pearson_regression) %>% select(pearson_regression, ligand)
## # A tibble: 131 x 2
##    pearson_regression ligand  
##                 <dbl> <chr>   
##  1              0.525 TNC     
##  2              0.497 TFPI    
##  3              0.491 SEMA5A  
##  4              0.488 ANXA1   
##  5              0.473 TNFSF13B
##  6              0.462 IBSP    
##  7              0.449 HDGF    
##  8              0.443 HSP90B1 
##  9              0.431 CALM3   
## 10              0.428 CXCL11  
## # ... with 121 more rows
```

Visualize the relation between ligand activity
and the cell's property score of interest

``` r
inner_join(cell_scores_tbl,normalized_ligand_activities) %>% ggplot(aes(score,TNC)) + geom_point() + geom_smooth(method = "lm")
```

![](ligand_activity_single_cell_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

### References

<div id="refs" class="references">

<div id="ref-puram_single-cell_2017">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
