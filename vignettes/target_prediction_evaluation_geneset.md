Assess how well top-ranked ligands can predict a gene set of interest
================
Robin Browaeys
2019-02-19

<!-- github markdown built using
rmarkdown::render("vignettes/target_prediction_evaluation_geneset.Rmd", output_format = "github_document")
-->

This vignette shows how NicheNet can be used to to predict which ligands
might regulate a given set of genes and how well they do this
prediction. For this analysis, you need to define:

  - a set of genes of which expression in a “receiver cell” is possibly
    affected by extracellular protein signals (ligands) (e.g. genes
    differentially expressed upon cell-cell interaction )
  - a set of potentially active ligands (e.g. ligands expressed by
    interacting “sender cells”)

Therefore, you often first need to process expression data of
interacting cells to define both.

In this example, we will use data from Puram et al. to explore
intercellular communication in the tumor microenvironment in head and
neck squamous cell carcinoma (HNSCC) (See Puram et al. 2017). More
specifically, we will look at which ligands expressed by
cancer-associated fibroblasts (CAFs) can induce a specific gene program
in neighboring malignant cells. This program, a partial
epithelial-mesenschymal transition (p-EMT) program, could be linked by
Puram et al. to metastasis.

For this analysis, we will first assess the ligand activity of each
ligand, or in other words, we will assess how well each CAF-ligand can
predict the p-EMT gene set compared to the background of expressed
genes. This allows us to prioritize p-EMT-regulating ligands. Then, we
will assess how well the prioritized ligands together can predict
whether genes belong to the gene set of interest or not.

The used ligand-target matrix and example expression data of interacting
cells can be downloaded from Zenodo.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)

### Load packages required for this vignette

``` r
library(nichenetr)
library(tidyverse)
```

### Read in expression data of interacting cells

First, we will read in the publicly available single-cell data from CAF
and malignant cells from HNSCC
tumors.

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
##   [3]
## names
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

CAF_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% pull(cell)
##   <environment: namespace:dplyr>
## .data
## var
##   <environment: namespace:tidyselect>
## vars
## var
##   <environment: namespace:rlang>
## x
##   <environment: namespace:rlang>
## x
##   <environment: namespace:rlang>
## x
## as_type
##   <environment: namespace:rlang>
## x
## n
## finite
malignant_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
##   <environment: namespace:base>
## x
## trim
## na.rm
## ...
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
```

### Load the ligand-target model we want to use

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
## dim
## dimnames
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##   <environment: namespace:knitr>
## x
## ...
## inline
##   <environment: namespace:evaluate>
## x
## ...
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
```

### Load the gene set of interest and background of genes

As gene set of interest, we consider the genes of which the expression
is possibly affected due to communication with other cells.

Because we here want to investigate how CAF regulate the expression of
p-EMT genes in malignant cells, we will use the p-EMT gene set defined
by Puram et al. as gene set of interset and use all genes expressed in
malignant cells as background of
genes.

``` r
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
##   <environment: namespace:readr>
## na
## quoted_na
## quote
## comment
## trim_ws
## skip_empty_rows
##   <environment: namespace:readr>
## delim
## quote
## na
## quoted_na
## comment
## trim_ws
## escape_double
## escape_backslash
## skip_empty_rows
##   <environment: namespace:readr>
## file
## tokenizer
## col_names
## col_types
## locale
## skip
## skip_empty_rows
## comment
## n_max
## guess_max
## progress
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## path
## input
##   <environment: namespace:readr>
## path
## skip
## skip_empty_rows
## comment
## [1]
## [2]
## [3]
##   <environment: namespace:readr>
## con
##   <environment: namespace:base>
## con
## ...
##   <environment: namespace:base>
## con
## open
## blocking
## ...
##   <environment: namespace:readr>
## con
## filename
## chunk_size
##   <environment: namespace:base>
## con
## what
## n
## size
## signed
## endian
##   <environment: namespace:base>
## e
## f
## onexit
##   <environment: namespace:readr>
## path
## skip
## skip_empty_rows
## comment
## ...
##   <environment: namespace:readr>
## path
##   <environment: namespace:readr>
## type
## x
## skip
## skip_empty_rows
## comment
## ...
##   <environment: namespace:readr>
## file
## col_names
## col_types
## guessed_types
## comment
## skip
## skip_empty_rows
## guess_max
## tokenizer
## locale
## drop_skipped_names
## [1]
## [2]
## [3]
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## col_types
## default
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
##   <environment: namespace:readr>
## type
## ...
##   <environment: namespace:readr>
## file
## skip
## skip_empty_rows
## comment
##   <environment: namespace:readr>
## datasource
## tokenizer
## locale
## guess_max
## max_limit
##   <environment: namespace:readr>
## guess_max
## max_limit
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## sourceSpec
## tokenizerSpec
## locale_
## n
##   <environment: namespace:readr>
##   <environment: namespace:readr>
## date_names
## date_format
## time_format
## decimal_mark
## grouping_mark
## tz
## encoding
## asciify
##   <environment: namespace:readr>
## language
##   [1]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [2]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [3]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [4]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [5]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [8]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [9]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [10]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [11]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [12]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [13]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [14]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [15]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [16]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [17]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [18]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [19]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [20]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [21]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [22]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [23]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [24]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [25]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [26]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [27]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [28]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [29]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [30]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [31]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [32]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [33]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [34]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [35]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [36]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [37]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [38]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [39]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [40]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [41]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [42]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [43]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [44]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [45]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [46]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [47]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [48]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [49]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [50]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [51]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [52]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [53]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [54]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [55]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [56]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [57]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [58]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [59]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [60]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [61]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [62]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [63]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [64]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [65]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [66]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [67]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [68]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [69]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [70]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [71]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [72]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [73]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [74]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [75]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [76]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [77]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [78]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [79]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [80]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [81]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [82]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [83]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [84]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [85]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [86]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [87]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [88]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [89]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [90]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [91]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [92]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [93]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [94]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [95]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [96]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [97]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [98]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [99]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [100]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [101]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [102]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [103]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [104]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [105]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [106]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [107]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [108]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [109]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [110]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [111]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [112]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [113]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [114]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [115]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [116]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [117]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [118]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [119]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [120]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [121]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [122]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [123]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [124]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [125]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [126]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [127]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [128]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [129]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [130]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [131]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [132]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [133]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [134]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [135]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [136]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [137]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [138]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [139]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [140]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [141]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [142]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [143]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [144]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [145]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [146]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [147]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [148]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [149]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [150]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [151]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [152]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [153]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [154]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [155]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [156]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [157]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [158]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [159]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [160]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [161]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [162]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [163]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [164]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [165]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [166]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [167]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [168]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [169]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [170]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [171]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [172]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [173]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [174]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [175]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [176]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [177]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [178]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [179]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [180]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [181]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [182]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [183]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [184]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [185]
## [1]
## [2]
## [3]
## [4]
## [5]
## names
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## x
##   <environment: namespace:base>
## tzdir
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## name
##   <environment: namespace:readr>
##   <environment: namespace:readr>
## spec
## n
##   <environment: namespace:readr>
## x
## n
## condense
## colour
## ...
## [1]
## [2]
## [3]
##   <environment: namespace:readr>
## a
## b
##   <environment: namespace:readr>
## expr
## ...
## sep
##   <environment: namespace:readr>
## cols
## colourise
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## data
## tokenizer
## col_specs
## col_names
## locale_
## n_max
## progress
##   <environment: namespace:readr>
## sourceSpec
## tokenizerSpec
## colSpecs
## colNames
## locale_
## n_max
## progress
##   <environment: namespace:readr>
##   <environment: namespace:readr>
## x
## all_colnames
## name
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## x
##   <environment: namespace:readr>
## x
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
head(pemt_geneset)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"
```

### Perform NicheNet’s ligand activity analysis on the gene set of interest

In a first step, we will define a set of potentially active ligands. As
potentially active ligands, we will use ligands that are 1) expressed by
CAFs and 2) can bind a (putative) receptor expressed by malignant cells.
Putative ligand-receptor links were gathered from NicheNet’s
ligand-receptor data
sources.

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
##   [1]
##   [2]
##   [3]
##   [4]
## row.names
## class
## names

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_CAFs)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"
```

Now perform the ligand activity analysis: infer how well NicheNet’s
ligand-target potential scores can predict whether a gene belongs to the
p-EMT program or
not.

``` r
ligand_activities = predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
##   <environment: namespace:nichenetr>
## geneset
## background_expressed_genes
## ligand_target_matrix
## potential_ligands
## single
## ...
##   <environment: namespace:nichenetr>
## gene_list
## name
## ligands_oi
## background
##   <environment: namespace:nichenetr>
## settings
## all_ligands
## validation
## single
##   <environment: namespace:nichenetr>
## setting
## test_ligand
##   <environment: namespace:nichenetr>
## setting
## ligand_target_matrix
## ligands_position
## known
##   <environment: namespace:nichenetr>
## setting
## ligand_target_matrix
## ligands_position
##   <environment: namespace:nichenetr>
## response
## prediction
## continuous
## prediction_response_df
##   <environment: namespace:nichenetr>
## prediction
## response
## iregulon
##   <environment: namespace:ROCR>
## predictions
## labels
## label.ordering
##   <environment: namespace:base>
## x
## ...
##   <environment: namespace:base>
## e1
## e2
##   <environment: namespace:base>
## e1
## e2
## [1]
## [2]
## [3]
##   <environment: namespace:ROCR>
## predictions
## labels
##   <environment: namespace:base>
## x
## ...
## drop
##   <environment: namespace:methods>
## object
## name
##   <environment: namespace:ROCR>
## prediction.obj
## measure
## x.measure
## ...
##   <environment: namespace:ROCR>
##   <environment: namespace:ROCR>
## arglist
## args.to.select
## complement
##   <environment: namespace:ROCR>
## arglist
## ...
##   <environment: namespace:ROCR>
## predictions
## labels
## cutoffs
## fp
## tp
## fn
## tn
## n.pos
## n.neg
## n.pos.pred
## n.neg.pred
##   <environment: namespace:ROCR>
## predictions
## labels
## cutoffs
## fp
## tp
## fn
## tn
## n.pos
## n.neg
## n.pos.pred
## n.neg.pred
##   <environment: namespace:ROCR>
## p.obj.1
## p.obj.2
##   <environment: namespace:stats>
## x
## y
## method
## yleft
## yright
## rule
## f
## ties
## [1]
## [2]
## [3]
##   <environment: namespace:stats>
## x
## y
## ties
## warn.collapsing
##   <environment: namespace:grDevices>
## x
## y
## xlab
## ylab
## log
## recycle
## setLab
##   <environment: namespace:base>
## x
## na.rm
## strictly
##   <environment: namespace:stats>
## x
## y
## v
## method
## yleft
## yright
## f
##   <environment: namespace:ROCR>
## predictions
## labels
## cutoffs
## fp
## tp
## fn
## tn
## n.pos
## n.neg
## n.pos.pred
## n.neg.pred
##   <environment: namespace:tidyr>
## data
## replace
## ...
##   <environment: namespace:tidyr>
## x
## var
##   <environment: namespace:tidyr>
## x
##   <environment: namespace:base>
## x
## i
## j
## value
##   <environment: namespace:caTools>
## x
## y
##   <environment: namespace:nichenetr>
## observed
## known
##   <environment: namespace:ROCR>
## predictions
## labels
## cutoffs
## fp
## tp
## fn
## tn
## n.pos
## n.neg
## n.pos.pred
## n.neg.pred
##   <environment: namespace:stats>
## x
## y
## use
## method
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:base>
## x
## na.last
## ties.method
##   <environment: namespace:stats>
## x
## y
## alternative
## method
## exact
## conf.level
## continuity
## ...
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:stats>
## ...
##   <environment: namespace:stats>
## p
## mean
## sd
## lower.tail
## log.p
##   <environment: namespace:stats>
## q
## df
## ncp
## lower.tail
## log.p
##   <environment: namespace:limma>
## index
## statistics
## ...
##   <environment: namespace:limma>
## index
## statistics
## alternative
## type
## ranks.only
## nsim
## [1]
## [2]
## [3]
##   <environment: namespace:limma>
## index
## statistics
## correlation
## df
##   <environment: namespace:nichenetr>
## prior
## response
##   <environment: namespace:nichenetr>
## x
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:base>
## x
##   <environment: namespace:data.table>
## ...
## keep.rownames
## check.names
## key
## stringsAsFactors
## [1]
## [2]
## [3]
##   <environment: namespace:data.table>
##   <environment: namespace:data.table>
## x
## name
## value
##   <environment: namespace:data.table>
## DT
## n
## verbose
##   <environment: namespace:data.table>
## ...
##   <environment: namespace:data.table>
## x
## table
##   <environment: namespace:data.table>
## x
## keep.rownames
## ...
##   <environment: namespace:data.table>
## x
## [1]
## [2]
## [3]
##   <environment: namespace:data.table>
## x
##   <environment: namespace:data.table>
## x
## class
##   <environment: namespace:data.table>
## x
## table
## nomatch
##   <environment: namespace:data.table>
## x
##   <environment: namespace:data.table>
## x
##   <environment: namespace:data.table>
## x
##   <environment: namespace:data.table>
## x
##   <environment: namespace:data.table>
## x
## subset
## select
## ...
##   <environment: namespace:data.table>
## x
## i
## j
## by
## keyby
## with
## nomatch
## mult
## roll
## rollends
## which
## .SDcols
## verbose
## allow.cartesian
## drop
## on
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:data.table>
## n
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:data.table>
## e
##   <environment: namespace:base>
## expr
## functions
## max.names
## unique
##   <environment: namespace:data.table>
## isub
## x
## enclos
## notjoin
## verbose
##   <environment: namespace:nichenetr>
## oneRanking
## aucThreshold
## maxAUC
##   <environment: namespace:data.table>
## x
## ...
##   <environment: namespace:data.table>
## x
## cols
##   <environment: namespace:data.table>
## x
## cols
## retain.key
## unlock
```

Now, we want to rank the ligands based on their ligand activity. In our
validation study, we showed that the pearson correlation between a
ligand’s target predictions and the observed transcriptional response
was the most informative measure to define ligand activity. Therefore,
we will rank the ligands based on their pearson correlation coefficient.

``` r
ligand_activities %>% arrange(-pearson)
##   <environment: namespace:pillar>
## x
## ...
## sigfig
##   <environment: namespace:pillar>
## x
## sigfig
## ...
##   <environment: namespace:pillar>
## x
## sigfig
## scientific
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## digits
##   <environment: namespace:pillar>
## x
## sigfig
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## sigfig
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## s
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## .x
## .f
## ...
##   <environment: namespace:pillar>
## x
## negative
## significant
##   <environment: namespace:pillar>
## x
## negative
##   <environment: namespace:pillar>
## s
##   <environment: namespace:pillar>
## s
##   <environment: namespace:pillar>
## s
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:tibble>
## message
## n
##   <environment: namespace:tibble>
## x
## # A tibble: 131 x 4
##    test_ligand auroc   aupr pearson
##    <chr>       <dbl>  <dbl>   <dbl>
##  1 PTHLH       0.667 0.0720   0.128
##  2 CXCL12      0.680 0.0507   0.123
##  3 AGT         0.676 0.0581   0.120
##  4 TGFB3       0.689 0.0454   0.117
##  5 IL6         0.693 0.0510   0.115
##  6 INHBA       0.695 0.0502   0.113
##  7 ADAM17      0.672 0.0526   0.113
##  8 TNC         0.700 0.0444   0.109
##  9 CTGF        0.680 0.0473   0.108
## 10 FN1         0.679 0.0505   0.108
## # ... with 121 more rows
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
##   <environment: namespace:rlang>
## x
## finite
head(best_upstream_ligands)
## [1] "PTHLH"  "CXCL12" "AGT"    "TGFB3"  "IL6"    "INHBA"
```

For the top 20 ligands, we will now build a multi-ligand model that uses
all top-ranked ligands to predict whether a gene belongs to the p-EMT
program of not. This classification model will be trained via
cross-validation and returns a probability for every
gene.

``` r
# change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
k = 3 # 3-fold
n = 2 # 2 rounds

pemt_gene_predictions_top20_list = seq(n) %>% lapply(assess_rf_class_probabilities, folds = k, geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligands_oi = best_upstream_ligands, ligand_target_matrix = ligand_target_matrix)
##   <environment: namespace:nichenetr>
## round
## folds
## geneset
## background_expressed_genes
## ligands_oi
## ligand_target_matrix
##   <environment: namespace:base>
## seed
## kind
## normal.kind
## sample.kind
##   <environment: namespace:nichenetr>
## i
## affected_genes_grouped
## strict_background_expressed_genes_grouped
## best_upstream_ligands
## ligand_target_matrix
##   <environment: namespace:nichenetr>
## setting
## ligand_target_matrix
## ligands_position
## ntrees
## mtry
## continuous
## known
## filter_genes
##   <environment: namespace:rlang>
## x
##   <environment: namespace:rlang>
## x
## ns
## private
##   <environment: namespace:rlang>
## x
## ns
## private
##   <environment: namespace:tidyselect>
## expr
##   <environment: namespace:randomForest>
## x
## y
## xtest
## ytest
## ntree
## mtry
## replace
## classwt
## cutoff
## strata
## sampsize
## nodesize
## maxnodes
## importance
## localImp
## nPerm
## proximity
## oob.prox
## norm.votes
## do.trace
## keep.forest
## corr.bias
## keep.inbag
## ...
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:base>
## x
##   <environment: namespace:randomForest>
## x
##   <environment: namespace:base>
## frame
## rownames.force
##   <environment: namespace:base>
## length
##   <environment: namespace:base>
## x
##   <environment: namespace:base>
## x
## i
## j
## ...
## drop
##   <environment: namespace:nichenetr>
## setting
## rf_model
## ligand_target_matrix
## ligands_position
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
##   [8]
## names
## class
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
##   [8]
## names
## class
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
##   [8]
## names
## class
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
##   [8]
## names
## class
##   [1]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
##   [2]
##   [3]
##   [4]
##   [5]
##   [6]
##   [7]
##   [8]
##   [9]
##   [10]
## [1]
##   [11]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
##   [8]
## names
## class
##   [1]
## names
##   [1]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
## [11]
## [12]
## [13]
## [14]
## [15]
## [16]
## [17]
## [18]
## [19]
## [20]
## [21]
## [22]
## [23]
## [24]
## [25]
## [26]
## [27]
## [28]
## [29]
## [30]
## [31]
## [32]
## [33]
## [34]
## [35]
## [36]
## [37]
## [38]
## [39]
## [40]
## [41]
## [42]
## [43]
## [44]
## [45]
## [46]
## [47]
## [48]
## [49]
## [50]
## [51]
## [52]
## [53]
## [54]
## [55]
## [56]
## [57]
## [58]
## [59]
## [60]
## [61]
## [62]
## [63]
## [64]
## [65]
## [66]
## [67]
## [68]
## [69]
## [70]
## [71]
## [72]
## [73]
## [74]
## [75]
## [76]
## [77]
## [78]
## [79]
## [80]
## [81]
## [82]
## [83]
## [84]
## [85]
## [86]
## [87]
## [88]
## [89]
## [90]
## [91]
## [92]
## [93]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## names
##   <environment: namespace:e1071>
## x
## ...
##   <environment: namespace:e1071>
## x
## ...
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
##   [8]
## names
## class
##   <environment: namespace:randomForest>
## object
## newdata
## type
## norm.votes
## predict.all
## proximity
## nodes
## cutoff
## ...
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:base>
## x
## table
## nomatch
##   <environment: namespace:base>
## x
## MARGIN
## STATS
## FUN
## check.margin
## ...
```

Evaluate now how well the target gene probabilies accord to the gene set
assignments

``` r
# get performance: auroc-aupr-pearson
target_prediction_performances_cv = pemt_gene_predictions_top20_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>% bind_rows() %>% mutate(round=seq(1:nrow(.)))
##   <environment: namespace:nichenetr>
## response_prediction_tibble
```

What is the AUROC, AUPR and PCC of this model (averaged over
cross-validation rounds)?

``` r
target_prediction_performances_cv$auroc %>% mean()
## [1] 0.7295863
target_prediction_performances_cv$aupr %>% mean()
## [1] 0.07603073
target_prediction_performances_cv$pearson %>% mean()
## [1] 0.1660327
```

Evaluate now whether genes belonging to the gene set are more likely to
be top-predicted. We will look at the top 5% of predicted targets
here.

``` r
# get performance: how many p-EMT genes and non-p-EMT-genes among top 5% predicted targets
target_prediction_performances_discrete_cv = pemt_gene_predictions_top20_list %>% lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>% bind_rows() %>% ungroup() %>% mutate(round=rep(1:length(pemt_gene_predictions_top20_list), each = 2))
##   <environment: namespace:nichenetr>
## affected_gene_predictions
## quantile_cutoff
##   <environment: namespace:stats>
## x
## probs
## na.rm
## names
## type
## ...
##   <environment: namespace:stats>
## x
## digits
## probability
## use.fC
## ...
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:dplyr>
## .data
##   <environment: namespace:dplyr>
## data
##   <environment: namespace:dplyr>
## x
## wt
## sort
## name
##   <environment: namespace:rlang>
## quo
##   <environment: namespace:dplyr>
## x
## name
##   <environment: namespace:dplyr>
## .data
## ...
##   <environment: namespace:dplyr>
## df
## dots
## frame
## caller_env
##   <environment: namespace:dplyr>
## ...
## .drop
##   <environment: namespace:dplyr>
## .data
## ...
##   <environment: namespace:dplyr>
## .tbl
## [1]
## [2]
## [3]
##   <environment: namespace:dplyr>
## index
## mask_proxy_xp
## [1]
## [2]
## [3]
##   <environment: namespace:dplyr>
## idx
## mask_proxy_xp
```

What is the fraction of p-EMT genes that belongs to the top 5% predicted
targets?

``` r
target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean()
## [1] 0.25
```

What is the fraction of non-p-EMT genes that belongs to the top 5%
predicted
targets?

``` r
target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean()
## [1] 0.04769076
```

We see that the p-EMT genes are enriched in the top-predicted target
genes. To test this, we will now apply a Fisher’s exact test for every
cross-validation round and report the average
p-value.

``` r
target_prediction_performances_discrete_fisher = pemt_gene_predictions_top20_list %>% lapply(calculate_fraction_top_predicted_fisher, quantile_cutoff = 0.95) 
##   <environment: namespace:nichenetr>
## affected_gene_predictions
## quantile_cutoff
## p_value_output
##   <environment: namespace:dplyr>
## grouping_data
## frame
##   <environment: namespace:stats>
## x
## y
## workspace
## hybrid
## hybridPars
## control
## or
## alternative
## conf.int
## conf.level
## simulate.p.value
## B
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:base>
## input
## target
## nomatch
##   <environment: namespace:stats>
## x
## m
## n
## k
## log
##   <environment: namespace:stats>
## q
## m
## n
## k
## lower.tail
## log.p
##   <environment: namespace:stats>
## f
## interval
## ...
## lower
## upper
## f.lower
## f.upper
## extendInt
## check.conv
## tol
## maxiter
## trace
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:base>
## ...
## na.rm
target_prediction_performances_discrete_fisher %>% unlist() %>% mean()
## [1] 5.647773e-10
```

Finally, we will look at which p-EMT genes are well-predicted in every
cross-validation round.

``` r
# get top predicted genes
top_predicted_genes = seq(length(pemt_gene_predictions_top20_list)) %>% lapply(get_top_predicted_genes,pemt_gene_predictions_top20_list) %>% reduce(full_join, by = c("gene","true_target"))
##   <environment: namespace:nichenetr>
## round
## gene_prediction_list
## quantile_cutoff
##   <environment: namespace:purrr>
## .x
## .f
## ...
## .init
## .dir
## .acc
## [1]
## [2]
## [3]
##   <environment: namespace:rlang>
## arg
## values
##   <environment: namespace:purrr>
## x
## init
## left
##   <environment: namespace:purrr>
## x
## init
## left
##   <environment: namespace:purrr>
## start
## end
##   <environment: namespace:dplyr>
## x
## y
## by
## copy
## suffix
## ...
## na_matches
##   <environment: namespace:dplyr>
## x
## y
## by_x
## by_y
## aux_x
## aux_y
## na_match
## frame
##   <environment: namespace:rlang>
## x
## empty
top_predicted_genes %>% filter(true_target)
##   <environment: namespace:pillar>
## widths
## tier_widths
##   <environment: namespace:pillar>
## widths
## tier_widths
## widths_offset
## tier_widths_offset
##   <environment: namespace:pillar>
## x
## ...
##   <environment: namespace:fansi>
## x
## # A tibble: 28 x 4
##    gene    true_target predicted_top_target_roun~ predicted_top_target_rou~
##    <chr>   <lgl>       <lgl>                      <lgl>                    
##  1 COL1A1  TRUE        TRUE                       TRUE                     
##  2 MMP2    TRUE        TRUE                       TRUE                     
##  3 MMP1    TRUE        TRUE                       TRUE                     
##  4 PLAU    TRUE        TRUE                       TRUE                     
##  5 TIMP3   TRUE        TRUE                       TRUE                     
##  6 MT2A    TRUE        TRUE                       TRUE                     
##  7 INHBA   TRUE        TRUE                       TRUE                     
##  8 COL4A2  TRUE        TRUE                       TRUE                     
##  9 MMP10   TRUE        TRUE                       TRUE                     
## 10 COL17A1 TRUE        TRUE                       TRUE                     
## # ... with 18 more rows
```

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
