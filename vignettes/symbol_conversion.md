Converting NicheNet’s model from human to mouse symbols
================
Robin Browaeys
2019-07-31

<!-- github markdown built using
rmarkdown::render("vignettes/symbol_conversion.Rmd", output_format = "github_document")
-->

In this vignette, we show how to convert NicheNet’s ligand-target matrix
model from human to mouse gene symbols. This is necessary if you want to
apply NicheNet to mouse expression data, because the NicheNet prior
information is given in human gene symbols (because most data sources at
the basis of NicheNet are based on human data). One-to-one orthologs
were gathered from NCBI HomoloGene and also from ENSEMBL via biomaRt.

### Load required packages

``` r
library(nichenetr)
library(tidyverse)
```

### Load NicheNet’s ligand-target model:

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05
dim(ligand_target_matrix)
## [1] 25345   688
```

### Convert the ligand-target model from human to mouse symbols.

Because not all human genes have a mouse one-to-one ortholog, these
genes will be removed from the mouse
model.

``` r
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

dim(ligand_target_matrix)
## [1] 17330   644
```

Show the top 10 targets of TNF (in mouse
symbols):

``` r
top_targets = extract_top_n_targets("Tnf",10,ligand_target_matrix) %>% names()
top_targets
##  [1] "Hacd4"  "P3h2"   "Sele"   "Vcam1"  "Ubd"    "Ccl19"  "Muc5ac"
##  [8] "Cxcl9"  "Crp"    "Icam1"
```

If you want to convert mouse to human symbols, you can use:

``` r
top_targets %>% convert_mouse_to_human_symbols()
##    Hacd4     P3h2     Sele    Vcam1      Ubd    Ccl19   Muc5ac    Cxcl9 
##  "HACD4"   "P3H2"   "SELE"  "VCAM1"    "UBD"  "CCL19" "MUC5AC"  "CXCL9" 
##      Crp    Icam1 
##    "CRP"  "ICAM1"
```
