Perform NicheNet analysis starting from a Seurat object
================
Robin Browaeys
2019-11-08

<!-- github markdown built using 
rmarkdown::render("vignettes/seurat_wrapper.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform a basic NicheNet analysis
on a Seurat v3 object. Such a NicheNet analysis can help you to generate
hypotheses about an intercellular communication process of interest for
which you have single-cell gene expression data as a Seurat object.
Specifically, NicheNet can predict 1) which ligands from one or more
cell population(s) (“sender/niche”) are most likely to affect target
gene expression in an interacting cell population (“receiver/target”)
and 2) which specific target genes are affected by which of these
predicted ligands.

Because NicheNet studies how ligands affect gene expression in
putatively neighboring/interacting cells, you need to have data about
this effect in gene expression you want to study. So, there need to be
‘some kind of’ differential expression in a receiver cell population,
caused by ligands from one of more interacting sender cell populations.

In this vignette, we demonstrate the use of NicheNet on a Seurat Object.
The wrapper function we will show consists of the same different steps
that are discussed in detail in the main, basis, NicheNet vignette
[NicheNet’s ligand activity analysis on a gene set of interest: predict
active ligands and their target
genes](ligand_activity_geneset.md):`vignette("ligand_activity_geneset", package="nichenetr")`.
Make sure you understand the different steps in a NicheNet analysis that
are described in that vignette before proceeding with this vignette and
performing a real NicheNet analysis on your data. In another vignette
[Perform NicheNet analysis starting from a Seurat object: step-by-step
analysis](seurat_steps.md):`vignette("seurat_steps", package="nichenetr")`,
we also show the execution of these steps one for one, but in contrast
to the main vignette now specifically for a Seurat Object. This allows
users to adapt specific steps of the pipeline to make them more
appropriate for their data (recommended).

As example expression data of interacting cells, we will use mouse
NICHE-seq data from Medaglia et al. to explore intercellular
communication in the T cell area in the inguinal lymph node before and
72 hours after lymphocytic choriomeningitis virus (LCMV) infection (See
Medaglia et al. 2017). We will NicheNet to explore immune cell crosstalk
in response to this LCMV infection.

In this dataset, differential expression is observed between CD8 T cells
in steady-state and CD8 T cells after LCMV infection. NicheNet can be
applied to look at how several immune cell populations in the lymph node
(i.e., monocytes, dendritic cells, NK cells, B cells, CD4 T cells) can
regulate and induce these observed gene expression changes. NicheNet
will specifically prioritize ligands from these immune cells and their
target genes that change in expression upon LCMV infection.

The used NicheNet networks, ligand-target matrix and example expression
data of interacting cells can be downloaded from Zenodo. The NicheNet
networks and ligand-target matrix at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758)
and the Seurat object of the processed NICHE-seq single-cell data at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3531889.svg)](https://doi.org/10.5281/zenodo.3531889).

# Prepare NicheNet analysis

## Load required packages, read in the Seurat object with processed expression data of interacting cells and NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks.

The NicheNet ligand-receptor network and weighted networks are necessary
to define and show possible ligand-receptor interactions between two
cell populations. The ligand-target matrix denotes the prior potential
that particular ligands might regulate the expression of particular
target genes. This matrix is necessary to prioritize possible
ligand-receptor interactions based on observed gene expression effects
(i.e. NicheNet’s ligand activity analysis) and infer affected target
genes of these prioritized ligands.

### Load Packages:

``` r
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
```

If you would use and load other packages, we recommend to load these 3
packages after the others.

### Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                 CXCL1        CXCL2        CXCL3        CXCL5         PPBP
## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04 2.628388e-04
## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04 1.231819e-04
## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04 3.211944e-04
## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04 4.590521e-04
## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05 5.120439e-05

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)
## # A tibble: 6 x 4
##   from  to    source         database
##   <chr> <chr> <chr>          <chr>   
## 1 CXCL1 CXCR2 kegg_cytokines kegg    
## 2 CXCL2 CXCR2 kegg_cytokines kegg    
## 3 CXCL3 CXCR2 kegg_cytokines kegg    
## 4 CXCL5 CXCR2 kegg_cytokines kegg    
## 5 PPBP  CXCR2 kegg_cytokines kegg    
## 6 CXCL6 CXCR2 kegg_cytokines kegg

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  ABCC6  0.422 
## 2 A1BG  ACE2   0.101 
## 3 A1BG  ADAM10 0.0970
## 4 A1BG  AGO1   0.0525
## 5 A1BG  AKT1   0.0855
## 6 A1BG  ANXA7  0.457
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
## # A tibble: 6 x 3
##   from  to     weight
##   <chr> <chr>   <dbl>
## 1 A1BG  A2M    0.0294
## 2 AAAS  GFAP   0.0290
## 3 AADAC CYP3A4 0.0422
## 4 AADAC IRF8   0.0275
## 5 AATF  ATM    0.0330
## 6 AATF  ATR    0.0355
```

### Read in the expression data of interacting cells:

The dataset used here is publicly available single-cell data from immune
cells in the T cell area of the inguinal lymph node. The data was
processed and aggregated by applying the Seurat alignment pipeline. The
Seurat object contains this aggregated data. Note that this should be a
Seurat v3 object and that gene should be named by their official
mouse/human gene symbol.

``` r
seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj@meta.data %>% head()
##         nGene nUMI orig.ident aggregate res.0.6 celltype nCount_RNA nFeature_RNA
## W380370   880 1611      LN_SS        SS       1    CD8 T       1607          876
## W380372   541  891      LN_SS        SS       0    CD4 T        885          536
## W380374   742 1229      LN_SS        SS       0    CD4 T       1223          737
## W380378   847 1546      LN_SS        SS       1    CD8 T       1537          838
## W380379   839 1606      LN_SS        SS       0    CD4 T       1603          836
## W380381   517  844      LN_SS        SS       0    CD4 T        840          513
```

Visualize which cell populations are present: CD4 T cells (including
regulatory T cells), CD8 T cells, B cells, NK cells, dendritic cells
(DCs) and inflammatory monocytes

``` r
seuratObj@meta.data$celltype %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application
## .
##     B CD4 T CD8 T    DC  Mono    NK  Treg 
##   382  2562  1645    18    90   131   199
DimPlot(seuratObj, reduction = "tsne")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Visualize the data to see to which condition cells belong. The metadata
dataframe column that denotes the condition (steady-state or after LCMV
infection) is here called ‘aggregate.’

``` r
seuratObj@meta.data$aggregate %>% table()
## .
## LCMV   SS 
## 3886 1141
DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Perform the NicheNet analysis

In this case study, we want to apply NicheNet to predict which ligands
expressed by all immune cells in the T cell area of the lymph node are
most likely to have induced the differential expression in CD8 T cells
after LCMV infection.

As described in the main vignette, the pipeline of a basic NicheNet
analysis consist of the following steps:

-   1.  Define a “sender/niche” cell population and a “receiver/target”
        cell population present in your expression data and determine
        which genes are expressed in both populations

-   1.  Define a gene set of interest: these are the genes in the
        “receiver/target” cell population that are potentially affected
        by ligands expressed by interacting cells (e.g. genes
        differentially expressed upon cell-cell interaction)

-   1.  Define a set of potential ligands: these are ligands that are
        expressed by the “sender/niche” cell population and bind a
        (putative) receptor expressed by the “receiver/target”
        population

-   1.  Perform NicheNet ligand activity analysis: rank the potential
        ligands based on the presence of their target genes in the gene
        set of interest (compared to the background set of genes)

-   1.  Infer receptors and top-predicted target genes of ligands that
        are top-ranked in the ligand activity analysis

All these steps are contained in one of three following similar single
functions: `nichenet_seuratobj_aggregate`,
`nichenet_seuratobj_cluster_de` and
`nichenet_seuratobj_aggregate_cluster_de`.

## NicheNet analysis on Seurat object: explain differential expression between two conditions

In this case study, the receiver cell population is the ‘CD8 T’ cell
population, whereas the sender cell populations are ‘CD4 T,’ ‘Treg,’
‘Mono,’ ‘NK,’ ‘B’ and ‘DC.’ The above described functions will consider
a gene to be expressed when it is expressed in at least a predefined
fraction of cells in one cluster (default: 10%).

The gene set of interest are the genes differentially expressed in CD8 T
cells after LCMV infection. The condition of interest is thus ‘LCMV,’
whereas the reference/steady-state condition is ‘SS.’ The notion of
conditions can be extracted from the metadata column ‘aggregate,’ the
method to calculate the differential expression is the standard Seurat
Wilcoxon test.

The number of top-ranked ligands that are further used to predict active
target genes and construct an active ligand-receptor network is 20 by
default.

To perform the NicheNet analysis with these specifications, run the
following:

``` r
# indicated cell types should be cell class identities
# check via: 
# seuratObj %>% Idents() %>% table()
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "CD8 T", 
  condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", 
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
```

### Interpret the NicheNet analysis output

#### Ligand activity analysis results

A first thing NicheNet does, is prioritizing ligands based on predicted
ligand activity. To see the ranking of these ligands, run the following
command:

``` r
nichenet_output$ligand_activities
## # A tibble: 44 x 6
##    test_ligand auroc  aupr pearson  rank bona_fide_ligand
##    <chr>       <dbl> <dbl>   <dbl> <dbl> <lgl>           
##  1 Ebi3        0.638 0.234  0.197      1 FALSE           
##  2 Il15        0.582 0.163  0.0961     2 TRUE            
##  3 Crlf2       0.549 0.163  0.0758     3 FALSE           
##  4 App         0.499 0.141  0.0655     4 TRUE            
##  5 Tgfb1       0.494 0.140  0.0558     5 TRUE            
##  6 Ptprc       0.536 0.149  0.0554     6 TRUE            
##  7 H2-M3       0.525 0.157  0.0528     7 TRUE            
##  8 Icam1       0.543 0.142  0.0486     8 TRUE            
##  9 Cxcl10      0.531 0.141  0.0408     9 TRUE            
## 10 Adam17      0.517 0.137  0.0359    10 TRUE            
## # ... with 34 more rows
```

The different ligand activity measures (auroc, aupr, pearson correlation
coefficient) are a measure for how well a ligand can predict the
observed differentially expressed genes compared to the background of
expressed genes. In our validation study, we showed that the pearson
correlation coefficient between a ligand’s target predictions and the
observed transcriptional response was the most informative measure to
define ligand activity. Therefore, NicheNet ranks the ligands based on
their pearson correlation coefficient. This allows us to prioritize
ligands inducing the antiviral response in CD8 T cells.

The column ‘bona\_fide\_ligand’ indicates whether the ligand is part of
ligand-receptor interactions that are documented in public databases
(‘bona\_fide\_ligand = TRUE’) and not of ligand-receptor interactions
that we predicted based on annotation as ligand/receptor and
protein-protein interaction databases (‘bona\_fide\_ligand = FALSE’).

To get a list of the 20 top-ranked ligands: run the following command

``` r
nichenet_output$top_ligands
##  [1] "Ebi3"   "Il15"   "Crlf2"  "App"    "Tgfb1"  "Ptprc"  "H2-M3"  "Icam1"  "Cxcl10" "Adam17"
## [11] "Cxcl11" "Cxcl9"  "H2-T23" "Sema4d" "Ccl5"   "C3"     "Cxcl16" "Itgb1"  "Anxa1"  "Sell"
```

These ligands are expressed by one or more of the input sender cells. To
see which cell population expresses which of these top-ranked ligands,
you can run the following:

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

As you can see, most op the top-ranked ligands seem to be mainly
expressed by dendritic cells and monocytes.

It could also be interesting to see whether some of these ligands are
differentially expressed after LCMV infection.

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), split.by = "aggregate") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
VlnPlot(seuratObj, features = c("Il15", "Cxcl10","Cxcl16"), split.by = "aggregate",    pt.size = 0, combine = FALSE)
## [[1]]
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

    ## 
    ## [[2]]

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

    ## 
    ## [[3]]

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

#### Inferred active ligand-target links

NicheNet also infers active target genes of these top-ranked ligands. To
see which top-ranked ligands are predicted to have regulated the
expression of which differentially expressed genes, you can run
following command for a heatmap visualization:

``` r
nichenet_output$ligand_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

This is a normal ggplot object that can be adapted likewise. For example
if you want to change the color code to blue instead of purple, change
the axis ticks of the legend, and change the axis labels of the heatmap,
you can do the following:

``` r
nichenet_output$ligand_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + xlab("anti-LCMV response genes in CD8 T cells") + ylab("Prioritized immmune cell ligands")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

If you want, you can also extract the ligand-target links and their
regulatory potential scores in matrix or data frame format (e.g. for
visualization in other ways or output to a csv file).

``` r
nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
##              Cd274 Cd53       Ddit4         Id3 Ifit3        Irf1
## Sell   0.000000000    0 0.001290863 0.001222706     0 0.001095100
## Itgb1  0.000000000    0 0.001162142 0.001214922     0 0.001069406
## C3     0.000000000    0 0.001105490 0.000000000     0 0.000000000
## Ccl5   0.000000000    0 0.001281096 0.001228147     0 0.001155790
## Sema4d 0.000000000    0 0.001103465 0.001179496     0 0.000000000
## H2.T23 0.000000000    0 0.001112018 0.001110184     0 0.000000000
## Adam17 0.002280965    0 0.001760241 0.001546186     0 0.001637201
## Cxcl10 0.000000000    0 0.001354334 0.001372142     0 0.001393116
## Icam1  0.000000000    0 0.001325195 0.001314746     0 0.001375860
## H2.M3  0.000000000    0 0.001436893 0.001506164     0 0.001329158
```

``` r
nichenet_output$ligand_target_df # weight column = regulatory potential
## # A tibble: 155 x 3
##    ligand target  weight
##    <chr>  <chr>    <dbl>
##  1 Ebi3   Cd274  0.00325
##  2 Ebi3   Cd53   0.00321
##  3 Ebi3   Ddit4  0.00335
##  4 Ebi3   Id3    0.00373
##  5 Ebi3   Ifit3  0.00320
##  6 Ebi3   Irf1   0.00692
##  7 Ebi3   Irf7   0.00312
##  8 Ebi3   Irf9   0.00543
##  9 Ebi3   Parp14 0.00336
## 10 Ebi3   Pdcd4  0.00335
## # ... with 145 more rows
```

To get a list of the top-predicted target genes of the 20 top-ranked
ligands: run the following command

``` r
nichenet_output$top_targets
##  [1] "Cd274"  "Cd53"   "Ddit4"  "Id3"    "Ifit3"  "Irf1"   "Irf7"   "Irf9"   "Parp14" "Pdcd4" 
## [11] "Pml"    "Psmb9"  "Rnf213" "Stat1"  "Stat2"  "Tap1"   "Ubc"    "Zbp1"   "Cd69"   "Gbp4"  
## [21] "Basp1"  "Casp8"  "Cxcl10" "Nlrc5"  "Vim"    "Actb"   "Ifih1"  "Myh9"   "B2m"    "H2-T23"
## [31] "Rpl13a" "Cxcr4"
```

You can visualize the expression of these as well. Because we only focus
on CD8 T cells as receiver cells, we will only show expression in these
cells. To emphasize that these target genes are differentially
expressed, we split cells up in steadys-state cells and cells after
response to LCMV infection.

``` r
DotPlot(seuratObj %>% subset(idents = "CD8 T"), features = nichenet_output$top_targets %>% rev(), split.by = "aggregate") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
VlnPlot(seuratObj %>% subset(idents = "CD8 T"), features = c("Zbp1","Ifit3","Irf7"), split.by = "aggregate",    pt.size = 0, combine = FALSE)
## [[1]]
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

    ## 
    ## [[2]]

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

    ## 
    ## [[3]]

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

To visualize both ligand activities and target genes of ligands, run the
following command

``` r
nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

#### Inferred ligand-receptor interactions for top-ranked ligands

NicheNet also infers the receiver cell receptors of these top-ranked
ligands. You can run following command for a heatmap visualization of
the ligand-receptor links:

``` r
nichenet_output$ligand_receptor_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

If you want, you can also extract the ligand-receptor links and their
interaction confidence scores in matrix or data frame format (e.g. for
visualization in other ways or output to a csv file).

``` r
nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]
##           Cxcl9    Cxcl16    Cxcl11      Ccl5    Cxcl10       App
## Cxcr6 0.3629049 0.6598705 0.2255185 0.2627207 0.4001071 0.2255185
## Ccr7  0.2217117 0.2217117 0.3567789 0.3933531 0.2582858 0.2217117
## Ccr9  0.1357118 0.2374693 0.2374693 0.1357118 0.1357118 0.1357118
## Gpr18 0.1374828 0.1374828 0.1374828 0.1374828 0.1374828 0.0000000
## S1pr1 0.1263826 0.1263826 0.1263826 0.1263826 0.1263826 0.0000000
## Itga4 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Cd47  0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Ptk2b 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Il2rb 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## Il2rg 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
```

``` r
nichenet_output$ligand_receptor_df # weight column accords to number of data sources that document this interaction
## # A tibble: 61 x 3
##    ligand receptor weight
##    <chr>  <chr>     <dbl>
##  1 Adam17 Notch1    0.482
##  2 Anxa1  Ccr7      0.222
##  3 Anxa1  Ccr9      0.237
##  4 Anxa1  Cxcr6     0.226
##  5 Anxa1  Itga4     0.201
##  6 App    Ccr7      0.222
##  7 App    Ccr9      0.136
##  8 App    Cxcr6     0.226
##  9 App    Notch1    0.354
## 10 App    Tgfbr2    0.441
## # ... with 51 more rows
```

To get a list of the receptors of the 20 top-ranked ligands: run the
following command

``` r
nichenet_output$top_receptors
##  [1] "Notch1" "Ccr7"   "Ccr9"   "Cxcr6"  "Itga4"  "Tgfbr2" "Itgb2"  "Gpr18"  "S1pr1"  "Il7r"  
## [11] "Il27ra" "Cd8a"   "Klrd1"  "Il2rg"  "Itgal"  "Spn"    "Il2rb"  "Cd47"   "Ptk2b"  "Cd2"   
## [21] "Cd28"   "Selplg" "Ptprc"
```

You can visualize the expression of these as well. Because we only focus
on CD8 T cells as receiver cells, we will only show expression in these
cells.

``` r
DotPlot(seuratObj %>% subset(idents = "CD8 T"), features = nichenet_output$top_receptors %>% rev(), split.by = "aggregate") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

You also can just show ‘bona fide’ ligand-receptor links that are
described in the literature and not predicted based on protein-protein
interactions:

``` r
nichenet_output$ligand_receptor_heatmap_bonafide
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
nichenet_output$ligand_receptor_matrix_bonafide
##            H2.M3     H2.T23        C3      Icam1     Tgfb1    Cxcl16      Il15
## Il2rb  0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.0000000 0.8714269
## Il2rg  0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.0000000 0.8587859
## Itgal  0.0000000 0.00000000 0.0000000 0.06542904 0.0000000 0.0000000 0.0000000
## Itgb2  0.0000000 0.00000000 0.2916032 0.06113009 0.0000000 0.0000000 0.0000000
## Tgfbr2 0.0000000 0.00000000 0.0000000 0.00000000 0.7665905 0.0000000 0.0000000
## Cxcr6  0.0000000 0.00000000 0.0000000 0.00000000 0.0000000 0.6598705 0.0000000
## Klrd1  0.8334165 0.05478448 0.0000000 0.00000000 0.0000000 0.0000000 0.0000000
nichenet_output$ligand_receptor_df_bonafide
## # A tibble: 9 x 3
##   ligand receptor weight
##   <chr>  <chr>     <dbl>
## 1 C3     Itgb2    0.292 
## 2 Cxcl16 Cxcr6    0.660 
## 3 H2-T23 Klrd1    0.0548
## 4 H2-M3  Klrd1    0.833 
## 5 Icam1  Itgal    0.0654
## 6 Icam1  Itgb2    0.0611
## 7 Il15   Il2rb    0.871 
## 8 Il15   Il2rg    0.859 
## 9 Tgfb1  Tgfbr2   0.767
```

If you are interested in checking which geneset (and background set of
genes) was used during the ligand activity analysis:

``` r
nichenet_output$geneset_oi
##   [1] "Irf7"          "Stat1"         "Ifit3"         "Ifit1"         "Bst2"         
##   [6] "B2m"           "Rnf213"        "Plac8"         "Isg15"         "Shisa5"       
##  [11] "Zbp1"          "Isg20"         "Samhd1"        "Usp18"         "H2-T23"       
##  [16] "Gbp2"          "Ifi203"        "Tmsb4x"        "Rsad2"         "Ly6e"         
##  [21] "Rtp4"          "Ifit2"         "Xaf1"          "Smchd1"        "Daxx"         
##  [26] "Alb"           "Samd9l"        "Actb"          "Parp9"         "Gbp4"         
##  [31] "Lgals3bp"      "Mx1"           "Gbp7"          "Cmpk2"         "Dtx3l"        
##  [36] "Slfn5"         "Oasl1"         "Herc6"         "Ifih1"         "Rpsa"         
##  [41] "P2ry13"        "Irgm2"         "Tapbp"         "Rps8"          "Stat2"        
##  [46] "Ifi44"         "Rpl8"          "Psmb8"         "Igfbp4"        "Ddx58"        
##  [51] "Rac2"          "Trafd1"        "Pml"           "Oas2"          "Psme1"        
##  [56] "Apoe"          "Basp1"         "Rps27a"        "Znfx1"         "Rpl13"        
##  [61] "Oas3"          "Nt5c3"         "Rnf114"        "Tap1"          "Rps28"        
##  [66] "Rplp0"         "Ddx60"         "Vim"           "Ifi35"         "Itm2b"        
##  [71] "Ctss"          "Pabpc1"        "Parp14"        "Hspa8"         "Tor3a"        
##  [76] "Rpl23"         "Tmbim6"        "Thy1"          "Ncoa7"         "Dhx58"        
##  [81] "Rps10"         "Rps19"         "Psmb9"         "Il2rg"         "Etnk1"        
##  [86] "Irf9"          "1600014C10Rik" "Parp12"        "Eif2ak2"       "Eef1b2"       
##  [91] "Eef2"          "Npc2"          "Rps2"          "Rps3"          "Sp110"        
##  [96] "Ube2l6"        "Nmi"           "Uba7"          "Psmb10"        "Cxcl10"       
## [101] "Rpl13a"        "Nhp2"          "Tbrg1"         "Usp25"         "Tor1aip2"     
## [106] "Adar"          "Gzma"          "Cd53"          "Hspa5"         "Cfl1"         
## [111] "Crip1"         "Slco3a1"       "Tlr7"          "Trim21"        "Rpl10"        
## [116] "Mycbp2"        "Rps16"         "Nlrc5"         "Rplp2"         "Acadl"        
## [121] "Trim12c"       "Rps4x"         "Irf1"          "Psma2"         "Nme2"         
## [126] "Zcchc11"       "Snord12"       "Phip"          "Ifitm3"        "Sp140"        
## [131] "Dusp2"         "Mrpl30"        "H2-M3"         "Gbp3"          "Dtx1"         
## [136] "Eef1g"         "Rbl1"          "Xpo1"          "Gm9844"        "Rpl35"        
## [141] "Rps26"         "Cxcr4"         "Eif3m"         "Treml2"        "Rpl35a"       
## [146] "Pdcd4"         "Arrb2"         "Ubc"           "Clic4"         "Rpl10a"       
## [151] "Lcp1"          "Cd274"         "Ddit4"         "Cnn2"          "Nampt"        
## [156] "Ascc3"         "Cd47"          "Snord49b"      "D17Wsu92e"     "Fam26f"       
## [161] "Hcst"          "Myh9"          "Rps27"         "Mov10"         "Arf4"         
## [166] "Arhgdib"       "Ppib"          "Trim25"        "Tspo"          "Id3"          
## [171] "Snord35a"      "Rnf8"          "Casp8"         "Ptpn7"         "Itk"          
## [176] "Cd69"          "Nop10"         "Anxa6"         "Hk1"           "Prkcb"        
## [181] "Iqgap1"        "Keap1"         "Rpl7"          "Parp10"
nichenet_output$background_expressed_genes %>% length()
## [1] 1487
```

### Rerun the NicheNet analysis with different sender cell definition

Instead of focusing on multiple sender cell types, it is possible that
you are only interested in doing the analyis for one sender cell type,
such as dendritic cells in this case.

``` r
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seuratObj, receiver = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "DC", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"

nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

Instead of focusing on one or multiple predefined sender cell types, it
is also possible that you want to consider all cell types present as
possible sender cell. This also includes the receiver cell type, making
that you can look at autocrine signaling as well.

``` r
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seuratObj, receiver = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "all", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"

nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

In some cases, it could be possible that you don’t have data of
potential sender cells. If you still want to predict possible upstream
ligands that could have been responsible for the observed differential
expression in your cell type, you can do this by following command. This
will consider all possible ligands in the NicheNet databases for which a
receptor is expressed by the receiver cell of interest.

``` r
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seuratObj, receiver = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "undefined", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"

nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

As you can see in this analysis result, many genes DE in CD8 T cells
after LCMV infection are strongly predicted type I interferon targets.
The presence of a type I interferon signature in the receiver cell type,
but the absence of expression of type I interferons in sender cell
types, might indicate that type I interferons are expressed by a
different, non-profiled cell type, or at a time point before sampling.
The latter could make sense, because there always is a time delay
between expression of a ligand-encoding gene and the effect of the
ligand on a target/receiver cell (i.e. expression of target genes).

### Run multiple NicheNet analyses on different receiver cell populations

In some cases, you might be interested in multiple target/receiver cell
populations. You can decide to run this for every cell type separately,
or in one line of code as demonstrated here (results are the same). As
example, we could have been interested in explaining DE between
steady-state and LCMV infection in both CD8 and CD4 T cells.

``` r
receiver_celltypes_oi = c("CD4 T", "CD8 T")
# receiver_celltypes_oi = seuratObj %>% Idents() %>% unique() # for all celltypes in the dataset: use only when this would make sense biologically

nichenet_output = receiver_celltypes_oi %>% lapply(nichenet_seuratobj_aggregate, seurat_obj = seuratObj, condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"

names(nichenet_output) = receiver_celltypes_oi
```

Check which ligands were top-ranked for both CD8T and CD4T and which
ligands were more cell-type specific

``` r
common_ligands = intersect(nichenet_output$`CD4 T`$top_ligands, nichenet_output$`CD8 T`$top_ligands)
print("common ligands are: ")
## [1] "common ligands are: "
print(common_ligands)
##  [1] "Ebi3"   "Il15"   "Crlf2"  "H2-M3"  "App"    "Ptprc"  "Icam1"  "Ccl5"   "Cxcl10" "Tgfb1" 
## [11] "Cxcl11" "Sema4d" "Cxcl9"  "H2-T23" "Cxcl16" "C3"     "Itgb1"

cd4_ligands = nichenet_output$`CD4 T`$top_ligands %>% setdiff(nichenet_output$`CD8 T`$top_ligands)
cd8_ligands = nichenet_output$`CD8 T`$top_ligands %>% setdiff(nichenet_output$`CD4 T`$top_ligands)

print("Ligands specifically regulating DE in CD4T: ")
## [1] "Ligands specifically regulating DE in CD4T: "
print(cd4_ligands)
## [1] "Cd274" "Hmgb1" "Cd28"

print("Ligands specifically regulating DE in CD8T: ")
## [1] "Ligands specifically regulating DE in CD8T: "
print(cd8_ligands)
## [1] "Adam17" "Anxa1"  "Sell"
```

## NicheNet analysis on Seurat object: explain differential expression between two cell populations

Previously, we demonstrated the use of a wrapper function for applying
NicheNet to explain differential expression between two conditions in
one cell type. However, also differential expression between two cell
populations might sometimes be (partially) caused by communication with
cells in the neighborhood. For example, differentiation from a
progenitor cell to the differentiated cell might be induced by niche
cells. A concrete example is discussed in this paper: [Stellate Cells,
Hepatocytes, and Endothelial Cells Imprint the Kupffer Cell Identity on
Monocytes Colonizing the Liver Macrophage
Niche](https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1).

Therefore, we will now also demonstrate the use of another Seurat
wrapper function that can be used in the case of explaining differential
expression between cell populations. But keep in mind that the
comparison that you make should be biologically relevant. It is possible
to use NicheNet to explain differential expression beween any two cell
populations in your dataset, but in most cases, differential expression
between cell populations will be a result of cell-intrinisc properties
(i.e. different cell types have a different gene expression profile) and
not of intercellular communication processes. In such a case, it does
not make any sense to use NicheNet.

For demonstration purposes, we will here first change the seuratObject
of the data described above, such that it can be used in this setting.

``` r
seuratObj@meta.data$celltype = paste(seuratObj@meta.data$celltype,seuratObj@meta.data$aggregate, sep = "_")

seuratObj@meta.data$celltype %>% table()
## .
##     B_LCMV       B_SS CD4 T_LCMV   CD4 T_SS CD8 T_LCMV   CD8 T_SS    DC_LCMV      DC_SS 
##        344         38       1961        601       1252        393         14          4 
##  Mono_LCMV    Mono_SS    NK_LCMV      NK_SS  Treg_LCMV    Treg_SS 
##         75         15         94         37        146         53

seuratObj = SetIdent(seuratObj,value = "celltype")
```

Now perform the NicheNet analysis to explain differential expression
between the ‘affected’ cell population ‘CD8 T cells after LCMV
infection’ and the reference cell population ‘CD8 T cells in
steady-state’ by ligands expressed by monocytes and DCs after LCMV
infection.

``` r
nichenet_output = nichenet_seuratobj_cluster_de(
  seurat_obj = seuratObj, 
  receiver_reference = "CD8 T_SS", receiver_affected = "CD8 T_LCMV", 
  sender = c("DC_LCMV","Mono_LCMV"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis between two receiver cell clusters"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
```

Check the top-ranked ligands and their target genes

``` r
nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Check the expression of the top-ranked ligands

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

It could be interested to check which top-ranked ligands are
differentially expressed in monocytes after LCMV infection

``` r
Mono_upregulated_ligands = FindMarkers(seuratObj, ident.1 = "Mono_LCMV", ident.2 = "Mono_SS") %>% rownames_to_column("gene") %>% filter(avg_log2FC > 0.25 & p_val_adj <= 0.05) %>% pull(gene) %>% intersect(nichenet_output$top_ligands)

print("Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-StSt and CD8T-LCMV are: ")
## [1] "Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-StSt and CD8T-LCMV are: "
print(Mono_upregulated_ligands)
## [1] "Cxcl10"
```

# Remarks

1.  Top-ranked ligands and target genes shown here differ from the
    predictions shown in the respective case study in the NicheNet paper
    because a different definition of expressed genes was used.
2.  Differential expression is here done via the classical Wilcoxon test
    used in Seurat to define marker genes of a cell cluster by comparing
    it to other clusters. This is not optimal if you would have repeated
    samples for your conditions. In such a case, we recommend to follow
    the vignette [Perform NicheNet analysis starting from a Seurat
    object: step-by-step
    analysis](seurat_steps.md):`vignette("seurat_steps", package="nichenetr")`
    and tweak the differential expression step there (and perform the
    analysis e.g. as discussed in <https://github.com/HelenaLC/muscat>).

# References

Bonnardel et al., 2019, Immunity 51, 1–17, [Stellate Cells, Hepatocytes,
and Endothelial Cells Imprint the Kupffer Cell Identity on Monocytes
Colonizing the Liver Macrophage
Niche](https://doi.org/10.1016/j.immuni.2019.08.017)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-medaglia_spatial_2017" class="csl-entry">

Medaglia, Chiara, Amir Giladi, Liat Stoler-Barak, Marco De Giovanni,
Tomer Meir Salame, Adi Biram, Eyal David, et al. 2017. “Spatial
Reconstruction of Immune Niches by Combining Photoactivatable Reporters
and <span class="nocase">scRNA</span>-Seq.” *Science*, December,
eaao4277. <https://doi.org/10.1126/science.aao4277>.

</div>

</div>
