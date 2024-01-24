Perform NicheNet analysis starting from a Seurat object
================
Robin Browaeys
2023-10-02

<!-- github markdown built using 
rmarkdown::render("vignettes/seurat_wrapper.Rmd", output_format = "github_document")
-->

In this vignette, you can learn how to perform a basic NicheNet analysis
on a Seurat (v3-v5) object containing single-cell expression data.
**Assuming you have captured the changes in gene expression resulting
from your cell-cell communication (CCC) process of interest,** a
NicheNet analysis can help you to generate hypotheses about the CCC
process. Specifically, NicheNet can predict 1) which ligands from the
microenvironment or cell population(s) (“sender/niche”) are most likely
to affect target gene expression in an interacting cell population
(“receiver/target”) and 2) which specific target genes are affected by
which of these predicted ligands.

The wrapper function we will show consists of the same different steps
that are discussed in detail in the main NicheNet vignette [Perform
NicheNet analysis starting from a Seurat object: step-by-step
analysis](seurat_steps.md).Please make sure you understand the different
steps described in this vignette before performing a real NicheNet
analysis on your data. We generally recommend the step-by-step analysis
as it allows users to adapt specific steps of the pipeline to make them
more appropriate for their data.

To perform a NicheNet analysis, three features are extracted from the
input data: the potential ligands, the gene set of interest, and the
background gene set. This vignette will extract each feature as
described in this flowchart:

<img src="figure2.svg" style="width:50.0%" />

As example expression data of interacting cells, we will use mouse
NICHE-seq data to explore intercellular communication in the T cell area
in the inguinal lymph node before and 72 hours after lymphocytic
choriomeningitis virus (LCMV) infection (Medaglia et al. 2017). We will
focus on CD8 T cells as the receiver population, and as this dataset
contains two conditions (before and after LCMV infection), the
differentially expressed genes between these two conditions in CD8 T
cells will be used as our gene set of interest. We will then prioritize
which ligands from the microenvironment (sender-agnostic approach) and
from specific immune cell populations like monocytes, dendritic cells,
NK cells, B cells, and CD4 T cells (sender-focused approach) can
regulate and induce these observed gene expression changes.

The [ligand-target matrix](https://doi.org/10.5281/zenodo.7074290) and
the [Seurat object of the processed NICHE-seq single-cell
data](https://doi.org/10.5281/zenodo.3531889) can be downloaded from
Zenodo.

# Prepare NicheNet analysis

### Load packages

``` r
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
```

If you would use and load other packages, we recommend to load these 3
packages after the others.

### Read in NicheNet’s networks

The ligand-target prior model, ligand-receptor network, and weighted
integrated networks are needed for this vignette. The ligand-target
prior model is a matrix describing the potential that a ligand may
regulate a target gene, and it is used to run the ligand activity
analysis. The ligand-receptor network contains information on potential
ligand-receptor bindings, and it is used to identify potential ligands.
Finally, the weighted ligand-receptor network contains weights
representing the potential that a ligand will bind to a receptor, and it
is used for visualization.

``` r
organism <- "mouse"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

}

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
## # A tibble: 6 × 2
##   from          to   
##   <chr>         <chr>
## 1 2300002M23Rik Ddr1 
## 2 2610528A11Rik Gpr15
## 3 9530003J23Rik Itgal
## 4 a             Atrn 
## 5 a             F11r 
## 6 a             Mc1r
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##               2300002M23Rik 2610528A11Rik 9530003J23Rik            a          A2m
## 0610005C13Rik  0.000000e+00  0.000000e+00  1.311297e-05 0.000000e+00 1.390053e-05
## 0610009B22Rik  0.000000e+00  0.000000e+00  1.269301e-05 0.000000e+00 1.345536e-05
## 0610009L18Rik  8.872902e-05  4.977197e-05  2.581909e-04 7.570125e-05 9.802264e-05
## 0610010F05Rik  2.194046e-03  1.111556e-03  3.142374e-03 1.631658e-03 2.585820e-03
## 0610010K14Rik  2.271606e-03  9.360769e-04  3.546140e-03 1.697713e-03 2.632082e-03

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
## # A tibble: 6 × 3
##   from          to     weight
##   <chr>         <chr>   <dbl>
## 1 0610010F05Rik App    0.110 
## 2 0610010F05Rik Cat    0.0673
## 3 0610010F05Rik H1f2   0.0660
## 4 0610010F05Rik Lrrc49 0.0829
## 5 0610010F05Rik Nicn1  0.0864
## 6 0610010F05Rik Srpk1  0.123
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
## # A tibble: 6 × 3
##   from          to            weight
##   <chr>         <chr>          <dbl>
## 1 0610010K14Rik 0610010K14Rik 0.121 
## 2 0610010K14Rik 2510039O18Rik 0.121 
## 3 0610010K14Rik 2610021A01Rik 0.0256
## 4 0610010K14Rik 9130401M01Rik 0.0263
## 5 0610010K14Rik Alg1          0.127 
## 6 0610010K14Rik Alox12        0.128
```

### Read in the expression data of interacting cells

We processed and aggregated the original dataset by using the Seurat
alignment pipeline. As we created this object using Seurat v3, it has to
be updated with `UpdateSeuratObject`. Note that genes should be named by
their official mouse/human gene symbol. If your expression data has the
older gene symbols, you may want to use our alias conversion function to
avoid the loss of gene names.

``` r
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))

# For newer Seurat versions, you may need to run the following
seuratObj <- UpdateSeuratObject(seuratObj)

# Convert gene names
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

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
(DCs) and inflammatory monocytes.

``` r
# Note that the number of cells of some cell types is very low and should preferably be higher for a real application
seuratObj@meta.data$celltype %>% table() 
## .
##     B CD4 T CD8 T    DC  Mono    NK  Treg 
##   382  2562  1645    18    90   131   199
DimPlot(seuratObj, reduction = "tsne")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Visualize the data to see to which condition cells belong. The metadata
column that denotes the condition (steady-state or after LCMV infection)
is here called ‘aggregate’.

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
expressed by the microenvironment (sender-agnostic) and immune cells in
the T cell area of the lymph node (sender-focused) are most likely to
have induced the differential expression in CD8 T cells after LCMV
infection. In contrary to NicheNet v1 where we only used the
“sender-focused” approach, we now recommend users to run both the
“sender-agnostic” approach and “sender-focused” approach. These
approaches only affect the list of potential ligands that are considered
for prioritization. As described in the flowchart above, we do not
define any sender populations in the ‘sender agnostic’ approach but
consider all ligands for which its cognate receptor is expressed in the
receiver population. The sender-focused approach will then filter the
list of ligands to ones where the ligands are expressed in the sender
cell population(s).

As described in the main vignette, the pipeline of a basic NicheNet
analysis consist of the following steps: \* 1. Define a set of potential
ligands for both the sender-agnostic and sender-focused approach \* 2.
Define the gene set of interest: these are the genes in the
“receiver/target” cell population that are potentially affected by
ligands expressed by interacting cells (e.g. genes differentially
expressed upon cell-cell interaction) \* 3. Define the background genes
\* 4. Perform NicheNet ligand activity analysis: rank the potential
ligands based on the presence of their target genes in the gene set of
interest (compared to the background set of genes) \* 5. Infer target
genes and receptors of top-ranked ligands

All these steps are contained in one of three wrapper functions:
`nichenet_seuratobj_aggregate`, `nichenet_seuratobj_cluster_de` and
`nichenet_seuratobj_aggregate_cluster_de`. These functions differ on how
the gene set of interest is calculated, as follows:

| **Function**                            | **Gene set of interest**                           | **Background genes**                   |
|-----------------------------------------|----------------------------------------------------|----------------------------------------|
| nichenet_seuratobj_aggregate            | DE between two conditions of the same cell type    | All expressed genes in the cell type   |
| nichenet_seuratobj_cluster_de           | DE between two cell types                          | All expressed genes in both cell types |
| nichenet_seuratobj_aggregate_cluster_de | DE between two cell types from specific conditions | All expressed genes in both cell types |

**Note:** Cell types should be the identities of the seurat object
(check using `table(Idents(seuratObj))`)

## `nichenet_seuratobj_aggregate`: explain differential expression between two conditions

For the sender-agnostic approach the sender is set to ‘undefined’. The
receiver cell population is the ‘CD8 T’ cell population, and the gene
set of interest are the genes differentially expressed in CD8 T cells
after LCMV infection. Thus, the condition of interest is ‘LCMV’, whereas
the reference/steady-state condition is ‘SS’. The column containing
condition information is ‘aggregate’. The method to calculate
differential expression is the standard Seurat Wilcoxon test. To use
other methods, users will have to go through step-by-step analysis. The
number of top-ranked ligands that are further used to predict active
target genes and construct an active ligand-receptor network is 30
(`top_n_ligands`). The number of target genes to consider per ligand
when performing the target gene inference is 200 (`top_n_targets`). We
only retain ligands and receptors that are expressed in at least a
predefined fraction of cells in one cluster (`expression_pct`, default:
10%).

``` r
nichenet_output_agnostic <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  sender = "undefined",
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
```

For the sender-focused approach, simply provide one or more sender
populations:

``` r
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
```

**Note:** It is also possible that you want to consider all cell types
present as possible sender cell by defining `sender = "all"`. This also
includes the receiver cell type, making that you can look at autocrine
signaling as well.

### Interpret the NicheNet analysis output

We will investigate the output of the sender-focused approach.

``` r
names(nichenet_output)
##  [1] "ligand_activities"                      "top_ligands"                            "top_targets"                           
##  [4] "top_receptors"                          "ligand_target_matrix"                   "ligand_target_heatmap"                 
##  [7] "ligand_target_df"                       "ligand_expression_dotplot"              "ligand_differential_expression_heatmap"
## [10] "ligand_activity_target_heatmap"         "ligand_receptor_matrix"                 "ligand_receptor_heatmap"               
## [13] "ligand_receptor_df"                     "geneset_oi"                             "background_expressed_genes"
```

#### Ligand activity analysis results

To see the ranking of ligands based on the predicted ligand activity:

``` r
nichenet_output$ligand_activities
## # A tibble: 127 × 6
##    test_ligand auroc  aupr aupr_corrected pearson  rank
##    <chr>       <dbl> <dbl>          <dbl>   <dbl> <dbl>
##  1 Il27        0.682 0.391         0.316    0.445     1
##  2 Ebi3        0.666 0.264         0.189    0.256     2
##  3 Tnf         0.671 0.205         0.131    0.249     3
##  4 Ptprc       0.660 0.198         0.124    0.168     4
##  5 H2-Eb1      0.656 0.195         0.120    0.182     5
##  6 Vsig10      0.649 0.194         0.119    0.170     6
##  7 H2-M3       0.632 0.192         0.118    0.185     7
##  8 Clcf1       0.637 0.175         0.101    0.162     8
##  9 H2-M2       0.634 0.174         0.0989   0.146    11
## 10 H2-T-ps     0.634 0.174         0.0989   0.146    11
## # ℹ 117 more rows
```

Ligands are ranked based on the area under the precision-recall curve
(AUPR) between a ligand’s target predictions and the observed
transcriptional response. Although other metrics like the AUROC and
pearson correlation coefficient are also computed, we demonstrated in
our validation study that the AUPRwas the most informative measure to
define ligand activity (this was the Pearson correlation for v1). The
vignette on how we performed the validation can be found at [Evaluation
of NicheNet’s ligand-target predictions](model_evaluation.md).

To get a list of the top 30 ligands:

``` r
nichenet_output$top_ligands
##  [1] "Il27"    "Ebi3"    "Tnf"     "Ptprc"   "H2-Eb1"  "Vsig10"  "H2-M3"   "Clcf1"   "H2-M2"   "H2-T-ps" "H2-T10"  "H2-T22"  "H2-T24"  "H2-T23" 
## [15] "H2-K1"   "H2-Q4"   "H2-Q6"   "H2-Q7"   "H2-D1"   "H2-Oa"   "Il18bp"  "Sirpa"   "Cd48"    "App"     "Ccl22"   "Siglech" "Ccl5"    "Siglec1"
## [29] "Cd320"   "Adam17"
```

Below we will show visualizations that are in the output object. In some
cases (including this one), not all top ligands that are present in
`top_ligands` will be shown in the plot. The left-out ligands are
ligands that don’t have target genes with high enough regulatory
potential scores, and therefore did not survive the used cutoffs (in the
functions `get_weighted_ligand_target_links` and
`prepare_ligand_target_visualization` that are run internally). To
include them, you can increase the number of target genes considered or
be less stringent in the used cutoffs (`top_n_targets` and
`cutoff_visualization` , respectively). In this case, CCl22 (ranked
25th) is missing from the plots.

To see which sender cell population expresses which of the top-ranked
ligands:

``` r
nichenet_output$ligand_expression_dotplot
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

As you can see, most op the top-ranked ligands seem to be mainly
expressed by dendritic cells and monocytes.

It could also be interesting to see whether some of these ligands are
differentially expressed after LCMV infection.

``` r
nichenet_output$ligand_differential_expression_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Although this ligand differential expression is not used for
prioritization and ranking of the ligands (the ranking is only
determined based on enrichment of target genes among DE genes in the
receiver, CD8T cells), most of the top-ranked ligands also seem to be
upregulated themselves in monocytes after viral infection. This is nice
additional “evidence” that these ligands might indeed be important.

#### Inferred active ligand-target links

NicheNet also infers active target genes of these top-ranked ligands,
best visualized with the following heatmap showing which top-ranked
ligands are predicted to have regulated the expression of which
differentially expressed genes:

``` r
nichenet_output$ligand_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

This is a normal ggplot object that can be adapted accordingly. For
example if you want to change the color code to blue instead of purple,
change the axis ticks of the legend, and change the axis labels of the
heatmap, you can do the following:

``` r
nichenet_output$ligand_target_heatmap +
  scale_fill_gradient2(low = "whitesmoke",high = "royalblue") +
  xlab("anti-LCMV response genes in CD8 T cells") + ylab("Prioritized immmune cell ligands")
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

If you want, you can also extract the ligand-target links and their
regulatory potential scores in matrix or data frame format (e.g. for
visualization in other ways or output to a csv file).

``` r
nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
##         Adar         B2m        Bst2 Calhm6      Cd274      Cxcl10
## Adam17     0 0.000000000 0.008167279      0 0.06549177 0.011094196
## Cd320      0 0.000000000 0.000000000      0 0.00000000 0.000000000
## Siglec1    0 0.000000000 0.000000000      0 0.00000000 0.000000000
## Ccl5       0 0.000000000 0.000000000      0 0.00000000 0.008424993
## Siglech    0 0.008857572 0.011974948      0 0.01257584 0.008780173
## App        0 0.000000000 0.000000000      0 0.04432138 0.000000000
## Cd48       0 0.000000000 0.000000000      0 0.00000000 0.000000000
## Sirpa      0 0.000000000 0.000000000      0 0.00000000 0.007796006
## Il18bp     0 0.000000000 0.000000000      0 0.00000000 0.007808540
## H2.Oa      0 0.000000000 0.000000000      0 0.00000000 0.008143571
```

``` r
nichenet_output$ligand_target_df # weight column = regulatory potential
## # A tibble: 656 × 3
##    ligand target weight
##    <chr>  <chr>   <dbl>
##  1 Il27   Adar    0.163
##  2 Il27   B2m     0.170
##  3 Il27   Bst2    0.111
##  4 Il27   Calhm6  0.129
##  5 Il27   Cd274   0.111
##  6 Il27   Cxcl10  0.178
##  7 Il27   Cxcr4   0.178
##  8 Il27   Ddx58   0.227
##  9 Il27   Ddx60   0.160
## 10 Il27   Dtx3l   0.150
## # ℹ 646 more rows
```

To get a list of the top-predicted target genes of the 30 top-ranked
ligands:

``` r
nichenet_output$top_targets
##  [1] "Adar"     "B2m"      "Bst2"     "Calhm6"   "Cd274"    "Cxcl10"   "Cxcr4"    "Ddx58"    "Ddx60"    "Dtx3l"    "Eif2ak2"  "Gbp2"     "Gbp3"    
## [14] "Gbp7"     "H2-D1"    "H2-K1"    "H2-M3"    "H2-Q6"    "H2-Q7"    "H2-T10"   "H2-T22"   "H2-T23"   "Ifi203"   "Ifi206"   "Ifi208"   "Ifi209"  
## [27] "Ifi213"   "Ifi35"    "Ifi44"    "Ifih1"    "Ifit1bl1" "Ifit2"    "Ifit3"    "Ifit3b"   "Ifitm3"   "Irf1"     "Irf7"     "Irf9"     "Lgals3bp"
## [40] "Ly6e"     "Mndal"    "Mx1"      "Mx2"      "Nampt"    "Nlrc5"    "Nmi"      "Oas2"     "Oas3"     "Parp12"   "Parp14"   "Parp9"    "Pml"     
## [53] "Psmb8"    "Psmb9"    "Psme1"    "Psme2b"   "Rnf213"   "Samhd1"   "Sp110"    "Stat1"    "Stat2"    "Tap1"     "Tapbp"    "Tnfsf10"  "Trafd1"  
## [66] "Ube2l6"   "Xaf1"     "Ddit4"    "Dhx58"    "Gzmb"     "Isg15"    "Lcp1"     "Oas1a"    "Oas1g"    "Rsad2"    "Zbp1"     "Cd47"     "Ctss"    
## [79] "Trim21"   "Cd69"     "H3f3b"    "Id3"      "Vim"      "Isg20"    "Oasl1"    "Hspa5"    "Ifit1"    "Nt5c3"    "Usp18"    "Basp1"    "Plac8"   
## [92] "Sp100"    "Sp140"    "Ubc"
```

You can visualize the expression of these target genes as well (only the
top 50 are shown here). Because we only focus on CD8 T cells as receiver
cells, we will only show expression in these cells. To emphasize that
these target genes are differentially expressed, we split cells up in
steady-state cells and cells after response to LCMV infection.

``` r
DotPlot(seuratObj %>% subset(idents = "CD8 T"), features = nichenet_output$top_targets[1:50] %>%
          rev(), split.by = "aggregate") + coord_flip()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
VlnPlot(seuratObj %>% subset(idents = "CD8 T"), features = c("Ptprc", "H2-M3", "Cxcl10"), split.by = "aggregate", pt.size = 0, combine = TRUE)
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

The display the combined plot of ligand activities, expression,
differential expression and target genes of ligands:

``` r
nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

**Important: the above figure can be considered as one of the most
important summary figures of the NicheNet analysis. Here you can see
which ligand-receptor pairs have both high differential expression and
ligand activity (=target gene enrichment). These are very interesting
predictions as key regulators of your intercellular communication
process of interest!**

#### Inferred ligand-receptor interactions for top-ranked ligands

NicheNet also infers the receiver cell receptors of these top-ranked
ligands. You can run following command for a heatmap visualization of
the ligand-receptor links:

``` r
nichenet_output$ligand_receptor_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

If you want, you can also extract the ligand-receptor links and their
interaction confidence scores in matrix or data frame format (e.g. for
visualization in other ways or output to a csv file).

``` r
nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]
##        H2.T24 H2.T23 H2.T22 H2.T10 H2.T.ps H2.Q7
## Il6ra       0      0      0      0       0     0
## Itgb1       0      0      0      0       0     0
## Notch1      0      0      0      0       0     0
## Ptprc       0      0      0      0       0     0
## Spn         0      0      0      0       0     0
## Cd47        0      0      0      0       0     0
## Cd69        0      0      0      0       0     0
## Ccr7        0      0      0      0       0     0
## Dpp4        0      0      0      0       0     0
## Cd247       0      0      0      0       0     0
```

``` r
nichenet_output$ligand_receptor_df # weight column accords to number of data sources that document this interaction
## # A tibble: 54 × 3
##    ligand receptor weight
##    <chr>  <chr>     <dbl>
##  1 Adam17 Il6ra     0.447
##  2 Adam17 Itgb1     0.454
##  3 Adam17 Notch1    1.05 
##  4 App    Cd74      0.670
##  5 App    Sorl1     0.922
##  6 Ccl22  Ccr7      0.679
##  7 Ccl22  Dpp4      0.717
##  8 Ccl5   Cxcr3     0.848
##  9 Cd320  Jaml      0.507
## 10 Cd320  Tmem167   0.432
## # ℹ 44 more rows
```

To get a list of the receptors of the 30 top-ranked ligands:

``` r
nichenet_output$top_receptors
##  [1] "Il6ra"    "Itgb1"    "Notch1"   "Cd74"     "Sorl1"    "Ccr7"     "Dpp4"     "Cxcr3"    "Jaml"     "Tmem167"  "Cd2"      "Il6st"    "Il27ra"  
## [14] "Cd8a"     "Klrd1"    "Cd4"      "Cd247"    "Cd47"     "Ptprc"    "Spn"      "Cd69"     "Tnfrsf1b"
```

You can visualize the expression of these as well. Because we only focus
on CD8 T cells as receiver cells, we will only show expression in these
cells.

``` r
DotPlot(seuratObj %>% subset(idents = "CD8 T"),
        features = nichenet_output$top_receptors %>% rev(), split.by = "aggregate") +
  coord_flip()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

If you are interested in checking which geneset (and background set of
genes) was used during the ligand activity analysis:

``` r
nichenet_output$geneset_oi
##   [1] "Ifi27l2b"      "Irf7"          "Ly6a"          "Stat1"         "Ly6c2"         "Ifit3"         "Ifit1"         "Ly6c1"         "Bst2"         
##  [10] "B2m"           "Rnf213"        "Ifit1bl1"      "Plac8"         "Slfn1"         "Ifi209"        "Isg15"         "Igtp"          "Ifi206"       
##  [19] "Shisa5"        "Ms4a4c"        "H2-K1"         "Zbp1"          "Oasl2"         "Isg20"         "Samhd1"        "Ifi208"        "Ms4a6b"       
##  [28] "Trim30a"       "Usp18"         "Mndal"         "H2-T23"        "Slfn8"         "Gbp2"          "Ifi203"        "Iigp1"         "Tmsb4x"       
##  [37] "H2-T22"        "Rsad2"         "Ly6e"          "Rtp4"          "Ifit3b"        "Zfas1"         "Ifit2"         "Phf11b"        "Xaf1"         
##  [46] "Smchd1"        "Daxx"          "Alb"           "Samd9l"        "Actb"          "Parp9"         "Gbp4"          "Lgals3bp"      "Mx1"          
##  [55] "Ifi213"        "Irgm1"         "2410006H16Rik" "Gbp7"          "Cmpk2"         "Dtx3l"         "Slfn5"         "H2-D1"         "Oasl1"        
##  [64] "Herc6"         "Ifih1"         "Rpsa"          "P2ry13"        "Apoa2"         "Irgm2"         "Tapbp"         "Rps8"          "Stat2"        
##  [73] "Ifi44"         "Phf11c"        "Rpl8"          "Psmb8"         "Gm12250"       "Igfbp4"        "Rplp2-ps1"     "Ddx58"         "Rac2"         
##  [82] "Trafd1"        "Sp100"         "Gbp9"          "Pml"           "Oas2"          "Slfn2"         "Psme1"         "Apoe"          "Gas5"         
##  [91] "H2-Q7"         "Basp1"         "Ms4a4b"        "Rps27a"        "Cd52"          "Znfx1"         "Rpl13"         "Ahsg"          "Oas3"         
## [100] "Nt5c3"         "Rnf114"        "Tap1"          "Rps28"         "Oas1a"         "Rplp0"         "Ddx60"         "Vim"           "Gbp6"         
## [109] "Ifi35"         "Itm2b"         "Ctss"          "Tgtp1"         "Trf"           "Pabpc1"        "H2-Q6"         "Parp14"        "Hspa8"        
## [118] "Tor3a"         "Rpl23"         "Mx2"           "Tmbim6"        "Thy1"          "Ncoa7"         "Dhx58"         "Rps10"         "Rps19"        
## [127] "Psmb9"         "Il2rg"         "Etnk1"         "Irf9"          "Rps3a1"        "Gbp10"         "1600014C10Rik" "Parp12"        "Trim30d"      
## [136] "Eif2ak2"       "Eef1b2"        "Eef2"          "Ncf2"          "Npc2"          "Rps2"          "Rps3"          "Sp110"         "Ube2l6"       
## [145] "Nmi"           "Uba7"          "Psmb10"        "Cxcl10"        "Rpl13a"        "Trim30c"       "Nhp2"          "Tbrg1"         "Jaml"         
## [154] "Usp25"         "Tor1aip2"      "Adar"          "Gzma"          "Gm2000"        "Rps18-ps5"     "Cd53"          "Phf11"         "Hspa5"        
## [163] "Cfl1"          "Crip1"         "Slco3a1"       "Tlr7"          "Trim21"        "Gbp8"          "Rpl10"         "Mycbp2"        "Rps16"        
## [172] "Nlrc5"         "Rplp2"         "Acadl"         "Trim12c"       "Rps4x"         "Irf1"          "Psma2"         "Nme2"          "Tut4"         
## [181] "Apobec3"       "Snord12"       "Phip"          "Gzmb"          "Ifitm3"        "Sp140"         "Dusp2"         "Mrpl30"        "Malat1"       
## [190] "H2-M3"         "Gbp3"          "Tmsb10"        "Dtx1"          "Tmem184b"      "Eef1g"         "Rbl1"          "Epb41l4aos"    "Xpo1"         
## [199] "Rgcc"          "Gm9844"        "Rpl35"         "Rps26"         "Il18bp"        "Sdc3"          "Cxcr4"         "Eif3m"         "Treml2"       
## [208] "Rpl35a"        "Lgals8"        "Pdcd4"         "Arrb2"         "Ubc"           "Clic4"         "H2-T10"        "Rpl10a"        "Lcp1"         
## [217] "Cd274"         "Ddit4"         "Cnn2"          "Nampt"         "Ascc3"         "Ms4a6d"        "Cd47"          "Ogfrl1"        "Snord49b"     
## [226] "Ilrun"         "Calhm6"        "Psme2b"        "Hcst"          "Myh9"          "Rps27"         "Mov10"         "Gm15772"       "Arf4"         
## [235] "Arhgdib"       "Ppib"          "Ubb"           "Trim25"        "Tspo"          "Id3"           "Snord35a"      "Zup1"          "Oas1g"        
## [244] "Ms4a6c"        "Rnf8"          "Casp8"         "Tnfsf10"       "Ptpn7"         "Itk"           "Rps27rt"       "Cd69"          "H3f3b"        
## [253] "Nop10"         "Anxa6"         "Hk1"           "Prkcb"         "Iqgap1"        "Keap1"         "Rpl7"          "Parp10"
nichenet_output$background_expressed_genes %>% length()
## [1] 3476
```

### Results of the sender-agnostic approach

``` r
# There is no log-fold change or expression plot because we did not define cell types
nichenet_output_agnostic$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

As you can see in this analysis result, many genes DE in CD8 T cells
after LCMV infection are strongly predicted type I interferon targets.
The presence of a type I interferon signature in the receiver cell type,
but the absence of expression of type I interferons in sender cell
types, might indicate that type I interferons are expressed by a
different, non-profiled cell type, or at a time point before sampling.
The latter could make sense, because there always is a time delay
between expression of a ligand-encoding gene and the effect of the
ligand on a target/receiver cell (i.e. expression of target genes).

#### Running multiple NicheNet analyses on different receiver cell populations

In some cases, you might be interested in multiple target/receiver cell
populations. You can decide to run this for every cell type separately,
or in one line of code as demonstrated here (results are the same). As
example, we could have been interested in explaining DE between
steady-state and LCMV infection in both CD8 and CD4 T cells.

``` r
# To run with  all celltypes in the dataset (only when this would make sense biologically!)
# receiver_celltypes_oi <- seuratObj %>% Idents() %>% unique()

receiver_celltypes_oi = c("CD4 T", "CD8 T")

nichenet_output <- receiver_celltypes_oi %>% lapply(nichenet_seuratobj_aggregate,
                                                    seurat_obj = seuratObj,
                                                    condition_colname = "aggregate",
                                                    condition_oi = "LCMV",
                                                    condition_reference = "SS",
                                                    sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
                                                    ligand_target_matrix = ligand_target_matrix,
                                                    lr_network = lr_network,
                                                    weighted_networks = weighted_networks)
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

names(nichenet_output) <- receiver_celltypes_oi
```

Check which ligands were top-ranked for both CD8T and CD4T and which
ligands were more cell-type specific

``` r
common_ligands <- intersect(nichenet_output$`CD4 T`$top_ligands, nichenet_output$`CD8 T`$top_ligands)
print("Common ligands:")
## [1] "Common ligands:"
print(common_ligands)
##  [1] "Ebi3"   "Ptprc"  "H2-M3"  "H2-M2"  "H2-T10" "H2-T22" "H2-T23" "Sirpa"  "H2-K1"  "H2-Q4"  "H2-Q6"  "H2-Q7"  "H2-D1"  "Ccl22"  "Cd48"   "App"   
## [17] "Tgfb1"  "Selplg" "Icam1"  "Btla"   "Cd72"   "B2m"    "Hp"     "Itgb2"

cd4_ligands <- nichenet_output$`CD4 T`$top_ligands %>% setdiff(nichenet_output$`CD8 T`$top_ligands)
cd8_ligands <- nichenet_output$`CD8 T`$top_ligands %>% setdiff(nichenet_output$`CD4 T`$top_ligands)

print("Ligands specifically regulating DE in CD4T:")
## [1] "Ligands specifically regulating DE in CD4T:"
print(cd4_ligands)
## [1] "H2-Eb1"  "H2-Oa"   "Il16"    "Fn1"     "H2-DMb1" "H2-DMb2"

print("Ligands specifically regulating DE in CD8T:")
## [1] "Ligands specifically regulating DE in CD8T:"
print(cd8_ligands)
## [1] "Cxcl10" "Adam17" "Cxcl11" "Tgm2"   "Cxcl9"  "Vcan"
```

## `nichenet_seuratobj_cluster_de`: explain differential expression between two cell types

Unlike the case above where we applied NicheNet to explain differential
expression between two conditions in one cell type, here we try to
explain differential expression between two cell populations. DE between
cell populations are sometimes (partially) caused by communication with
cells in the neighborhood, e.g., the differentiation from a progenitor
cell to a differentiated cell might be induced by niche cells. A
concrete example is discussed in the paper by Bonnardel et al. (2019):
[Stellate Cells, Hepatocytes, and Endothelial Cells Imprint the Kupffer
Cell Identity on Monocytes Colonizing the Liver Macrophage
Niche](https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1).

However, keep in mind that the comparison that you make should be
biologically relevant. as in most cases, differential expression between
cell populations will be a result of cell-intrinsic properties
(i.e. different cell types have a different gene expression profile) and
not of an intercellular communication processes. In such a case, it does
not make any sense to use NicheNet.

For demonstration purposes, we will change the Seurat object of the same
dataset such that it can be used in this setting.

``` r
seuratObj <- SetIdent(seuratObj, value = paste(seuratObj$celltype, seuratObj$aggregate, sep = "_"))
Idents(seuratObj) %>% table()
## .
##   CD8 T_SS   CD4 T_SS    Treg_SS       B_SS      NK_SS    Mono_SS      DC_SS CD8 T_LCMV CD4 T_LCMV     B_LCMV  Treg_LCMV    NK_LCMV  Mono_LCMV 
##        393        601         53         38         37         15          4       1252       1961        344        146         94         75 
##    DC_LCMV 
##         14
```

Now perform the NicheNet analysis to explain differential expression
between the ‘affected’ cell population ‘CD8 T cells after LCMV
infection’ and the reference cell population ‘CD8 T cells in
steady-state’ by ligands expressed by monocytes and DCs after LCMV
infection.

``` r
nichenet_output <- nichenet_seuratobj_cluster_de(
  seurat_obj = seuratObj, 
  receiver_reference = "CD8 T_SS",
  receiver_affected = "CD8 T_LCMV", 
  sender = c("DC_LCMV", "Mono_LCMV"), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis between two receiver cell clusters"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
```

Check the top-ranked ligands and their target genes:

``` r
nichenet_output$ligand_activity_target_heatmap
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Check the expression of the top-ranked ligands:

``` r
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") +
  RotatedAxis()
```

![](seurat_wrapper_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

It could be interesting to check which top-ranked ligands are
differentially expressed in monocytes after LCMV infection:

``` r
Mono_upregulated_ligands <- FindMarkers(seuratObj, ident.1 = "Mono_LCMV", ident.2 = "Mono_SS") %>% 
  rownames_to_column("gene") %>% filter(avg_log2FC > 0.25 & p_val_adj <= 0.05) %>%
  pull(gene) %>% intersect(nichenet_output$top_ligands)

print("Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-SS and CD8T-LCMV are: ")
## [1] "Monocyte ligands upregulated after LCMV infection and explaining DE between CD8T-SS and CD8T-LCMV are: "
print(Mono_upregulated_ligands)
## [1] "B2m"    "H2-D1"  "Cxcl10"
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
    object: step-by-step analysis](seurat_steps.md) and tweak the
    differential expression step there (and perform the analysis e.g.,
    as discussed in <https://github.com/HelenaLC/muscat>).

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
##  [8] ggplot2_3.4.4      tidyverse_1.3.1    SeuratObject_5.0.1 Seurat_4.4.0       nichenetr_2.0.4    testthat_3.2.1    
## 
## loaded via a namespace (and not attached):
##   [1] IRanges_2.34.1              R.methodsS3_1.8.1           progress_1.2.3              ParamHelpers_1.14.1         nnet_7.3-19                
##   [6] goftest_1.2-3               vctrs_0.6.5                 spatstat.random_3.2-2       digest_0.6.33               png_0.1-8                  
##  [11] shape_1.4.6                 proxy_0.4-27                OmnipathR_3.9.6             ggrepel_0.9.4               deldir_2.0-2               
##  [16] parallelly_1.36.0           MASS_7.3-60                 reprex_2.0.1                reshape2_1.4.4              httpuv_1.6.13              
##  [21] foreach_1.5.2               BiocGenerics_0.46.0         withr_2.5.2                 ggrastr_1.0.2               xfun_0.41                  
##  [26] ggpubr_0.6.0                ellipsis_0.3.2              survival_3.5-7              memoise_2.0.1               ggbeeswarm_0.7.2           
##  [31] zoo_1.8-12                  GlobalOptions_0.1.2         pbapply_1.7-2               R.oo_1.24.0                 prettyunits_1.2.0          
##  [36] Formula_1.2-5               promises_1.2.1              httr_1.4.7                  rstatix_0.7.2               globals_0.16.2             
##  [41] fitdistrplus_1.1-11         rstudioapi_0.15.0           miniUI_0.1.1.1              generics_0.1.3              base64enc_0.1-3            
##  [46] dir.expiry_1.8.0            curl_5.2.0                  S4Vectors_0.38.1            zlibbioc_1.46.0             ScaledMatrix_1.8.1         
##  [51] polyclip_1.10-6             BBmisc_1.13                 randomForest_4.7-1.1        GenomeInfoDbData_1.2.10     xtable_1.8-4               
##  [56] desc_1.4.3                  doParallel_1.0.17           evaluate_0.23               S4Arrays_1.2.0              hms_1.1.3                  
##  [61] GenomicRanges_1.52.0        irlba_2.3.5.1               colorspace_2.1-0            filelock_1.0.2              visNetwork_2.1.2           
##  [66] ROCR_1.0-11                 smoof_1.6.0.3               reticulate_1.34.0           readxl_1.4.3                spatstat.data_3.0-3        
##  [71] magrittr_2.0.3              lmtest_0.9-40               later_1.3.2                 mlrMBO_1.1.5.1              lattice_0.21-9             
##  [76] spatstat.geom_3.2-7         future.apply_1.11.0         lhs_1.1.6                   scuttle_1.10.2              scattermore_1.2            
##  [81] shadowtext_0.1.2            cowplot_1.1.2               matrixStats_1.2.0           RcppAnnoy_0.0.21            class_7.3-22               
##  [86] Hmisc_5.1-0                 pillar_1.9.0                nlme_3.1-163                emoa_0.5-0.2                iterators_1.0.14           
##  [91] beachmat_2.16.0             caTools_1.18.2              compiler_4.3.2              stringi_1.7.6               gower_1.0.1                
##  [96] tensor_1.5                  SummarizedExperiment_1.30.2 lubridate_1.9.3             devtools_2.4.3              plyr_1.8.9                 
## [101] crayon_1.5.2                abind_1.4-5                 parallelMap_1.5.1           locfit_1.5-9.8              haven_2.4.3                
## [106] sp_2.1-2                    modelr_0.1.8                fastmatch_1.1-4             codetools_0.2-19            recipes_1.0.7              
## [111] BiocSingular_1.16.0         e1071_1.7-14                GetoptLong_1.0.5            plotly_4.10.0               mime_0.12                  
## [116] splines_4.3.2               circlize_0.4.15             Rcpp_1.0.11                 basilisk_1.12.1             sparseMatrixStats_1.12.2   
## [121] dbplyr_2.1.1                cellranger_1.1.0            knitr_1.45                  utf8_1.2.4                  clue_0.3-64                
## [126] fs_1.6.3                    listenv_0.9.0               checkmate_2.3.1             DelayedMatrixStats_1.22.5   logger_0.2.2               
## [131] pkgbuild_1.4.3              ggsignif_0.6.4              Matrix_1.6-4                statmod_1.5.0               tzdb_0.4.0                 
## [136] tweenr_2.0.2                pkgconfig_2.0.3             DiceKriging_1.6.0           tools_4.3.2                 cachem_1.0.8               
## [141] viridisLite_0.4.2           rvest_1.0.2                 DBI_1.1.3                   fastmap_1.1.1               rmarkdown_2.11             
## [146] scales_1.3.0                grid_4.3.2                  usethis_2.2.2               ica_1.0-3                   liana_0.1.12               
## [151] broom_0.7.12                patchwork_1.1.3             dotCall64_1.1-1             carData_3.0-5               RANN_2.6.1                 
## [156] rpart_4.1.21                farver_2.1.1                yaml_2.3.8                  MatrixGenerics_1.12.3       DiagrammeR_1.0.10          
## [161] foreign_0.8-85              cli_3.6.2                   stats4_4.3.2                leiden_0.3.9                lifecycle_1.0.4            
## [166] caret_6.0-94                uwot_0.1.16                 Biobase_2.60.0              mlr_2.19.1                  bluster_1.10.0             
## [171] lava_1.7.3                  sessioninfo_1.2.2           backports_1.4.1             BiocParallel_1.34.2         timechange_0.2.0           
## [176] gtable_0.3.4                rjson_0.2.21                ggridges_0.5.5              progressr_0.14.0            parallel_4.3.2             
## [181] pROC_1.18.5                 limma_3.56.2                edgeR_3.42.4                jsonlite_1.8.8              bitops_1.0-7               
## [186] assertthat_0.2.1            brio_1.1.4                  Rtsne_0.17                  BiocNeighbors_1.18.0        spatstat.utils_3.0-4       
## [191] highr_0.10                  metapod_1.8.0               dqrng_0.3.2                 R.utils_2.11.0              timeDate_4032.109          
## [196] lazyeval_0.2.2              shiny_1.7.1                 htmltools_0.5.7             sctransform_0.4.0           rappdirs_0.3.3             
## [201] basilisk.utils_1.12.1       glue_1.6.2                  spam_2.10-0                 XVector_0.40.0              RCurl_1.98-1.12            
## [206] rprojroot_2.0.4             scran_1.28.2                gridExtra_2.3               igraph_1.2.11               R6_2.5.1                   
## [211] SingleCellExperiment_1.22.0 fdrtool_1.2.17              labeling_0.4.3              cluster_2.1.4               pkgload_1.3.3              
## [216] GenomeInfoDb_1.36.1         ipred_0.9-14                vipor_0.4.5                 DelayedArray_0.26.7         tidyselect_1.2.0           
## [221] htmlTable_2.4.1             ggforce_0.4.1               xml2_1.3.6                  car_3.1-2                   future_1.33.0              
## [226] ModelMetrics_1.2.2.2        rsvd_1.0.5                  munsell_0.5.0               KernSmooth_2.23-22          data.table_1.14.10         
## [231] htmlwidgets_1.6.2           ComplexHeatmap_2.16.0       RColorBrewer_1.1-3          rlang_1.1.2                 spatstat.sparse_3.0-3      
## [236] spatstat.explore_3.2-1      remotes_2.4.2               ggnewscale_0.4.9            fansi_1.0.6                 hardhat_1.3.0              
## [241] beeswarm_0.4.0              prodlim_2023.08.28
```

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
