Differential NicheNet
================
Robin Browaeys
2019-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/differential_nichenet.Rmd", output_format = "github_document")
-->

Remark: this is a beta version of a new extension of NicheNet, namely
Differential NicheNet. Short-term improvements will include scalability,
visualization and documentation of this vignette and the underlying
functions (end of May 2021).

The goal of Differential NicheNet is to predict ligand-receptors pairs
that are both differential expressed and active between different niches
of interest.

This vignette guides you in detail through all the steps of a
Differntial NicheNet analysis. As example expression data of interacting
cells, we will use data from Puram et al. to explore intercellular
communication in the tumor microenvironment in head and neck squamous
cell carcinoma (HNSCC) See (Puram et al. 2017). More specifically, we
will look at cell-cell communication differences between pEMT-high and
pEMT-low tumors (pEMT = partial epithelial-mesenschymal transition).

The used ligand-receptor network and ligand-target matrix can be
downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3260758.svg)](https://doi.org/10.5281/zenodo.3260758).
The Seurat object containing expression data of interacting cells in
HNSCC can also be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4675430.svg)](https://doi.org/10.5281/zenodo.4675430).

# 0. Read in the expression data of interest, and the NicheNet ligand-receptor network and ligand-target matrix

## Load in packages

``` r
library(circlize)
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(RColorBrewer)
library(tidyverse)
source("../R/differential_nichenet.R")
source("../R/differential_nichenet_plotting.R")
```

## Read in the expression data

In this case study, we want to study differences in cell-cell
communication patterns between pEMT-high and pEMT-low tumors. The meta
data columns that indicate the pEMT status of tumors are ‘pEMT,’ and the
cell type is indicated in the ‘celltype’ column.

``` r
seurat_obj = readRDS(url("https://zenodo.org/record/4675430/files/seurat_obj_hnscc.rds"))
DimPlot(seurat_obj, group.by = "celltype")
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
DimPlot(seurat_obj, group.by = "pEMT")
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

We will now also check the number of cells per cell type condition
combination

``` r
table(seurat_obj@meta.data$celltype, seurat_obj@meta.data$pEMT) # cell types vs conditions
##                
##                 High  Low
##   CAF            396  104
##   Endothelial    105   53
##   Malignant     1093  549
##   Myeloid         92    7
##   myofibroblast  382   61
##   T.cell         689    3
```

For the Differential NicheNet, we need to compare at least 2 niches or
conditions to each other. In this case, the 2 niches are the LCMV-niche
and the SS-niche. We will adapt the names of the cell types based on
their niche of origin.

``` r
seurat_obj@meta.data$celltype_aggregate = paste(seurat_obj@meta.data$celltype, seurat_obj@meta.data$pEMT,sep = "_")
DimPlot(seurat_obj, group.by = "celltype_aggregate")
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
seurat_obj@meta.data$celltype_aggregate %>% table() %>% sort(decreasing = TRUE)
## .
##     Malignant_High        T.cell_High      Malignant_Low           CAF_High myofibroblast_High   Endothelial_High            CAF_Low       Myeloid_High  myofibroblast_Low    Endothelial_Low        Myeloid_Low         T.cell_Low 
##               1093                689                549                396                382                105                104                 92                 61                 53                  7                  3
```

``` r
celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])
```

## Read in the NicheNet ligand-receptor network and ligand-target matrix

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

``` r
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

head(lr_network)
## # A tibble: 6 x 3
##   ligand receptor bonafide
##   <chr>  <chr>    <lgl>   
## 1 CXCL1  CXCR2    TRUE    
## 2 CXCL2  CXCR2    TRUE    
## 3 CXCL3  CXCR2    TRUE    
## 4 CXCL5  CXCR2    TRUE    
## 5 PPBP   CXCR2    TRUE    
## 6 CXCL6  CXCR2    TRUE
```

Note: if your data is of mouse origin: convert human gene symbols to
their one-to-one orthologs

``` r
organism = "human"
```

``` r
if(organism == "mouse"){
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
}
```

# 1. Define the niches/microenvironments of interest

Each niche should have at least one “sender/niche” cell population and
one “receiver/target” cell population (present in your expression data)

In this case study, we are interested to find differences in cell-cell
interactions to malignant cells between pEMT high and pEMT low tumors.
The receiver cell population in the pEMT-High niche is thus the
“Malignant\_High” cell type, and in the pEMT-Low niche this is
“Malignant\_Low.” The sender cell populations of interest are
myofibroblasts, Endothelial, CAF, T.cell, and Myeloid. Importantly, we
only include T.Cell and Myeloid in the pEMT-High niche, because there
are too few cells of these populations present in the pEMT-low niche.
Hereby, we demonstrate the possibility to include a condition-specific
cell type in the analysis - which is possible because we calculate DE
compared to all sender cells of the other niche, and not only to the
pEMT-low group of cells of the same cell type.

``` r
niches = list(
  "pEMT_High_niche" = list(
    "sender" = c("myofibroblast_High", "Endothelial_High", "CAF_High", "T.cell_High", "Myeloid_High"),
    "receiver" = c("Malignant_High")),
  "pEMT_Low_niche" = list(
    "sender" = c("myofibroblast_Low",  "Endothelial_Low", "CAF_Low"),
    "receiver" = c("Malignant_Low"))
  )
```

# 2. Calculate differential expression between the niches

In this step, we will determine DE between the different niches for both
senders and receivers to define the DE of L-R pairs.

### Calculate DE

The method to calculate the differential expression is here the standard
Seurat Wilcoxon test, but this can be replaced if wanted by the user
(only requirement: output tables `DE_sender_processed` and
`DE_receiver_processed` should be in the same format as shown here).

DE will be calculated for each pairwise sender (or receiver) cell type
comparision between the niches (so across niches, not within niche). In
our case study, this means that DE of myofibroblast\_High ligands will
be calculated by DE analysis of myofibroblast\_High vs
myofibroblast\_Low, myofibroblast\_High vs Endothelial\_Low, and
myofibroblast\_High vs CAF\_Low. We split the cells per cell type
instead of merging all cells from the other niche to avoid that the DE
analysis will be driven by the most abundant cell types.

Because of all these pairwise comparisons, this can take a long time to
run. Scalability will be improved in the near future

``` r
assay_oi = "SCT" # other possibilities: RNA,...
DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important
## [1] "Calculate Sender DE between: myofibroblast_High and myofibroblast_Low" "Calculate Sender DE between: myofibroblast_High and Endothelial_Low"   "Calculate Sender DE between: myofibroblast_High and CAF_Low"          
## [1] "Calculate Sender DE between: Endothelial_High and myofibroblast_Low" "Calculate Sender DE between: Endothelial_High and Endothelial_Low"   "Calculate Sender DE between: Endothelial_High and CAF_Low"          
## [1] "Calculate Sender DE between: CAF_High and myofibroblast_Low" "Calculate Sender DE between: CAF_High and Endothelial_Low"   "Calculate Sender DE between: CAF_High and CAF_Low"          
## [1] "Calculate Sender DE between: T.cell_High and myofibroblast_Low" "Calculate Sender DE between: T.cell_High and Endothelial_Low"   "Calculate Sender DE between: T.cell_High and CAF_Low"          
## [1] "Calculate Sender DE between: Myeloid_High and myofibroblast_Low" "Calculate Sender DE between: Myeloid_High and Endothelial_Low"   "Calculate Sender DE between: Myeloid_High and CAF_Low"          
## [1] "Calculate Sender DE between: myofibroblast_Low and myofibroblast_High" "Calculate Sender DE between: myofibroblast_Low and Endothelial_High"   "Calculate Sender DE between: myofibroblast_Low and CAF_High"           "Calculate Sender DE between: myofibroblast_Low and T.cell_High"       
## [5] "Calculate Sender DE between: myofibroblast_Low and Myeloid_High"      
## [1] "Calculate Sender DE between: Endothelial_Low and myofibroblast_High" "Calculate Sender DE between: Endothelial_Low and Endothelial_High"   "Calculate Sender DE between: Endothelial_Low and CAF_High"           "Calculate Sender DE between: Endothelial_Low and T.cell_High"       
## [5] "Calculate Sender DE between: Endothelial_Low and Myeloid_High"      
## [1] "Calculate Sender DE between: CAF_Low and myofibroblast_High" "Calculate Sender DE between: CAF_Low and Endothelial_High"   "Calculate Sender DE between: CAF_Low and CAF_High"           "Calculate Sender DE between: CAF_Low and T.cell_High"        "Calculate Sender DE between: CAF_Low and Myeloid_High"
DE_receiver = calculate_niche_de(seurat_obj = seurat_obj, niches = niches, type = "receiver", assay_oi = assay_oi) # all genes important, because also the goal to infer DE targets
## [1] "Calculate Receiver DE between: Malignant_High and Malignant_Low"
## [1] "Calculate Receiver DE between: Malignant_Low and Malignant_High"
```

### Process DE results:

``` r
expression_pct = 0.10
DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
```

### Combine sender-receiver DE based on L-R pairs:

As mentioned above, DE of ligands from one sender cell type is
determined be calculating DE between that cell type, and all the sender
cell types of the other niche. To summarize the DE of ligands of that
cell type we have several options: we could take the average LFC, but
also the minimum LFC compared to the other niche. We recommend using the
minimum LFC, because this is the strongest specificity measure of ligand
expression, because a high min LFC means that a ligand is more strongly
expressed in the cell type of niche 1 compared to all cell types of
niche 2 (in contrast to a high average LFC, which could occur if one or
more cell types in niche 2 also strongly express that ligand).

``` r
specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
```

# 3. Optional: Calculate differential expression between the different spatial regions

To improve the cell-cell interaction predictions, you can consider
spatial information if possible and applicable. Spatial information can
come from microscopy data, or from spatial transcriptomics data such as
Visium.

There are several ways to incorporate spatial information in the
Differential NicheNet pipeline. First, you can only consider cell types
as belonging to the same niche if they are in the same spatial location.
Another way is including spatial differential expression of
ligand-receptor pairs within one cell type in the prioritization
framework.

For example: We have a cell type X, located in regions A and B, and we
want to study cell-cell communication in region A. We first add only
celltypeX of regionA in the niche definition, and then calculate DE
between celltypeX-regionA and celltypeX-regionB to give higher
prioritization weight to regionA-specific ligands.

In this case study, our region of interest is the tumor leading edge,
since Puram et al defined this region as important regarding the pEMT
process. Puram et al also defined CAFs as the fibroblasts that are close
to leading edge, whereas the other fibroblasts (myofibroblasts) were not
preferentially located in the tumor leading edge. We can thus now
prioritize fibroblast ligands further by looking at ligands that are DE
between leading-edge fibroblasts (=CAFs) and non-leading-edge
fibroblasts (myofibroblasts).

We do this as follows, by first defining a ‘spatial info’ dataframe. If
no spatial information in your data: set the following two parameters to
FALSE, and make a mock ‘spatial\_info’ data frame.

``` r
include_spatial_info_sender = TRUE # if not spatial info to include: put this to false
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true
```

``` r
spatial_info = tibble(celltype_region_oi = c("CAF_High"), celltype_other_region = c("myofibroblast_High")) %>% mutate(niche =  "pEMT_High_niche", celltype_type = "sender") %>% 
  bind_rows(tibble(celltype_region_oi = c("CAF_Low"), celltype_other_region = c("myofibroblast_Low")) %>% mutate(niche =  "pEMT_Low_niche", celltype_type = "sender"))
specificity_score_spatial = "lfc"
```

``` r
# this is how this should be defined if you don't have spatial info
# mock spatial info
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
} 
```

``` r
if(include_spatial_info_sender == TRUE){
  sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj, spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral zonation score for sender celltypes in which the zonation is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)

  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_zonation = scale_quantile_adapted(ligand_score_zonation))

} else {
  # # add a neutral zonation score for all sender celltypes (for none of them, zonation is relevant in this case)
  sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_zonation = scale_quantile_adapted(ligand_score_zonation))  

}
## [1] "Calculate Spatial DE between: CAF_High and myofibroblast_High"
## [1] "Calculate Spatial DE between: CAF_Low and myofibroblast_Low"
```

``` r
if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj, spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral zonation score for receiver celltypes in which the zonation is not known / not of importance
  receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)

  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_zonation = scale_quantile_adapted(receptor_score_zonation))

} else {
    # # add a neutral zonation score for all receiver celltypes (for none of them, zonation is relevant in this case)
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_zonation = scale_quantile_adapted(receptor_score_zonation))
}
```

# 4. Calculate ligand activities and infer active ligand-target links

In this step, we will predict ligand activities of each ligand for each
of the receiver cell types across the different niches. This is similar
to the ligand activity analysis done in the normal NicheNet pipeline.

To calculate ligand activities, we first need to define a geneset of
interest for each niche. In this case study, the geneset of interest for
the pEMT-high niche are the genes upregulated in pEMT-high tumors
compared to pEMT-low tumors, and vice versa.

Note that you can also define these geneset of interest in a different
way! (eg pathway-based geneset etc)

Ligand-target links are inferred in the same way as described in the
basic NicheNet vignettes.

``` r
lfc_cutoff = 0.15
top_n_target = 250

specificity_score_targets = "min_lfc"

DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
  
background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  
## too many genes - damn... - maybe change this...

# Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
# If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
##  [1] "ANXA8L2"       "PRKCDBP"       "IL8"           "PTRF"          "SEPP1"         "C1orf186"      "CCDC109B"      "C10orf54"      "LEPREL1"       "ZNF812"        "LOC645638"     "LOC401397"     "LINC00162"     "DFNA5"         "PLK1S1"        "ZMYM6NB"       "C19orf10"      "CTSL1"         "SQRDL"         "LOC375295"    
## [21] "WBP5"          "LOC100505633"  "AIM1"          "C1orf63"       "LOC100507463"  "GPR115"        "VIMP"          "SEP15"         "C1orf172"      "NAPRT1"        "LHFP"          "KRT16P1"       "C7orf10"       "PTPLA"         "GRAMD3"        "CPSF3L"        "MESDC2"        "C10orf10"      "KIAA1609"      "CCDC53"       
## [41] "TXLNG2P"       "NGFRAP1"       "ERO1L"         "FAM134A"       "LSMD1"         "TCEB2"         "B3GALTL"       "HN1L"          "LOC550643"     "KIAA0922"      "GLT25D1"       "FAM127A"       "C1orf151-NBL1" "SEPW1"         "GPR126"        "LOC100505806"  "LINC00478"     "TCEB1"         "GRAMD2"        "GNB2L1"       
## [61] "KIRREL"
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
##   [1] "LOC344887"    "AGPAT9"       "C1orf110"     "KIAA1467"     "LOC100292680" "EPT1"         "CT45A4"       "LOC654433"    "UPK3BL"       "LINC00340"    "LOC100128338" "FAM60A"       "CCDC144C"     "LOC401109"    "LOC286467"    "LEPREL4"      "LOC731275"    "LOC642236"    "LINC00516"    "LOC101101776" "SC5DL"       
##  [22] "PVRL4"        "LOC100130093" "LINC00338"    "LOC100132891" "PPAP2C"       "C6orf1"       "C2orf47"      "WHSC1L1"      "LOC100289019" "SETD8"        "KDM5B-AS1"    "SPG20"        "CXCR7"        "LOC100216479" "LOC100505761" "MGC57346"     "LPHN3"        "CENPC1"       "C11orf93"     "C14orf169"    "LOC100506060"
##  [43] "FLJ31485"     "LOC440905"    "MLF1IP"       "TMEM194A"     "RRP7B"        "REXO1L1"      "LOC100129269" "KIAA1715"     "CTAGE5"       "LOC202781"    "LOC100506714" "LOC401164"    "UTS2D"        "LOC146880"    "KIAA1804"     "C5orf55"      "C21orf119"    "PRUNE"        "LRRC16A"      "LOC339240"    "FLJ35024"    
##  [64] "C5orf28"      "LOC100505876" "MGC21881"     "LOC100133985" "PPAPDC2"      "FRG1B"        "CECR5"        "LOC100129361" "CCBL1"        "PTPLAD1"      "MST4"         "LOC550112"    "LOC389791"    "CCDC90A"      "KIAA0195"     "LOC100506469" "LOC100133161" "LOC646719"    "LOC728819"    "BRE"          "LOC284581"   
##  [85] "LOC441081"    "LOC728377"    "LOC100134229" "C3orf65"      "SMEK2"        "KIAA1737"     "C17orf70"     "PLEKHM1P"     "LOC338758"    "PCNXL2"       "LOC91948"     "C17orf89"     "LOC100505783" "SMCR7L"       "C8orf4"       "GPR56"        "ATHL1"        "LOC339535"    "PPAPDC1B"     "DAK"          "LOC100507173"
## [106] "CRHR1-IT1"    "PPAP2B"       "ADCK4"        "KIAA0146"     "GYLTL1B"      "LOC100272216" "LOC400027"    "WHSC1"        "LOC100130855" "C7orf55"      "C19orf40"     "ADCK3"        "C9orf142"     "SGOL1"        "LOC90834"     "PTPLAD2"      "KIAA1967"     "LOC100132352" "LOC100630918" "ADRBK2"       "LINC00263"   
## [127] "FAM64A"       "LOC401074"    "FAM179B"      "RP1-177G6.2"  "METTL21D"     "ERO1LB"       "FLJ45445"     "NADKD1"       "LOC100506233" "LOC100652772" "FAM175A"      "LINC00630"    "C11orf82"     "SETD5-AS1"    "SGK196"       "FLJ14186"     "CCDC104"      "FAM63A"       "NARG2"        "MTERFD1"      "CCDC74B-AS1" 
## [148] "LOC286186"    "WDR67"        "C12orf52"     "FLJ30403"     "KIAA2018"     "GCN1L1"       "FLJ43681"     "LOC152217"    "FONG"         "C18orf8"      "ALG1L9P"      "GTDC2"        "LOC100507217" "NBPF24"       "WBSCR27"      "C14orf1"      "LOC284889"    "KIAA0317"     "FAM65A"       "PMS2L2"       "LUST"        
## [169] "C15orf52"     "FAM195A"      "LOC399744"    "PYCRL"        "LOC338799"    "LOC100506190" "C9orf91"      "FLJ45340"     "LOC349196"    "LOC100128881" "TOMM70A"      "ALS2CR8"      "LDOC1L"       "HDGFRP3"      "ZNF767"       "LOC728558"    "LOC283693"    "LEPREL2"      "QTRTD1"       "SELM"         "C6orf25"     
## [190] "C1orf86"      "HNRPLL"       "LOC145820"    "LOC100289341" "C17orf85"     "C3orf72"      "C14orf64"     "C9orf9"       "LOC100506394"
  
niche_geneset_list = list(
  "pEMT_High_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "pEMT_Low_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )
  
ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
## [1] "Calculate Ligand activities for: Malignant_High"
## [1] "Calculate Ligand activities for: Malignant_Low"
```

# 5. Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)

In this step, we will calculate average (scaled) expression, and
fraction of expression, of ligands, receptors, and target genes across
all cell types of interest. Now this is here demonstrated via the
DotPlot function of Seurat, but this can also be done via other ways of
course.

``` r
features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target)
  
dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl = dotplot$data %>% as_tibble()
exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
```

``` r
exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
```

# 6. Expression fraction and receptor

In this step, we will score ligand-receptor interactions based on
expression strength of the receptor, in such a way that we give higher
scores to the most strongly expressed receptor of a certain ligand, in a
certain celltype. This will not effect the rank of individual ligands
later on, but will help in prioritizing the most important receptors per
ligand (next to other factors regarding the receptor - see later).

``` r
exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 
```

# 7. Prioritization of ligand-receptor and ligand-target links

``` r
output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
         ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
```

explanation per factor - TO DO TO DO TO DO

``` r
prioritizing_weights = c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 2,
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "scaled_avg_score_ligand_receptor" = 0,
                         "scaled_ligand_score_zonation" = 0, #2
                         "scaled_receptor_score_zonation" = 0,
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_activity" = 0,"scaled_activity_normalized" = 1,
                         "ligand_fraction" = 1, "receptor_fraction" = 1, 
                         "bona_fide" = 1)
```

``` r
prioritization_tables = get_prioritization_tables(output, prioritizing_weights)

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
## # A tibble: 10 x 37
##    niche   receiver  sender  ligand_receptor ligand  receptor bonafide ligand_score ligand_significa~ ligand_present ligand_expression ligand_expressio~ ligand_fraction ligand_score_zo~ receptor_score receptor_signif~ receptor_present receptor_expres~ receptor_express~ receptor_fracti~ receptor_score_z~ ligand_scaled_recep~
##    <chr>   <chr>     <chr>   <chr>           <chr>   <chr>    <lgl>           <dbl>             <dbl>          <dbl>             <dbl>             <dbl>           <dbl>            <dbl>          <dbl>            <dbl>            <dbl>            <dbl>             <dbl>            <dbl>             <dbl>                <dbl>
##  1 pEMT_H~ Malignan~ T.cell~ PTPRC--MET      PTPRC   MET      FALSE            3.22                 1              1              9.32              2.5            0.939                0        0.463                  1                1            1.02              1.76             0.526                 0                0.900
##  2 pEMT_H~ Malignan~ T.cell~ PTPRC--EGFR     PTPRC   EGFR     FALSE            3.22                 1              1              9.32              2.5            0.939                0        0.454                  1                1            1.21              1.21             0.581                 0                0.950
##  3 pEMT_H~ Malignan~ T.cell~ PTPRC--CD44     PTPRC   CD44     FALSE            3.22                 1              1              9.32              2.5            0.939                0        0.104                  1                1            3.27              0.327            0.905                 0                1    
##  4 pEMT_H~ Malignan~ T.cell~ PTPRC--ERBB2    PTPRC   ERBB2    FALSE            3.22                 1              1              9.32              2.5            0.939                0       -0.0286                 0                1            0.629             1.16             0.345                 0                0.799
##  5 pEMT_H~ Malignan~ T.cell~ PTPRC--IFNAR1   PTPRC   IFNAR1   FALSE            3.22                 1              1              9.32              2.5            0.939                0        0.248                  0                1            0.765            -0.658            0.446                 0                0.850
##  6 pEMT_H~ Malignan~ Myeloi~ SERPINA1--LRP1  SERPIN~ LRP1     TRUE             2.52                 1              1              4.83              2.5            0.761                0       -0.159                  1                1            0.312            -0.526            0.234                 0                1    
##  7 pEMT_H~ Malignan~ T.cell~ PTPRC--INSR     PTPRC   INSR     FALSE            3.22                 1              1              9.32              2.5            0.939                0       -0.0722                 0                1            0.582            -0.719            0.252                 0                0.749
##  8 pEMT_H~ Malignan~ T.cell~ PTPRC--LEPR     PTPRC   LEPR     FALSE            3.22                 1              1              9.32              2.5            0.939                0       -0.307                  1                1            0.279            -0.833            0.191                 0                0.699
##  9 pEMT_H~ Malignan~ T.cell~ GZMB--IGF2R     GZMB    IGF2R    TRUE             2.46                 1              1              4.57              2.5            0.374                0       -0.00836                0                1            0.318            -0.432            0.224                 0                0.7  
## 10 pEMT_H~ Malignan~ Myeloi~ ITGB2--ICAM1    ITGB2   ICAM1    TRUE             2.48                 1              1              5.26              2.13           0.848                0        0.462                  1                1            0.761            -0.714            0.256                 0                0.739
## # ... with 15 more variables: avg_score_ligand_receptor <dbl>, activity <dbl>, activity_normalized <dbl>, scaled_ligand_score <dbl>, scaled_ligand_expression_scaled <dbl>, scaled_receptor_score <dbl>, scaled_receptor_expression_scaled <dbl>, scaled_avg_score_ligand_receptor <dbl>, scaled_ligand_score_zonation <dbl>,
## #   scaled_receptor_score_zonation <dbl>, scaled_ligand_fraction_adapted <dbl>, scaled_receptor_fraction_adapted <dbl>, scaled_activity <dbl>, scaled_activity_normalized <dbl>, prioritization_score <dbl>
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
## # A tibble: 10 x 20
##    niche           receiver       sender      ligand_receptor ligand receptor bonafide target  target_score target_significant target_present target_expression target_expression_scaled target_fraction ligand_target_weight activity activity_normalized scaled_activity scaled_activity_normalized prioritization_score
##    <chr>           <chr>          <chr>       <chr>           <chr>  <chr>    <lgl>    <chr>          <dbl>              <dbl>          <dbl>             <dbl>                    <dbl>           <dbl>                <dbl>    <dbl>               <dbl>           <dbl>                      <dbl>                <dbl>
##  1 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    ACTB           0.185                  1              1             9.02                    -0.250           1                  0.00133    0.179                1.13           0.971                      0.873                0.925
##  2 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    ADM            0.288                  1              1             3.23                     1.38            0.693              0.00113    0.179                1.13           0.971                      0.873                0.925
##  3 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    AHNAK          0.236                  1              1             1.88                    -0.716           0.763              0.00136    0.179                1.13           0.971                      0.873                0.925
##  4 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    AJUBA          0.470                  1              1             0.751                    2.5             0.370              0.00117    0.179                1.13           0.971                      0.873                0.925
##  5 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    ATF3           0.510                  1              1             2.19                    -0.421           0.559              0.00125    0.179                1.13           0.971                      0.873                0.925
##  6 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    BHLHE40        0.384                  1              1             2.11                    -0.579           0.653              0.00143    0.179                1.13           0.971                      0.873                0.925
##  7 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    CCND1          0.314                  1              1             1.47                     1.54            0.570              0.00164    0.179                1.13           0.971                      0.873                0.925
##  8 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    CDKN1A         0.462                  1              1             2.68                    -0.777           0.763              0.00172    0.179                1.13           0.971                      0.873                0.925
##  9 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    CITED2         0.333                  1              1             0.392                   -0.891           0.202              0.00124    0.179                1.13           0.971                      0.873                0.925
## 10 pEMT_High_niche Malignant_High T.cell_High PTPRC--MET      PTPRC  MET      FALSE    CLIC1          0.239                  1              1             5.39                     0.729           0.962              0.00133    0.179                1.13           0.971                      0.873                0.925

prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
## # A tibble: 10 x 37
##    niche   receiver  sender   ligand_receptor ligand receptor bonafide ligand_score ligand_significa~ ligand_present ligand_expression ligand_expressio~ ligand_fraction ligand_score_zo~ receptor_score receptor_signif~ receptor_present receptor_expres~ receptor_express~ receptor_fracti~ receptor_score_z~ ligand_scaled_recep~
##    <chr>   <chr>     <chr>    <chr>           <chr>  <chr>    <lgl>           <dbl>             <dbl>          <dbl>             <dbl>             <dbl>           <dbl>            <dbl>          <dbl>            <dbl>            <dbl>            <dbl>             <dbl>            <dbl>             <dbl>                <dbl>
##  1 pEMT_L~ Malignan~ Endothe~ F8--LRP1        F8     LRP1     TRUE            0.952               1                1             2.17               2.5            0.528            0             0.159                 1                1            0.464           -0.404             0.346                 0                1    
##  2 pEMT_L~ Malignan~ CAF_Low  SLIT2--SDC1     SLIT2  SDC1     TRUE            0.494               1                1             0.846              2.5            0.288            0.749         0.0498                0                1            3.47             1.92              0.918                 0                1    
##  3 pEMT_L~ Malignan~ CAF_Low  SLIT2--GPC1     SLIT2  GPC1     TRUE            0.494               1                1             0.846              2.5            0.288            0.749         0.377                 1                1            1.49             2.21              0.719                 0                0.875
##  4 pEMT_L~ Malignan~ CAF_Low  RSPO3--LGR6     RSPO3  LGR6     TRUE            0.557               0.8              1             1.27               2.5            0.240            1.18          0.435                 1                1            0.384            2.39              0.164                 0                0.75 
##  5 pEMT_L~ Malignan~ CAF_Low  FGF10--FGFR2    FGF10  FGFR2    TRUE            0.385               0.8              1             1.07               2.46           0.25             0.934         0.154                 1                1            0.585            2.22              0.328                 0                1    
##  6 pEMT_L~ Malignan~ Endothe~ PLAT--LRP1      PLAT   LRP1     TRUE            0.913               1                1             2.70               2.12           0.509            0             0.159                 1                1            0.464           -0.404             0.346                 0                1    
##  7 pEMT_L~ Malignan~ Endothe~ IL33--IL1RAP    IL33   IL1RAP   FALSE           1.34                1                1             2.75               2.5            0.585            0            -0.582                 1                1            0.341            0.178             0.195                 0                1    
##  8 pEMT_L~ Malignan~ CAF_Low  COMP--SDC1      COMP   SDC1     TRUE            0.290               0.8              1             1.27               2.31           0.202            1.14          0.0498                0                1            3.47             1.92              0.918                 0                1    
##  9 pEMT_L~ Malignan~ CAF_Low  FGF14--FGFR2    FGF14  FGFR2    TRUE            0.200               0.4              1             0.221              2.32           0.106            0.219         0.154                 1                1            0.585            2.22              0.328                 0                1    
## 10 pEMT_L~ Malignan~ CAF_Low  SEMA3C--NRP2    SEMA3C NRP2     TRUE            0.652               1                1             1.73               2.07           0.423            1.21         -0.0634                0                1            0.565           -0.0958            0.293                 0                1    
## # ... with 15 more variables: avg_score_ligand_receptor <dbl>, activity <dbl>, activity_normalized <dbl>, scaled_ligand_score <dbl>, scaled_ligand_expression_scaled <dbl>, scaled_receptor_score <dbl>, scaled_receptor_expression_scaled <dbl>, scaled_avg_score_ligand_receptor <dbl>, scaled_ligand_score_zonation <dbl>,
## #   scaled_receptor_score_zonation <dbl>, scaled_ligand_fraction_adapted <dbl>, scaled_receptor_fraction_adapted <dbl>, scaled_activity <dbl>, scaled_activity_normalized <dbl>, prioritization_score <dbl>
prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
## # A tibble: 10 x 20
##    niche          receiver      sender          ligand_receptor ligand receptor bonafide target  target_score target_significant target_present target_expression target_expression_scaled target_fraction ligand_target_weight activity activity_normalized scaled_activity scaled_activity_normalized prioritization_score
##    <chr>          <chr>         <chr>           <chr>           <chr>  <chr>    <lgl>    <chr>          <dbl>              <dbl>          <dbl>             <dbl>                    <dbl>           <dbl>                <dbl>    <dbl>               <dbl>           <dbl>                      <dbl>                <dbl>
##  1 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     ASPSCR1        0.302                  1              1             0.800                   2.10             0.393             0.000855   0.0476               0.560           0.305                      0.820                0.813
##  2 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     AXIN2          0.195                  1              1             0.228                  -0.115            0.107             0.000749   0.0476               0.560           0.305                      0.820                0.813
##  3 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     BCL6           0.299                  1              1             1.64                   -0.342            0.689             0.000841   0.0476               0.560           0.305                      0.820                0.813
##  4 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     BDNF           0.275                  1              1             0.222                   2.5              0.144             0.000732   0.0476               0.560           0.305                      0.820                0.813
##  5 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     BRCA1          0.204                  1              1             0.317                   2.5              0.184             0.000909   0.0476               0.560           0.305                      0.820                0.813
##  6 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     CBX5           0.360                  1              1             1.41                    0.462            0.741             0.000770   0.0476               0.560           0.305                      0.820                0.813
##  7 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     CBX8           0.187                  1              1             0.239                   2.5              0.149             0.000789   0.0476               0.560           0.305                      0.820                0.813
##  8 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     CEBPA          0.195                  1              1             0.204                   1.01             0.131             0.000770   0.0476               0.560           0.305                      0.820                0.813
##  9 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     CTNNB1         0.304                  1              1             3.97                   -0.481            0.982             0.000746   0.0476               0.560           0.305                      0.820                0.813
## 10 pEMT_Low_niche Malignant_Low Endothelial_Low F8--LRP1        F8     LRP1     TRUE     DDIT3          0.382                  1              1             1.89                   -0.0166           0.599             0.000775   0.0476               0.560           0.305                      0.820                0.813
```

# 8. Visualization of the Differential NicheNet output

## Differential expression of ligand and expression

Visualization: minimum LFC compared to other niches

``` r
receiver_oi = "Malignant_High"
prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% group_by(sender) %>% top_n(15, prioritization_score) # top15 per sender
prioritized_tbl_oi2 = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% ungroup() %>% top_n(20, prioritization_score) # top10 in total
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

filtered_ligands = prioritized_tbl_oi %>% pull(ligand) %>% unique()

prioritized_tbl_lr_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(1, prioritization_score) %>% ungroup() 
prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup()  %>% filter(prioritization_score > prioritized_tbl_lr_oi$prioritization_score %>% min())

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

Show the spatialDE as additional information

``` r
lfc_plot = make_ligand_receptor_lfc_zonation_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

## Ligand expression, activity and target genes

Active target gene inference - cf Default NicheNet

Now: visualization of ligand activity and ligand-target links

``` r
prioritized_tbl_oi =  prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% group_by(sender) %>% top_n(5, prioritization_score) 
prioritized_tbl_oi2 =  prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% distinct() %>% group_by(ligand, receiver) %>% top_n(1, prioritization_score) %>% filter(receiver == receiver_oi) %>% ungroup() %>% top_n(20, prioritization_score) 
prioritized_tbl_oi = prioritized_tbl_oi %>% bind_rows(prioritized_tbl_oi2) %>% distinct()

exprs_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, plot_legend = FALSE, heights = NULL, widths = NULL)
exprs_plot$combined_plot
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

## Circos plot of prioritized ligand-receptor pairs

``` r
prioritized_tbl_all = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == receiver_oi) %>% group_by(ligand) %>% top_n(1, prioritization_score)  %>% ungroup() %>% top_n(20, prioritization_score)
prioritized_tbl_all_extended = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == receiver_oi) %>% filter(ligand %in% prioritized_tbl_all$ligand) %>% filter(prioritization_score >= min(prioritized_tbl_all$prioritization_score)) %>% group_by(ligand) %>% top_n(3, prioritization_score) %>% ungroup()

prioritized_tbl_oi = prioritized_tbl_all_extended %>% distinct()
prioritized_tbl_oi = prioritized_tbl_oi %>% inner_join(prioritization_tables$prioritization_tbl_ligand_receptor %>% distinct(sender, receiver, niche) )  %>% arrange(sender, prioritization_score) 

colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% length(), name = 'Spectral') %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver, cutoff = 0, scale = FALSE, transparency = NULL, circos_type = "normal", border = TRUE)
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->![](differential_nichenet_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

``` r
circos_output$p_circos
```

![](differential_nichenet_files/figure-gfm/unnamed-chunk-29-3.png)<!-- -->

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-puram_single-cell_2017" class="csl-entry">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
