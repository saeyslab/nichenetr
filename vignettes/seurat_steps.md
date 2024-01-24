Perform NicheNet analysis starting from a Seurat object: step-by-step
analysis
================
Robin Browaeys & Chananchida Sang-aram
2023-10-02

<!-- github markdown built using 
rmarkdown::render("vignettes/seurat_steps.Rmd", output_format = "github_document")
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

Please make sure you understand the different steps described in this
vignette before performing a real NicheNet analysis on your data. There
are also wrapper functions that perform the same steps as in this
vignette in [Perform NicheNet analysis starting from a Seurat
object](seurat_wrapper.md). However, in that case users will not be able
to adapt specific steps of the pipeline to make them more appropriate
for their data.

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

### Read in the expression data of interacting cells

We processed and aggregated the original dataset by using the Seurat
alignment pipeline. As we created this object using Seurat v3, it has to
be updated with `UpdateSeuratObject`. Note that genes should be named by
their official mouse/human gene symbol.

``` r
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))

seuratObj@meta.data %>% head()
##         nGene nUMI orig.ident aggregate res.0.6 celltype nCount_RNA nFeature_RNA
## W380370   880 1611      LN_SS        SS       1    CD8 T       1607          876
## W380372   541  891      LN_SS        SS       0    CD4 T        885          536
## W380374   742 1229      LN_SS        SS       0    CD4 T       1223          737
## W380378   847 1546      LN_SS        SS       1    CD8 T       1537          838
## W380379   839 1606      LN_SS        SS       0    CD4 T       1603          836
## W380381   517  844      LN_SS        SS       0    CD4 T        840          513

# For older Seurat objects, you may need to run this
seuratObj <- UpdateSeuratObject(seuratObj)
```

Additionally, if your expression data has the older gene symbols, you
may want to use our alias conversion function to avoid the loss of gene
names.

``` r
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")
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

![](seurat_steps_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

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

![](seurat_steps_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

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

organism = "mouse"

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

weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
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

# Perform the NicheNet analysis

In contrary to NicheNet v1, we now recommend users to run both the
“sender-agnostic” approach and “sender-focused” approach. These
approaches only affect the list of potential ligands that are considered
for prioritization. As described in the flowchart above, we do not
define any sender populations in the ‘sender agnostic’ approach but
consider all ligands for which its cognature receptor is expressed. The
sender-focused approach will then filter the list of ligands to ones
where the ligands are expressed in the sender cell population(s).

## 1. Define a set of potential ligands for both the sender-agnostic and sender-focused approach

We first define a “receiver/target” cell population and determine which
genes are expressed. Here, we will consider a gene to be expressed if it
is expressed in at least 5% of cells (by default this is set to 10%).
The receiver cell population can only consist of one cell type, so in
case of multiple receiver populations, you will have to rerun the
vignette separately for each one. We will only look at CD8 T cells in
this vignette.

``` r
receiver = "CD8 T"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)
```

Get a list of all receptors available in the ligand-receptor network,
and define expressed receptors as genes that are in the ligand-receptor
network and expressed in the receiver. Then, define the potential
ligands as all ligands whose cognate receptors are expressed.

``` r
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
```

For the sender-focused approach, define sender cell types (CD4 T, Treg,
Mono, NK, B, and DC) and expressed genes in all sender populations.
(Although we pool all ligands from all sender cell types together in
this step, later on during the interpretation of the output, we will
check which sender cell type expresses which ligand.) Then, filter
potential ligands to those that are expressed in sender cells.

``` r
sender_celltypes <- c("CD4 T", "Treg", "Mono", "NK", "B", "DC")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
```

## 2. Define the gene set of interest

The gene set of interest are genes within the receiver cell type that
are likely to be influenced by ligands from the CCC event. In typical
case-control studies like this one, we use the differentially expressed
(DE) genes between the two conditions in the receiver cell type,
assuming that the observed DE pattern is a result of the CCC event
(i.e., LCMV infection). The condition of interest is thus ‘LCMV’,
whereas the reference/steady-state condition is ‘SS’. The condition can
be extracted from the metadata column ‘aggregate’. The method to
calculate the differential expression is here the standard Seurat
Wilcoxon test, but this can be changed if necessary.

``` r
condition_oi <-  "LCMV"
condition_reference <- "SS"

seurat_obj_receiver <- subset(seuratObj, idents = receiver)
seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["aggregate", drop=TRUE]])

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
```

## 3. Define the background genes

All expressed genes in the receiver cell population (that are also in
the ligand-target matrix) is defined as the ‘background set’ for our
ligand prioritization procedure in the next step.

``` r
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
```

## 4. Perform NicheNet ligand activity analysis

This is the main step of NicheNet where the potential ligands are ranked
based on the presence of their target genes in the gene set of interest
(compared to the background set of genes). In this case, we prioritize
ligands that induce the antiviral response in CD8 T cells.

Ligands are ranked based on the area under the precision-recall curve
(AUPR) between a ligand’s target predictions and the observed
transcriptional response. Although other metrics like the AUROC and
pearson correlation coefficient are also computed, we demonstrated in
our validation study that the AUPRwas the most informative measure to
define ligand activity (this was the Pearson correlation for v1). The
vignette on how we performed the validation can be found at [Evaluation
of NicheNet’s ligand-target predictions](model_evaluation.md).

We will first show the results of the sender-agnostic approach.

``` r
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
## # A tibble: 483 × 6
##    test_ligand auroc  aupr aupr_corrected pearson  rank
##    <chr>       <dbl> <dbl>          <dbl>   <dbl> <dbl>
##  1 Ifna1       0.714 0.433          0.358   0.498     1
##  2 Ifnb1       0.711 0.401          0.327   0.433     2
##  3 Ifnl3       0.683 0.392          0.317   0.433     3
##  4 Il27        0.682 0.391          0.316   0.445     4
##  5 Ifng        0.732 0.382          0.307   0.451     5
##  6 Ifnk        0.671 0.282          0.207   0.272     6
##  7 Ifne        0.667 0.279          0.204   0.289     7
##  8 Ebi3        0.666 0.264          0.189   0.256     8
##  9 Ifnl2       0.658 0.252          0.177   0.246     9
## 10 Ifna2       0.669 0.247          0.172   0.205    10
## # ℹ 473 more rows
```

We will use the top 30 ligands to predict active target genes and
construct an active ligand-receptor network.

``` r
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
```

We can also visualize the ligand activity measure (AUPR) of these
top-ranked ligands:

``` r
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## 5. Infer target genes and receptors of top-ranked ligands

### Active target gene inference

First, infer the active target genes, which are defined as genes in the
gene set of interest have the highest regulatory potential for each
top-ranked ligand. The function get_weighted_ligand_target_links will
return genes that are in the gene set of interest and are the top `n`
targets of a ligand (default: `n = 200`, but there are too many target
genes here so we only considered the top 100).

``` r
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
```

Visualize this in a heatmap. The rows of the heatmap are ordered based
on the rankings of the ligands, and the columns are ordered
alphabetically.

Note that not all ligands from the top 30 are present in the heatmap.
The left-out ligands are ligands that don’t have target genes with high
enough regulatory potential scores. Therefore, they did not survive the
used cutoffs. To include them, you can be less stringent in the used
cutoffs or increase the number of target genes considered.

``` r
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-16-1.png)<!-- --> We see
a lot of interferons in the top ligands, which biologically make sense
as we are looking at response to a viral infection.

### Receptors of top-ranked ligands

Similar to above, we identify which receptors have the highest
interaction potential with the top-ranked ligands.

``` r
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
```

Then, we create a heatmap for ligand-receptor interactions. Here, both
the ligands and receptors are ordered by hierarchical clustering You can
choose to order only ligands or receptors hierarachically (with
`order_hclust = ligands` or `receptors`, respectively) or not at all
(`none`), in which case the ligands are ordered based on their rankings,
and the receptors are ordered alphabetically..

``` r
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-18-1.png)<!-- --> \## 6.
Sender-focused approach

To perform the sender-focused approach, simply subset the ligand
activities to only contain expressed ligands from all populations
(calculated in Step 1. We can then perform target gene and receptor
inference as above.

``` r
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()
```

``` r
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank())

p_ligand_aupr
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Here, we instead observe that the top-ranked ligands consist of many H2
genes (which encode MHC-II proteins), and not IFN genes as in the
sender-agnostic approach. This is because IFN genes are not expressed by
the sender cell populations, and it was already filtered out during
preprocessing for being too lowly expressed.

``` r
best_upstream_ligands_all %in% rownames(seuratObj) %>% table()
## .
## FALSE  TRUE 
##    23     7
```

### Visualizing expression and log-fold change in sender cells

For the sender-focused approach, we can also investigate further on
which sender cell populations are potentially the true sender of these
ligands. First, we can simply check which sender cell population
expresses which of these top-ranked ligands.

``` r
# Dotplot of sender-agnostic approach
p_dotplot <- DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
        features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

As you can see, most op the top-ranked ligands seem to be mainly
expressed by dendritic cells and monocytes.

Next, we can also check upregulation of ligands in sender cells by
computing the log-fold change between the two conditions. This can add a
useful extra layer of information next to the ligand activities, because
you can assume that some of the ligands inducing DE in receiver cells,
will be DE themselves in the sender cells. This is of course only
possible in some cases, such as case-control studies.

``` r

celltype_order <- levels(Idents(seuratObj)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seuratObj,
  condition_colname = "aggregate",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ]) 

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                "Prioritized ligands", "LFC in Sender",
                                low_color = "midnightblue", mid_color = "white",
                                mid = median(vis_ligand_lfc), high_color = "red",
                                legend_title = "LFC")

p_lfc
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Finally, you can also compare rankings between the sender-agnostic and
sender-focused approach. Here, the red sections of the left bar plot
indicates which ligands in the sender-agnostic approach are filtered out
in the sender-focused approach because they are not expressed.

``` r
(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
   theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

## 7. Summary visualizations of the NicheNet analysis

Finally, we can make a combined plot containing heatmap of ligand
activities, ligand expression, ligand log-fold change and the target
genes of the top-ranked ligands. As mentioned earlier, sometimes ligands
do not appear in the ligand-target heatmap because they don’t have
target genes with high enough regulatory potential scores. In this case,
CCl22 is present in other plots (ranked 25th) but is missing in the
rightmost plot. If users wish for these plots to be consistent, they may
use the variable `order_ligands` defined when creating the ligand-target
heatmap to subset other plots instead of `best_upstream_ligands`.

``` r
figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  p_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot
```

![](seurat_steps_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

## Other follow-up analyses:

As another follow-up analysis, you can infer possible signaling paths
between ligands and targets of interest. You can read how to do this in
the following vignette [Inferring ligand-to-target signaling
paths](ligand_target_signaling_path.md):`vignette("ligand_target_signaling_path", package="nichenetr")`.

Another follow-up analysis is getting a “tangible” measure of how well
top-ranked ligands predict the gene set of interest and assess which
genes of the gene set can be predicted well. You can read how to do this
in the following vignette [Assess how well top-ranked ligands can
predict a gene set of
interest](target_prediction_evaluation_geneset.md):`vignette("target_prediction_evaluation_geneset", package="nichenetr")`.

In case you want to visualize ligand-target links between multiple
interacting cells, you can make an appealing circos plot as shown in
vignette [Circos plot visualization to show active ligand-target links
between interacting
cells](circos.md):`vignette("circos", package="nichenetr")`.

``` r
sessionInfo()
## R version 4.3.1 (2023-06-16)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS Stream 8
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Brussels
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] nichenetr_2.0.4    testthat_3.2.1     forcats_1.0.0      stringr_1.5.0      dplyr_1.1.4        purrr_1.0.2        readr_2.1.2        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.4.4      tidyverse_1.3.1    SeuratObject_5.0.1 Seurat_4.4.0      
## 
## loaded via a namespace (and not attached):
##   [1] IRanges_2.34.1              progress_1.2.3              ParamHelpers_1.14.1         nnet_7.3-19                 goftest_1.2-3               vctrs_0.6.5                 spatstat.random_3.2-2       digest_0.6.33              
##   [9] png_0.1-8                   shape_1.4.6                 proxy_0.4-27                OmnipathR_3.9.6             ggrepel_0.9.4               deldir_2.0-2                parallelly_1.36.0           magick_2.7.5               
##  [17] MASS_7.3-60                 reprex_2.0.1                reshape2_1.4.4              httpuv_1.6.13               foreach_1.5.2               BiocGenerics_0.46.0         withr_2.5.2                 xfun_0.41                  
##  [25] ggpubr_0.6.0                ellipsis_0.3.2              survival_3.5-5              memoise_2.0.1               zoo_1.8-12                  GlobalOptions_0.1.2         V8_4.3.3                    pbapply_1.7-2              
##  [33] Formula_1.2-5               prettyunits_1.2.0           promises_1.2.1              httr_1.4.7                  rstatix_0.7.2               globals_0.16.2              fitdistrplus_1.1-11         rstudioapi_0.15.0          
##  [41] miniUI_0.1.1.1              generics_0.1.3              base64enc_0.1-3             dir.expiry_1.8.0            curl_5.2.0                  S4Vectors_0.38.1            zlibbioc_1.46.0             ScaledMatrix_1.8.1         
##  [49] BBmisc_1.13                 polyclip_1.10-6             randomForest_4.7-1.1        GenomeInfoDbData_1.2.10     desc_1.4.3                  xtable_1.8-4                doParallel_1.0.17           evaluate_0.23              
##  [57] S4Arrays_1.2.0              hms_1.1.3                   GenomicRanges_1.52.0        irlba_2.3.5.1               colorspace_2.1-0            filelock_1.0.2              visNetwork_2.1.2            ROCR_1.0-11                
##  [65] smoof_1.6.0.3               reticulate_1.34.0           readxl_1.4.3                spatstat.data_3.0-3         magrittr_2.0.3              lmtest_0.9-40               mlrMBO_1.1.5.1              later_1.3.2                
##  [73] lattice_0.21-8              spatstat.geom_3.2-7         future.apply_1.11.0         lhs_1.1.6                   scattermore_1.2             XML_3.99-0.16               scuttle_1.10.2              shadowtext_0.1.2           
##  [81] cowplot_1.1.2               matrixStats_1.2.0           RcppAnnoy_0.0.21            class_7.3-22                Hmisc_5.1-0                 pillar_1.9.0                nlme_3.1-162                emoa_0.5-0.2               
##  [89] iterators_1.0.14            caTools_1.18.2              compiler_4.3.1              beachmat_2.16.0             stringi_1.7.6               gower_1.0.1                 tensor_1.5                  SummarizedExperiment_1.30.2
##  [97] lubridate_1.9.3             devtools_2.4.3              plyr_1.8.9                  crayon_1.5.2                abind_1.4-5                 parallelMap_1.5.1           haven_2.4.3                 locfit_1.5-9.8             
## [105] sp_2.1-2                    modelr_0.1.8                fastmatch_1.1-4             codetools_0.2-19            recipes_1.0.7               BiocSingular_1.16.0         e1071_1.7-14                GetoptLong_1.0.5           
## [113] plotly_4.10.0               mime_0.12                   splines_4.3.1               circlize_0.4.15             DiagrammeRsvg_0.1           Rcpp_1.0.11                 basilisk_1.12.1             dbplyr_2.1.1               
## [121] sparseMatrixStats_1.12.2    cellranger_1.1.0            knitr_1.45                  utf8_1.2.4                  clue_0.3-64                 fs_1.6.3                    listenv_0.9.0               checkmate_2.3.1            
## [129] DelayedMatrixStats_1.22.5   logger_0.2.2                pkgbuild_1.4.3              ggsignif_0.6.4              Matrix_1.6-4                statmod_1.5.0               tzdb_0.4.0                  DiceKriging_1.6.0          
## [137] tweenr_2.0.2                pkgconfig_2.0.3             tools_4.3.1                 cachem_1.0.8                viridisLite_0.4.2           rvest_1.0.2                 DBI_1.1.3                   fastmap_1.1.1              
## [145] rmarkdown_2.11              scales_1.3.0                grid_4.3.1                  usethis_2.2.2               ica_1.0-3                   liana_0.1.12                broom_0.7.12                patchwork_1.1.3            
## [153] dotCall64_1.1-1             carData_3.0-5               RANN_2.6.1                  rpart_4.1.19                farver_2.1.1                yaml_2.3.8                  MatrixGenerics_1.12.3       DiagrammeR_1.0.10          
## [161] foreign_0.8-84              cli_3.6.2                   stats4_4.3.1                leiden_0.3.9                lifecycle_1.0.4             caret_6.0-94                uwot_0.1.16                 Biobase_2.60.0             
## [169] mlr_2.19.1                  bluster_1.10.0              lava_1.7.3                  sessioninfo_1.2.2           backports_1.4.1             BiocParallel_1.34.2         timechange_0.2.0            gtable_0.3.4               
## [177] rjson_0.2.21                ggridges_0.5.5              progressr_0.14.0            parallel_4.3.1              pROC_1.18.5                 limma_3.56.2                jsonlite_1.8.8              edgeR_3.42.4               
## [185] bitops_1.0-7                assertthat_0.2.1            brio_1.1.4                  Rtsne_0.17                  spatstat.utils_3.0-4        BiocNeighbors_1.18.0        highr_0.10                  metapod_1.8.0              
## [193] dqrng_0.3.2                 timeDate_4032.109           lazyeval_0.2.2              shiny_1.7.1                 htmltools_0.5.7             sctransform_0.4.0           rappdirs_0.3.3              basilisk.utils_1.12.1      
## [201] glue_1.6.2                  spam_2.10-0                 XVector_0.40.0              RCurl_1.98-1.12             rprojroot_2.0.4             scran_1.28.2                gridExtra_2.3               igraph_1.2.11              
## [209] R6_2.5.1                    SingleCellExperiment_1.22.0 fdrtool_1.2.17              labeling_0.4.3              cluster_2.1.4               pkgload_1.3.3               GenomeInfoDb_1.36.1         ipred_0.9-14               
## [217] DelayedArray_0.26.7         tidyselect_1.2.0            htmlTable_2.4.1             ggforce_0.4.1               xml2_1.3.6                  car_3.1-2                   future_1.33.0               ModelMetrics_1.2.2.2       
## [225] rsvd_1.0.5                  munsell_0.5.0               KernSmooth_2.23-21          data.table_1.14.10          htmlwidgets_1.6.2           ComplexHeatmap_2.16.0       RColorBrewer_1.1-3          rlang_1.1.2                
## [233] spatstat.sparse_3.0-3       spatstat.explore_3.2-1      remotes_2.4.2               ggnewscale_0.4.9            fansi_1.0.6                 hardhat_1.3.0               prodlim_2023.08.28
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-medaglia_spatial_2017" class="csl-entry">

Medaglia, Chiara, Amir Giladi, Liat Stoler-Barak, Marco De Giovanni,
Tomer Meir Salame, Adi Biram, Eyal David, et al. 2017. “Spatial
Reconstruction of Immune Niches by Combining Photoactivatable Reporters
and <span class="nocase">scRNA</span>-Seq.” *Science*, December,
eaao4277. <https://doi.org/10.1126/science.aao4277>.

</div>

</div>
