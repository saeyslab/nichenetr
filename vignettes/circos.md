Circos plot visualization to show active ligand-target links between
interacting cells
================
Robin Browaeys & Chananchida Sang-aram
2023-07-20

<!-- github markdown built using 
rmarkdown::render("vignettes/circos.Rmd", output_format = "github_document")
-->

This vignette shows how NicheNet can be used to predict active
ligand-target links between multiple interacting cells and how you can
make a circos plot to summarize the top-predicted links (via the
circlize package). This vignette starts in the same way as the main,
basis, NicheNet vignette [NicheNet’s ligand activity analysis on a gene
set of interest: predict active ligands and their target
genes](ligand_activity_geneset.md):`vignette("ligand_activity_geneset", package="nichenetr")`.
Make sure you understand the different steps described in that vignette
before proceeding with this vignette. In contrast to the basic vignette,
we will look communication between multiple cell types. More
specifically, we will predict which ligands expressed by both CAFs and
endothelial cells can induce the p-EMT program in neighboring malignant
cells (See Puram et al. 2017).

### Load packages required for this vignette

``` r
library(nichenetr)
library(tidyverse)
library(circlize)
```

### Read in expression data of interacting cells

First, we will read in the publicly available single-cell data from
CAFs, endothelial cells and malignant cells from HNSCC tumors.

``` r
hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells
```

Secondly, we will determine which genes are expressed in CAFs,
endothelial and malignant cells from high quality primary tumors.
Therefore, we wil not consider cells from tumor samples of less quality
or from lymph node metastases. To determine expressed genes, we use the
definition used by of Puram et al.

``` r
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% pull(cell)
endothelial_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "Endothelial") %>% pull(cell)
malignant_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_endothelial = expression[endothelial_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
```

### Load the ligand-target model we want to use

``` r
ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                     A2M        AANAT        ABCA1          ACE        ACE2
## A-GAMMA3'E 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
## A1BG       0.0018503922 0.0011108718 0.0014225077 0.0028594037 0.001139013
## A1BG-AS1   0.0007400797 0.0004677614 0.0005193137 0.0007836698 0.000375007
## A1CF       0.0024799266 0.0013026348 0.0020420890 0.0047921048 0.003273375
## A2M        0.0084693452 0.0040689323 0.0064256379 0.0105191365 0.005719199
```

### Load the gene set of interest and background of genes

As gene set of interest, we consider the genes of which the expression
is possibly affected due to communication with other cells.

Because we here want to investigate how CAFs and endothelial cells
regulate the expression of p-EMT genes in malignant cells, we will use
the p-EMT gene set defined by Puram et al. as gene set of interest and
use all genes expressed in malignant cells as background of genes.

``` r
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(pemt_geneset)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"

background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"
```

### Perform NicheNet’s ligand activity analysis on the gene set of interest

In a first step, we will define a set of potentially active ligands. As
potentially active ligands, we will use ligands that are 1) expressed by
CAFs and/or endothelial cells and 2) can bind a (putative) receptor
expressed by malignant cells. Putative ligand-receptor links were
gathered from NicheNet’s ligand-receptor data sources.

Note that we combine the ligands from CAFs and endothelial cells in one
ligand activity analysis now. Later on, we will look which of the
top-ranked ligands is mainly expressed by which of both cell types.

``` r
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands_CAFs = intersect(ligands,expressed_genes_CAFs)
expressed_ligands_endothelial = intersect(ligands,expressed_genes_endothelial)
expressed_ligands = union(expressed_ligands_CAFs, expressed_genes_endothelial)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "A2M"    "ACE"    "ADAM10" "ADAM12" "ADAM15" "ADAM17"
```

Now perform the ligand activity analysis: infer how well NicheNet’s
ligand-target potential scores can predict whether a gene belongs to the
p-EMT program or not.

``` r
ligand_activities = predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
```

Now, we want to rank the ligands based on their ligand activity. In our
validation study, we showed that the AUPR between a ligand’s target
predictions and the observed transcriptional response was the most
informative measure to define ligand activity. Therefore, we will rank
the ligands based on their AUPR. We will choose the top 20 ligands
here - as opposed to the top 30 in the main vignette - to avoid
overcrowding the circos plot.

``` r
ligand_activities %>% arrange(-aupr_corrected) 
## # A tibble: 232 × 5
##    test_ligand auroc   aupr aupr_corrected pearson
##    <chr>       <dbl>  <dbl>          <dbl>   <dbl>
##  1 TGFB2       0.768 0.123          0.107    0.199
##  2 CXCL12      0.708 0.0884         0.0721   0.144
##  3 BMP8A       0.770 0.0880         0.0718   0.177
##  4 INHBA       0.773 0.0866         0.0703   0.124
##  5 GDF3        0.758 0.0817         0.0654   0.156
##  6 LTBP1       0.722 0.0785         0.0622   0.163
##  7 ACE         0.711 0.0780         0.0617   0.151
##  8 TNXB        0.713 0.0737         0.0574   0.158
##  9 ENG         0.759 0.0732         0.0569   0.157
## 10 BMP5        0.745 0.0715         0.0552   0.150
## # … with 222 more rows
best_upstream_ligands = ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
head(best_upstream_ligands)
## [1] "TGFB2"  "CXCL12" "BMP8A"  "INHBA"  "GDF3"   "LTBP1"
```

We see here that the top-ranked ligands can predict the p-EMT genes
reasonably, this implies that ranking of the ligands might be accurate
as shown in our study. However, it is possible that for some gene sets,
the target gene prediction performance of the top-ranked ligands would
not be much better than random prediction. In that case, prioritization
of ligands will be less trustworthy.

Determine now which prioritized ligands are expressed by CAFs and or
endothelial cells

``` r
best_upstream_ligands %>% intersect(expressed_ligands_CAFs) 
##  [1] "TGFB2"  "CXCL12" "BMP8A"  "INHBA"  "LTBP1"  "TNXB"   "ENG"    "BMP5"   "VCAN"   "HGF"    "COMP"   "COL3A1" "MMP14"  "COL4A1" "FBN1"   "TIMP2"  "COL5A1"
best_upstream_ligands %>% intersect(expressed_ligands_endothelial)
##  [1] "CXCL12" "GDF3"   "LTBP1"  "ACE"    "TNXB"   "ENG"    "VCAN"   "HGF"    "COL4A1" "FBN1"   "TIMP2"  "EDN1"

# lot of overlap between both cell types in terms of expressed ligands
# therefore, determine which ligands are more strongly expressed in which of the two
ligand_expression_tbl = tibble(
  ligand = best_upstream_ligands, 
  CAF = expression[CAF_ids,best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
  endothelial = expression[endothelial_ids,best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}))

CAF_specific_ligands = ligand_expression_tbl %>% filter(CAF > endothelial + 2) %>% pull(ligand)
endothelial_specific_ligands = ligand_expression_tbl %>% filter(endothelial > CAF + 2) %>% pull(ligand)
general_ligands = setdiff(best_upstream_ligands,c(CAF_specific_ligands,endothelial_specific_ligands))

ligand_type_indication_df = tibble(
  ligand_type = c(rep("CAF-specific", times = CAF_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length()),
                  rep("Endothelial-specific", times = endothelial_specific_ligands %>% length())),
  ligand = c(CAF_specific_ligands, general_ligands, endothelial_specific_ligands))
```

### Infer target genes of top-ranked ligands and visualize in a circos plot

Now we will show how you can look at the regulatory potential scores
between ligands and target genes of interest. In this case, we will look
at links between top-ranked p-EMT-regulating ligands and p-EMT genes. In
this example, inferred target genes should belong to the p-EMT gene set
and to the 250 most strongly predicted targets of at least one of the
selected top-ranked ligands (the top 250 targets according to the
general prior model, so not the top 250 targets for this dataset).

Get first the active ligand-target links by looking which of the p-EMT
genes are among the top-predicted target genes for the prioritized
ligands:

``` r
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = pemt_geneset, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "p_emt") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
```

To avoid making a circos plots with too many ligand-target links, we
will show only links with a weight higher than a predefined cutoff:
links belonging to the 66% of lowest scores were removed. Not that this
cutoffs and other cutoffs used for this visualization can be changed
according to the user’s needs.

``` r
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
  
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
```

Prepare the circos visualization: give each segment of ligands and
targets a specific color and order

``` r
grid_col_ligand =c("General" = "lawngreen",
            "CAF-specific" = "royalblue",
            "Endothelial-specific" = "gold")
grid_col_target =c(
            "p_emt" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
```

Prepare the circos visualization: order ligands and targets

``` r
target_order = circos_links$target %>% unique()
ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)
```

Prepare the circos visualization: define the gaps between the different
segments

``` r
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CAF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Endothelial-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "p_emt") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
  )
```

Render the circos plot (all links same transparancy). Only the widths of
the blocks that indicate each target gene is proportional the
ligand-target regulatory potential (\~prior knowledge supporting the
regulatory interaction).

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 
```

![](circos_files/figure-gfm/unnamed-chunk-94-1.png)<!-- -->

``` r
circos.clear()
```

Render the circos plot (degree of transparancy determined by the
regulatory potential value of a ligand-target interaction)

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
```

![](circos_files/figure-gfm/unnamed-chunk-95-1.png)<!-- -->

``` r
circos.clear()
```

Save circos plot to an svg file

``` r
svg("ligand_target_circos.svg", width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()
## png 
##   2
```

### Visualize ligand-receptor interactions of the prioritized ligands in a circos plot

``` r
# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "p_emt_receptor") %>% inner_join(ligand_type_indication_df)
```

``` r
grid_col_ligand =c("General" = "lawngreen",
            "CAF-specific" = "royalblue",
            "Endothelial-specific" = "gold")
grid_col_receptor =c(
            "p_emt_receptor" = "darkred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
```

Prepare the circos visualization: order ligands and receptors

``` r
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)
```

Prepare the circos visualization: define the gaps between the different
segments

``` r
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  # width_ligand_receptor,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CAF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Endothelial-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_receptor,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "p_emt_receptor") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
  )
```

Render the circos plot (all links same transparancy). Only the widths of
the blocks that indicate each receptor is proportional the
ligand-receptor interaction weight (\~prior knowledge supporting the
interaction).

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1, order=order, link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
```

![](circos_files/figure-gfm/unnamed-chunk-101-1.png)<!-- -->

``` r
circos.clear()
```

Render the circos plot (degree of transparancy determined by the prior
interaction weight of the ligand-receptor interaction - just as the
widths of the blocks indicating each receptor)

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
```

![](circos_files/figure-gfm/unnamed-chunk-102-1.png)<!-- -->

``` r
circos.clear()
```

Save circos plot to an svg file

``` r
svg("ligand_receptor_circos.svg", width = 15, height = 15)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()
dev.off()
## png 
##   2
```

### Adding an outer track to the circos plot (ligand-receptor-target circos plot)

In the paper of Bonnardel, T’Jonck et al. [Stellate Cells, Hepatocytes,
and Endothelial Cells Imprint the Kupffer Cell Identity on Monocytes
Colonizing the Liver Macrophage
Niche](https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1), we
showed in Fig. 6B a ligand-receptor-target circos plot to visualize the
main NicheNet predictions. This “ligand-receptor-target” circos plot was
made by making first two separate circos plots: the ligand-target and
ligand-receptor circos plot. Then these circos plots were overlayed in
Inkscape (with the center of the two circles at the same location and
the ligand-receptor circos plot bigger than the ligand-target one). To
generate the combined circos plot as shown in Fig. 6B, we then manually
removed all elements of the ligand-receptor circos plot except the outer
receptor layer.

It is also possible to generate this kind of plot programmatically,
given that you are able to group your ligands. Below, we demonstrate two
approaches that uses either the `circlize::draw.sector` or
`circlize::highlight.sector` function. The former is more complicated
but gives you the flexibility to draw receptor arcs of different
lengths, while the `highlight.sector` function is constrained to the
widths of the targets in the inner track.

First, let’s rerun a code chunk from above to redefine `circos_links`
and the color scheme.

``` r
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
            "CAF-specific" = "royalblue",
            "Endothelial-specific" = "gold")
grid_col_target =c(
            "p_emt" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand, target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col = c(grid_ligand_color,grid_target_color)

ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
```

#### Using the `draw.sector` function

For demonstration purposes, we will use a subset of ligands and divide them into four groups that differ in the signaling pathways they interact with. Then, we assign targets and receptors to each ligand group based on their relative rankings (and summed weights in case the rankings are the same). The way we assigned targets and receptors to ligand groups here does not always lead to the most biologically meaningful results, as you will see in the final plot. Hence, it is best if you curate this list manually, e.g., through prior knowledge and by also looking at the ligand-target and ligand-receptor heatmap. Also keep in mind that there is no real correspondence between receptors and targets, as explained more [here](https://github.com/saeyslab/nichenetr/issues/20#issuecomment-611601039). 

``` r
groups <- list(group1 = c("TGFB2", "ENG"),
               group2 = grep("BMP|GDF|INHBA", best_upstream_ligands, value = TRUE),
               group3 = grep("COL|MMP|TIMP", best_upstream_ligands, value = TRUE),
               group4 = "CXCL12")

# Create list of targets for each group of ligand
targets <- lapply(names(groups), function(i) {
  # Rank each target for each ligand
  active_ligand_target_links_df %>% group_by(ligand) %>% mutate(target_rank = dense_rank(desc(weight))) %>%
    filter(ligand %in% groups[[i]]) %>%
    # Make two metrics for each target -> summed weight and avg_rank
    group_by(target) %>%
    summarise(n = n(), summed_weight=sum(weight), avg_rank = mean(target_rank)) %>%
    # Only keep targets that are connected to at least half of ligands in the group
    filter(n > (length(groups[[i]])/2)) %>% mutate(type = i)
}) %>% do.call(rbind, .)

# Check if any targets are distinct per group (groups 2 and 4 have some)
lapply(paste0("group", 1:5), function(group_name) setdiff(targets %>% filter(type == group_name) %>% pull(target), targets %>% filter(type != group_name) %>% pull(target)))
## [[1]]
## character(0)
## 
## [[2]]
## [1] "DKK3"   "GJA1"   "HTRA1"  "SEMA3C"
## 
## [[3]]
## character(0)
## 
## [[4]]
## [1] "LAMC2"   "MMP10"   "PRSS23"  "SLC31A2" "TPM4"   
## 
## [[5]]
## character(0)

# Assign target to a specific group, first based on ranking and second based on weight
targets_filtered <- targets %>% group_by(target) %>% filter(avg_rank == min(avg_rank)) %>%
  filter(summed_weight == max(summed_weight)) %>% ungroup()

# Do the same for receptors
receptor_colors <- c("#387D7A", "#9DA9A0", "#DD7373", "#725752") %>% setNames(names(groups))
receptors <- lapply(names(groups), function(i) {
  weighted_networks$lr_sig %>% filter(from %in% groups[[i]] & to %in% best_upstream_receptors) %>% 
    group_by(from) %>% mutate(receptor_rank = dense_rank(desc(weight))) %>%
    group_by(to) %>%  summarise(summed_weight=sum(weight),  avg_weight=mean(weight),  avg_rank = mean(receptor_rank)) %>%
    mutate(type=i, color = receptor_colors[i]) %>% rename(receptor = to)
}) %>% do.call(rbind, .)

# Check distinct receptors per group
lapply(paste0("group", 1:5), function(group_name) setdiff(receptors %>% filter(type == group_name) %>% pull(receptor), receptors %>% filter(type != group_name) %>% pull(receptor)))
## [[1]]
## [1] "MET"
## 
## [[2]]
## character(0)
## 
## [[3]]
## [1] "CD47"  "ITGA2" "ITGA3" "ITGB5" "ITGB6" "ITGB8"
## 
## [[4]]
## [1] "BDKRB2"
## 
## [[5]]
## character(0)
 
# Assign receptor to a specific group
receptors_filtered <-  receptors %>% group_by(receptor) %>% filter(avg_rank == min(avg_rank)) %>%
  filter(summed_weight == max(summed_weight)) %>% ungroup()
```

We will then have to redefine some variables.

``` r
# Filter out targets and ligands that are no longer present
links_circle_approach1 <- links_circle %>% filter(ligand %in% paste0(unlist(groups), "  "),
                                                  target %in% targets_filtered$target)

order <- c(ligand_order %>%.[. %in% paste0(unlist(groups), "  ")], targets_filtered$target)

# Redefine gaps between sectors
width_same_cell_same_ligand_type = 0.6
width_different_cell = 4.5
width_ligand_target = 12
width_same_cell_same_target_type = 0.6    # Added
width_different_target = 4.5              # Added

group_widths <- sapply(paste0("group", 1:4), function(group) {
  # Gap between targets of the same group
  paste0(rep(width_same_cell_same_target_type, times =(targets_filtered %>% filter(type==group) %>% nrow)-1), collapse=",")}) %>% 
  # Separate this with gap of different group
  paste0(., collapse=paste0(",",width_different_target,",")) %>%
  str_split(., ",") %>% .[[1]] %>% as.numeric

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CAF-specific",
                                                                         ligand %in% paste0(unlist(groups), "  "),
                                                                         target %in% targets_filtered$target) %>%
                                                   distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General",
                                                                         ligand %in% paste0(unlist(groups), "  "),
                                                                         target %in% targets_filtered$target) %>%
                                                   distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Endothelial-specific",
                                                                         ligand %in% paste0(unlist(groups), "  "),
                                                                         target %in% targets_filtered$target) %>%
                                                   distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  group_widths,
  width_ligand_target
  )

gaps = gaps %>% .[!is.na(.)]
```

Finally, we create the plot. What’s different here is we add an extra
layer in `preAllocateTracks`, and we add a `for` loop at the end to draw
the outer layer. As mentioned previously, the resulting plot will
require some further manual tweaking (e.g., ITG receptors from group 1
should be moved to group 3).

``` r
circos.par(gap.degree = gaps)

chordDiagram(links_circle_approach1, order=order, transparency=0, directional = 1, link.sort = TRUE, link.decreasing = FALSE,
             grid.col = grid_col, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands, annotationTrack = "grid",
    # Add extra track for outer layer
    preAllocateTracks = list(list(track.height = 0.025),
                             list(track.height = 0.25)))

# we go back to the first track and customize sector labels
circos.track(track.index = 2, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #


padding <- 0.05        # Gaps between receptor arcs
rou_adjustment <- 0.08 # Might need adjustment
for (group in unique(targets_filtered$type)){
  # Subset target and receptor
  targets_subset <- targets_filtered %>% filter(type == group)
  receptors_subset <- receptors_filtered %>% filter(type == group)
  
  # Get approximate position of outer ring
  pos <- circlize(c(0, 1), c(0, 1), track.index = 1)
  rou1 <- pos[1, "rou"]-rou_adjustment
  rou2 <- pos[2, "rou"]-rou_adjustment
  
  # Get range of angles from first to last target of the current group
  theta1 <- circlize:::get.sector.data((targets_subset %>% pull(target) %>% .[1]))["start.degree"]
  theta2 <- circlize:::get.sector.data((targets_subset %>% pull(target) %>% .[length(.)]))["end.degree"]
  
  # Scale the arc lengths according to the summed ligand-receptor weights
  receptors_subset_scaled <- receptors_subset %>% mutate(scaled_weight = summed_weight/sum(summed_weight)*(theta1-theta2))
  # For each receptor
  current_theta <- theta1
  for (i in 1:nrow(receptors_subset_scaled)){
    # Get end angle of the arc
    end_theta <- current_theta-(receptors_subset_scaled %>% slice(i) %>% pull(scaled_weight))
    d1 <- abs(end_theta-current_theta) # For gaps
    
    # Main function - we draw the arc here
    draw.sector(current_theta+(d1*-padding), end_theta-(d1*-padding),
                rou1, rou2,
                col = receptors_subset_scaled %>% slice(i) %>% pull(color),
                clock.wise = TRUE, border=NA)
    
    # Add text containing receptor name
    pos_text <- reverse.circlize((current_theta + end_theta)/2 + ifelse(current_theta < end_theta, 180, 0), (rou1 +  rou2)/2)
    # It's going to give a Note that point is out of plotting region
    suppressMessages(circos.text(pos_text[1,1], pos_text[1,2]+convert_y(7, "mm"),
                labels = receptors_subset_scaled %>% slice(i) %>% pull(receptor),
                niceFacing=TRUE, facing="clockwise", cex=0.6))
    
    current_theta <- end_theta
  }
  
}
```

![](circos_files/figure-gfm/unnamed-chunk-107-1.png)<!-- -->

``` r


circos.clear()
```

#### Using the `highlight.sector` function

With this function, it is not possible to draw receptor arcs that end at
different points from the target arcs. To demonstrate this, we will
randomly assign the target genes into one of three groups (Receptors A,
B, and C).

``` r
target_gene_groups <- sample(c("Receptor A", "Receptor B", "Receptor C"), length(unique(circos_links$target)), replace = TRUE) %>%
                          setNames(unique(circos_links$target))
target_gene_groups
##        ACTN1          C1S      COL17A1       COL1A1       COL4A2           F3        FSTL3       IGFBP3        ITGA5        LAMC2        MFAP2         MMP2         MYH9       PDLIM7 
## "Receptor B" "Receptor A" "Receptor A" "Receptor B" "Receptor C" "Receptor C" "Receptor B" "Receptor C" "Receptor C" "Receptor C" "Receptor B" "Receptor C" "Receptor C" "Receptor C" 
##        PSMD2        PTHLH     SERPINE1     SERPINE2        TAGLN        TGFBI          TNC         TPM1         MMP1        MMP10         MT2A       PRSS23      SLC31A2        THBS1 
## "Receptor A" "Receptor B" "Receptor C" "Receptor A" "Receptor C" "Receptor A" "Receptor C" "Receptor B" "Receptor C" "Receptor B" "Receptor C" "Receptor C" "Receptor C" "Receptor B" 
##         TPM4          APP       COL5A2         DKK3         GJA1        HTRA1         PLAU       SEMA3C          VIM         CAV1        ITGA6       MAGED1       MAGED2        PLOD2 
## "Receptor A" "Receptor B" "Receptor B" "Receptor B" "Receptor C" "Receptor B" "Receptor C" "Receptor C" "Receptor B" "Receptor A" "Receptor B" "Receptor A" "Receptor B" "Receptor A" 
##     SLC39A14        FSTL1       LGALS1        P4HA2         IL32         FHL2        ITGB1 
## "Receptor B" "Receptor C" "Receptor C" "Receptor C" "Receptor C" "Receptor B" "Receptor C"

target_gene_group_colors <- c("#387D7A", "#9DA9A0", "#704E2E") %>% setNames(unique(target_gene_groups))
```

Again, we will redefine some variables.

``` r
# Order targets according to receptor they belong to
order = c(ligand_order, target_gene_groups %>% sort %>% names)

# Redefine gaps between sectors
width_same_cell_same_ligand_type = 0.6
width_different_cell = 4.5
width_ligand_target = 12
width_same_cell_same_target_type = 0.6    # Added
width_different_target = 4.5              # Added

# Add this to circos_links
circos_links = circos_links %>% mutate(target_receptor = target_gene_groups[target])

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CAF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Endothelial-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  # Add code to define gaps between different target groups
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_receptor == "Receptor A") %>% distinct(target) %>% nrow() -1)),
  width_different_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_receptor == "Receptor B") %>% distinct(target) %>% nrow() -1)),
  width_different_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_receptor == "Receptor C") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
  )
```

The general idea here is similar - adding an extra layer in
`preAllocateTracks` and a `for` loop at the end to draw the outer
layer - but the function allows for much cleaner code.

``` r
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1, transparency=0, order=order,link.sort = TRUE, link.decreasing = FALSE,
             grid.col = grid_col, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
    # Add extra track for outer layer
    preAllocateTracks = list(list(track.height = 0.025),
                             list(track.height = 0.2)))

# we go back to the first track and customize sector labels
circos.track(track.index = 2, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #

# Add outer layer
for (target_gene_group in unique(target_gene_groups)){
  highlight.sector(target_gene_groups %>% .[. == target_gene_group] %>% names,
                   track.index = 1,
                 col = target_gene_group_colors[target_gene_group],
                 text = target_gene_group,
                 cex = 0.8, facing="bending.inside", niceFacing = TRUE, text.vjust = "5mm")
}
```

![](circos_files/figure-gfm/unnamed-chunk-110-1.png)<!-- -->

``` r

circos.clear()
```

### References

Bonnardel et al., 2019, Immunity 51, 1–17, [Stellate Cells, Hepatocytes,
and Endothelial Cells Imprint the Kupffer Cell Identity on Monocytes
Colonizing the Liver Macrophage
Niche](https://doi.org/10.1016/j.immuni.2019.08.017)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-puram_single-cell_2017" class="csl-entry">

Puram, Sidharth V., Itay Tirosh, Anuraag S. Parikh, Anoop P. Patel,
Keren Yizhak, Shawn Gillespie, Christopher Rodman, et al. 2017.
“Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor
Ecosystems in Head and Neck Cancer.” *Cell* 171 (7): 1611–1624.e24.
<https://doi.org/10.1016/j.cell.2017.10.044>.

</div>

</div>
