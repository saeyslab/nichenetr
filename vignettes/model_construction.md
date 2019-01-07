This vignette shows how ligand-target regulatory potential scores are inferred in the NicheNet framework. You can use the procedure shown here to develop your own model with inclusion of context-specific networks or removal of noisy irrelevant data sources. The networks at the basis of NicheNet can be downloaded from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484138.svg)](https://doi.org/10.5281/zenodo.1484138). Load the required packages and networks we will use to construct the model.

``` r
library(nichenetr)
library(dplyr)

# in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions
lr_network = readRDS(url("https://zenodo.org/record/1484138/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/1484138/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/1484138/files/gr_network.rds"))
```

Construct the weighted integrated ligand-signaling and gene regulatory network. In this vignette, we give every data source the same weight (as given by the `source_weights_df` data frame provided by default by the nichenetr package). The vignette showing how to use mlrMBO to optimize data source weights and other parameters will be written soon.

``` r
# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)

# downweigh the importance of signaling and gene regulatory hubs - the parameters used here are unoptimized
weighted_networks = apply_hub_corrections(weighted_networks, lr_sig_hub = 0.5, gr_hub = 0.5)
```

Infer ligand-target regulatory potential scores based on the weighted integrated networks

``` r
# in this example we will calculate target gene regulatory potential scores for TNF and the ligand combination TNF+IL6
ligands = list("TNF",c("TNF","IL6"))
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, algorithm = "PPR", damping_factor = 0.5)
```

Show some top target genes of the ligand TNF and the ligand combination TNF+IL6

``` r
extract_top_n_targets("TNF",10,ligand_target_matrix)
##     HACD4      P3H2       UBD      SELE     CCL19      IGHD    MUC5AC 
## 0.2236069 0.1666667 0.1209449 0.1156012 0.1136843 0.1091838 0.1070629 
##       CRP     CXCL9      COX1 
## 0.1049671 0.1044287 0.1021654
```

``` r
extract_top_n_targets("TNF-IL6",10,ligand_target_matrix)
##      BUD23      HACD4       IGHD        CRP       COX1       P3H2 
## 0.12500000 0.11180344 0.10925690 0.10597879 0.10255118 0.08342425 
##       SELE       IL11      CASP3       SELP 
## 0.07758420 0.07678502 0.07460599 0.07351416
```
