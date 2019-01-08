Construction of NicheNet's ligand-target model
================
Robin Browaeys
2018-11-12

<!-- github markdown built using 
rmarkdown::render("vignettes/model_construction.Rmd", output_format = "github_document")
-->
This vignette shows how ligand-target regulatory potential scores are inferred in the NicheNet framework. You can use the procedure shown here to develop your own model with inclusion of context-specific networks or removal of noisy irrelevant data sources. The networks at the basis of NicheNet can be downloaded from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484138.svg)](https://doi.org/10.5281/zenodo.1484138).

Load the required packages and networks we will use to construct the model.

``` r
library(nichenetr)
library(dplyr)

# in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions

# Smaller toy networks are included in the package and will be used in this vignette:
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
head(sig_network)
## # A tibble: 6 x 4
##   from  to     source                  database      
##   <chr> <chr>  <chr>                   <chr>         
## 1 IL6   SH3GL2 omnipath_directional    omnipath      
## 2 IL6   SH3GL2 vinayagam_ppi           vinayagam     
## 3 TGFB1 PIK3R1 omnipath_directional    omnipath      
## 4 IL4   IL4R   omnipath_directional    omnipath      
## 5 IL4   IL4R   evex_regulation_binding evex_signaling
## 6 IL4   IL4R   evex_phosphorylation    evex_signaling
head(gr_network)
## # A tibble: 6 x 4
##   from  to     source database
##   <chr> <chr>  <chr>  <chr>   
## 1 BMP4  CHST11 trrust trrust  
## 2 IL4   ITGB3  trrust trrust  
## 3 TGFB1 COL1A2 trrust trrust  
## 4 TGFB1 FOXP3  trrust trrust  
## 5 TGFB1 KLF6   trrust trrust  
## 6 TGFB1 MYCN   trrust trrust

# The complete networks can be downloaded from Zenodo
# lr_network = readRDS(url("https://zenodo.org/record/1484138/files/lr_network.rds"))
# sig_network = readRDS(url("https://zenodo.org/record/1484138/files/signaling_network.rds"))
# gr_network = readRDS(url("https://zenodo.org/record/1484138/files/gr_network.rds"))
```

Construct the weighted integrated ligand-signaling and gene regulatory network. In this vignette, we give every data source the same weight (as given by the `source_weights_df` data frame provided by default by the nichenetr package). The vignette showing how to use mlrMBO to optimize data source weights and other parameters will be written in the near future.

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
##        LTBR       TRPC1       SNX20        CUBN       EXOC6      FUNDC2 
## 0.500008609 0.500008609 0.288680105 0.008342635 0.008342635 0.008342635 
##        GDF3      GEMIN5        HGH1      PRDM11       RIMS4        TEX2 
## 0.008342635 0.008342635 0.008342635 0.008342635 0.008342635 0.008342635 
##      TMEM70         TNR     ZSCAN22 
## 0.008342635 0.008342635 0.008342635
```

``` r
extract_top_n_targets("TNF-IL6",10,ligand_target_matrix)
##       CRP       IL4     ABCC1     DNMT1     SFTPB     STAT3     ICAM1 
## 0.5000091 0.5000091 0.3535599 0.3535599 0.3535599 0.3535599 0.3063494 
##    BCL2L1        TF       LEP 
## 0.2890819 0.2890172 0.2887815
```
