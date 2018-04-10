# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Reload all code:              'Ctrl + Shift + L'!!!! #   Reload all code:              'Ctrl + Shift + L'!!!!
#   Reload Tall code:              'Ctrl + Shift + D'!!!! #   Document all code  devtools::document()


library(tidyverse)
# install.packages("devtools")

# devtools::use_package("tidyverse")

devtools::use_package("dplyr")
devtools::use_package("tidyr")
devtools::use_package("tibble")
devtools::use_package("igraph")
devtools::use_package("Matrix")
devtools::use_package("fdrtool")
devtools::use_package("ROCR")
devtools::use_package("caTools")
devtools::use_package("data.table")
devtools::use_package("limma")
devtools::use_package("readr")
devtools::use_package("Hmisc")
devtools::use_package("caret")
devtools::use_package("purrr")
devtools::use_package("randomForest")
devtools::use_package("limma")
devtools::use_package("Biobase")
devtools::use_package("BiocGenerics")
devtools::use_package("DiagrammeR")
devtools::use_package("ggplot2")
devtools::use_package("mlrMBO")
devtools::use_package("parallelMap")

devtools::use_package("doMC","Suggests")
devtools::use_package("parallel", "Suggests")


# test data
# x = sample(1000)
# y = sample(500)
# devtools::use_data(x,mtcars,overwrite = T)
# devtools::use_data(y,internal = TRUE,overwrite = T) # you can add all datasets in same dataset here!

# real data
lr_network = readRDS("../staticNicheNet/results/networks/conf_ccc_complete_human_links_nov2017_121orths.rds")
sig_network = readRDS("../staticNicheNet/results/networks/conf_signaling_net_dec2017_121orths.rds")
gr_network = readRDS("../staticNicheNet/results/networks/conf_grn_human_nov2017_121orths.rds")

lr_network = lr_network %>% rename(source = evidence) %>% select(from,to,source) %>% mutate(from = humanentrez2humansymbol[from], to = humanentrez2humansymbol[to]) %>% drop_na()
sig_network = sig_network %>% rename(source = type) %>% select(from,to,source) %>% mutate(from = humanentrez2humansymbol[from], to = humanentrez2humansymbol[to]) %>% drop_na()
gr_network = gr_network  %>% rename(source = type) %>% select(from,to,source) %>% mutate(from = humanentrez2humansymbol[from], to = humanentrez2humansymbol[to]) %>% drop_na()

devtools::use_data(lr_network,sig_network,overwrite = T,compress = "bzip2")
devtools::use_data(gr_network,overwrite = T,compress = "xz")

source_weights_df = tibble(source = c(lr_network$source %>% unique,sig_network$source %>% unique(), gr_network$source %>% unique()) %>% unique(), weight = 1)
devtools::use_data(source_weights_df,overwrite = T,compress = "bzip2")

expression_settings_validation = readRDS("../staticNicheNet/results/expression_settings_validation")
devtools::use_data(expression_settings_validation,overwrite = T, compress = "bzip2")


Exprs_lsec = readRDS("../staticNicheNet/data/expression/niche/lsec_E.rds")
voom_lsec = readRDS("../staticNicheNet/data/expression/niche/lsec_v.rds")
devtools::use_data(Exprs_lsec,overwrite = T, compress = "bzip2")
devtools::use_data(voom_lsec,overwrite = T, compress = "bzip2")

Exprs_mono_kc = readRDS("../staticNicheNet/data/expression/immune/E_mono_kc.rds")
devtools::use_data(Exprs_mono_kc,overwrite = T, compress = "bzip2")


ncitations = get_ncitations_genes()
devtools::use_data(ncitations,overwrite = T, compress = "bzip2")


library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(dplyr)
library(purrr)
mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
geneinfo = tibble(symbol=unlist(as.list(org.Mm.egSYMBOL)), entrez=names(unlist(as.list(org.Mm.egSYMBOL))))
geneinfo_human = tibble(symbol=unlist(as.list(org.Hs.egSYMBOL)), entrez=names(unlist(as.list(org.Hs.egSYMBOL))))
# load("../staticNichenet/data/homologs.RData")
homologs = annotationTools::getHOMOLOG(geneinfo$entrez, 9606, read.delim("../staticNichenet/data/preproc/homologene/homologene.data", header=F))
mouseentrez2humanentrez = map_chr(homologs, ~.[[1]]) %>% setNames(geneinfo$entrez) %>% keep(!is.na(.))
humanentrez2mouseentrez = names(mouseentrez2humanentrez) %>% setNames(mouseentrez2humanentrez)
entrez2symbol = mapper(geneinfo, "symbol", "entrez")
symbol2entrez = mapper(geneinfo, "entrez", "symbol")
geneinfo_human = geneinfo_human %>% mutate(entrez_mouse=humanentrez2mouseentrez[entrez]) %>% mutate(symbol_mouse = entrez2symbol[entrez_mouse])
devtools::use_data(geneinfo_human, overwrite = T, compress = "bzip2")




