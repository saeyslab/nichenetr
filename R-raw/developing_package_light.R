

library(tidyverse)

## load in real data and save real and toy data in different locations

expression_settings_validation_real = readRDS("../paper/evaluation/ligand_treatment_datasets/expression_settings")
# devtools::use_data(expression_settings_validation_real,overwrite = T, compress = "bzip2", internal = TRUE)
expression_settings_validation = expression_settings_validation_real[c(1,2,3,4,54)]
devtools::use_data(expression_settings_validation,overwrite = T, compress = "bzip2")

lr_network_real = readRDS("../paper/networks/data/ligand_receptor/lr_network")
sig_network_real = readRDS("../paper/networks/data/signaling/signaling_network")
gr_network_real = readRDS("../paper/networks/data/gene_regulatory/gr_network")

source_weights_df = tibble(source = c(lr_network_real$source %>% unique(),sig_network_real$source %>% unique(), gr_network_real$source %>% unique()) %>% unique(), weight = 1)
devtools::use_data(source_weights_df,overwrite = T,compress = "bzip2")

geneinfo_human = readRDS("../paper/networks/data/annotation/complete_geneinfo")
devtools::use_data(geneinfo_human, overwrite = T, compress = "bzip2")

get_ncitations_genes = function(file = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz"){

  if (!is.character(file))
    stop("file must be a character vector")

  requireNamespace("dplyr")
  requireNamespace("readr")
  read_tsv(file, col_names = c("taxid", "entrez", "pubmedid"), skip = 1, col_types = cols(
    taxid = col_integer(),
    entrez = col_character(),
    pubmedid = col_integer())) %>% filter(taxid==9606) %>% group_by(entrez) %>% summarize(ncitations=n()) %>% left_join(geneinfo_human, "entrez") %>% arrange(-ncitations)
}
ncitations = get_ncitations_genes()
devtools::use_data(ncitations,overwrite = T, compress = "bzip2")


# devtools::use_data(lr_network_real,overwrite = T,compress = "bzip2", internal = TRUE)
lr_network = lr_network_real
devtools::use_data(lr_network,overwrite = T,compress = "bzip2")

# devtools::use_data(sig_network_real,overwrite = T,compress = "bzip2", internal = TRUE)
set.seed(1111)
sig_network = bind_rows(sig_network_real %>% filter(from %in% (expression_settings_validation %>% lapply(function(x){x$from}) %>% unlist() %>% unique())),sig_network_real %>% sample_n(10000))
devtools::use_data(sig_network,overwrite = T,compress = "bzip2")

# devtools::use_data(gr_network_real,overwrite = T,compress = "xz", internal = TRUE)
set.seed(1111)
gr_network = bind_rows(gr_network_real %>% filter(from %in% (expression_settings_validation %>% lapply(function(x){x$from}) %>% unlist() %>% unique())),gr_network_real %>% sample_n(10000))
devtools::use_data(gr_network,overwrite = T,compress = "xz")


# devtools::use_data(gr_network_real, sig_network_real,expression_settings_validation_real,overwrite = T, compress = "bzip2",internal = TRUE)


# load in the annotation table of the data sources
annotation_data_sources = read_tsv("../paper/networks/results/databases_sources_annotation.txt")
devtools::use_data(annotation_data_sources,overwrite = T,compress = "bzip2")


list_parameter_settings = readRDS("../paper/evaluation/evaluation/results/list_parameter_settings_nocv")
optimized_source_weights_df = tibble(source = list_parameter_settings$source_weights %>% names(),  weight = list_parameter_settings$source_weights)
hyperparameter_list = list(lr_sig_hub = list_parameter_settings$lr_sig_hub,
                           gr_hub = list_parameter_settings$gr_hub,
                           damping_factor = list_parameter_settings$damping_factor,
                           ltf_cutoff = list_parameter_settings$ltf_cutoff)
devtools::use_data(optimized_source_weights_df,overwrite = T,compress = "bzip2")
devtools::use_data(hyperparameter_list,overwrite = T,compress = "bzip2")

