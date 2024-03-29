% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_nichenet.R
\name{get_non_spatial_de}
\alias{get_non_spatial_de}
\title{Makes a table similar to the output of `calculate_spatial_DE` and `process_spatial_de`, but now in case you don't have spatial information for the sender and/or receiver celltype. This is needed for comparability reasons.}
\usage{
get_non_spatial_de(niches, spatial_info, type, lr_network)
}
\arguments{
\item{niches}{a list of lists/niches giving the name, senders and receiver celltypes for each nice. Sender and receiver cell types should be part of Idents(seurat_obj).}

\item{spatial_info}{Tibble giving information about which celltypes should be compared to each other for defining spatial differential expression. Contains the columns "celltype_region_oi", "celltype_other_region", "niche", "celltype_type".}

\item{type}{For what type of cellype is the DE analysis: "sender" or "receiver"?}

\item{lr_network}{Ligand-Receptor Network in tibble format: ligand, receptor as columns}
}
\value{
A tibble of mock processed spatial DE information in case you don't have spatial information for the sender and/or receiver celltype.
}
\description{
\code{get_non_spatial_de} Makes a table similar to the output of `calculate_spatial_DE` and `process_spatial_de`, but now in case you don't have spatial information for the sender and/or receiver celltype. This is needed for comparability reasons.
}
\examples{
\dontrun{
niches = list(
"KC_niche" = list(
  "sender" = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
  "receiver" = c("KCs")),
"MoMac2_niche" = list(
  "sender" = c("Cholangiocytes","Fibroblast 2"),
  "receiver" = c("MoMac2")),
"MoMac1_niche" = list(
  "sender" = c("Capsule fibroblasts","Mesothelial cells"),
 "receiver" = c("MoMac1"))
seurat_obj = readRDS(url("https://zenodo.org/record/5840787/files/seurat_obj_subset_integrated_zonation.rds"))
spatial_info = tibble(celltype_region_oi = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
celltype_other_region = c("LSECs_central","Hepatocytes_central","Stellate cells_central")
) \%>\%
 mutate(niche =  "KC_niche", celltype_type = "sender")
get_non_spatial_de(niches, spatial_info, type = "receiver", lr_network)
}

}
