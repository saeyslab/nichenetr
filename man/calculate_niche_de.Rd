% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_nichenet.R
\name{calculate_niche_de}
\alias{calculate_niche_de}
\title{Calculate differential expression of cell types in one niche versus all other niches of interest.}
\usage{
calculate_niche_de(seurat_obj, niches, type, assay_oi = "SCT")
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{niches}{a list of lists/niches giving the name, senders and receiver celltypes for each nice. Sender and receiver cell types should be part of Idents(seurat_obj).}

\item{type}{For what type of cellype is the DE analysis: "sender" or "receiver"?}

\item{assay_oi}{Which assay need to be used for DE calculation via `FindMarkers`. Default SCT, alternatives: RNA.}
}
\value{
A tibble containing the DE results of the niches versus each other.
}
\description{
\code{calculate_niche_de} Calculate differential expression of cell types in one niche versus all other niches of interest. This is possible for sender cell types and receiver cell types.
}
\examples{
\dontrun{
seurat_obj = readRDS(url("https://zenodo.org/record/5840787/files/seurat_obj_subset_integrated_zonation.rds"))
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
)
calculate_niche_de(seurat_obj, niches, "sender")
}

}
