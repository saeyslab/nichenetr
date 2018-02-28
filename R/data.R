## description of data


#' Ligand-receptor network
#'
#' A data.frame/tibble describing ligand-receptor interactions and the data source describing this interaction.
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{from}{the ligand}
#'   \item{to}{receptor of the ligand}
#'   \item{source}{name of data source describing the interaction}
#'   }
#'
"lr_network"

#' Signaling network
#'
#' A data.frame/tibble describing gene-gene interactions such as common PPI, phosphorylation,...
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{from}{source gene}
#'   \item{to}{target gene}
#'   \item{source}{name of data source describing the interaction}
#'   }
#'
"sig_network"

#' Gene regulatory network
#'
#' A data.frame/tibble describing gene regulatory interactions such as TF-target interactions
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{from}{regulatory gene}
#'   \item{to}{target gene}
#'   \item{source}{name of data source describing the interaction}
#'   }
#'
"gr_network"
#' Data source weights
#'
#' A data.frame/tibble describing weight associated to each data source. These weights will be used for weighted aggregagion of source networks
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{source}{name of data source}
#'   \item{weight}{weight associated to a data source}
#'   }
#'
"source_weights_df"

#' Expression datasets for validation
#'
#' A list
#'
#' @format A list
#' \describe{
#'   \item{diffexp}{}
#'   \describe{
#'     \item{lfc}{Log fold change (treated vs untreated)}
#'     \item{qval}{Fdr corrected p-value}
#'     \item{gene}{Gene symbol}
#'   }
#'   \item{from}{Gene symbol(s) of ligand(s) by which the cells were treated}
#'   \item{name}{Name of the expression validation dataset}
#'   \item{type}{Type of dataset based on response time: "primary", "secondary", "primary + secondary"}
#'   }
#'
"expression_settings_validation"
#' Gene annotation information
#'
#' A data.frame/tibble describing weight associated to each data source. These weights will be used for weighted aggregagion of source networks
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{symbol}{human gene symbol}
#'   \item{entrez}{human gene entrez}
#'   \item{entrez_mouse}{mouse homolog gene entrez}
#'   \item{symbol_mouse}{mouse homolog gene symbol}
#'   }
#'
"geneinfo_human"
#' ExpressionSet object containing gene expression data from liver sinusoidal endothelial cells (LSECs) in wild type situation and also 12 hours and 36 hours after Kupffer cells were killed by giving Diphteria toxin to KC-DTR mice.
#'
#' Formal class 'ExpressionSet'
#'
"Exprs_lsec"
#' Limma voom object corresponding to the ExpressionSet object \code{Exprs_lsec}.
#'
#' Formal class 'EList'
#'
"voom_lsec"
#' ExpressionSet object containing gene expression data from bone marrow monocytes (data from 'van de Laar et al, 2016') and liver monocytes that differentiate into Kupffer cells after Kupffer cells were killed by giving Diphteria toxin to KC-DTR mice. Time points in differentiation process after treatment with toxin: 1 day, 3 days, 7 days, 15 days and 30 days.
#'
#' Formal class 'ExpressionSet'
#'
"Exprs_mono_kc"
