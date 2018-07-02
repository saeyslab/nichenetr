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
#'   \item{database}{name of the more comprehensive database the data source is part of}
#'   }
#'
"lr_network"

#' Signaling network
#'
#' A data.frame/tibble describing gene-gene interactions that are involved in signaling transduction (e.g. PPI, kinase-target,...)
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{from}{source gene}
#'   \item{to}{target gene}
#'   \item{source}{name of data source describing the interaction}
#'   \item{database}{name of the more comprehensive database the data source is part of}
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
#'   \item{database}{name of the more comprehensive database the data source is part of}
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
#' A data.frame/tibble describing HGNC human gene symbols, their entrez ids and the MGI mouse gene symbols and entrez ids of the one-one orthologs as determined via NCBI's homologene db and biomaRt ensembl db.
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
#' Number of citations for genes
#'
#' A data.frame/tibble describing the number of times a particular gene is cited in the pubmed literature (date: 13-03-2018)
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{entrez}{entrez}
#'   \item{ncitations}{ncitations}
#'   \item{symbol}{symbol}
#'   \item{entrez_mouse}{entrez_mouse}
#'   \item{symbol_mouse}{symbol_mouse}
#'   }
#'
"ncitations"
