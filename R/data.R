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
#' Annotation table of all data sources used in the NicheNet model
#'
#' A data.frame/tibble describing for each type of data source the "parent" database and type of interactions it countains.
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{source}{Name of the data source}
#'   \item{database}{Name of the comprehensive database to which the data source belongs}
#'   \item{type_db}{Type of data source}
#'   \item{type_interaction}{Interaction type}
#'   \item{network}{Network layer to which the data source belongs: ligand_receptor, signaling or gene_regulatory}
#'   }
#'
"annotation_data_sources"
#' Optimized data source weights
#'
#' A data.frame/tibble describing weight associated to each data source. These weights will be used for weighted aggregagion of source networks and were determined via parameter optimization (see NicheNet paper).
#'
#' @format A data frame/tibble
#' \describe{
#'   \item{source}{name of data source}
#'   \item{weight}{weight associated to a data source}
#'   }
#'
"optimized_source_weights_df"
#' Optimized hyperparameter values
#'
#' A list
#'
#' @format A list
#' \describe{
#'   \item{lr_sig_hub}{ligand-signaling network hub correction factor}
#'   \item{gr_hub}{gene regulatory network hub correction factor}
#'   \item{ltf_cutoff}{cutoff on PPR scores}
#'   \item{damping_factor}{damping factor used in the PPR algorithm}
#'   }
#'
"hyperparameter_list"
