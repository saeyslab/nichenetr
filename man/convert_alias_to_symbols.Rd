% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supporting_functions.R
\name{convert_alias_to_symbols}
\alias{convert_alias_to_symbols}
\title{Convert aliases to official gene symbols}
\usage{
convert_alias_to_symbols(aliases, organism, verbose = TRUE)
}
\arguments{
\item{aliases}{A character vector of symbols/aliases}

\item{organism}{"mouse" or "human"}

\item{verbose}{TRUE or FALSE}
}
\value{
A character vector of official gene symbols.
}
\description{
\code{convert_alias_to_symbols} Convert aliases to offcial gene symbols
}
\examples{
library(dplyr)
human_symbols = c("TNF","IFNG","IL-15")
aliases = human_symbols \%>\% convert_alias_to_symbols(organism = "human")
}
