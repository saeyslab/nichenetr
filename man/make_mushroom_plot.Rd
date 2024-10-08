% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/application_visualization.R
\name{make_mushroom_plot}
\alias{make_mushroom_plot}
\title{Make a "mushroom plot" of ligand-receptor interactions}
\usage{
make_mushroom_plot(
  prioritization_table,
  top_n = 30,
  show_rankings = FALSE,
  show_all_datapoints = FALSE,
  true_color_range = TRUE,
  use_absolute_rank = FALSE,
  size = "scaled_avg_exprs",
  color = "scaled_p_val_adapted",
  ligand_fill_colors = c("#DEEBF7", "#08306B"),
  receptor_fill_colors = c("#FEE0D2", "#A50F15"),
  unranked_ligand_fill_colors = c(alpha("#FFFFFF", alpha = 0.2), alpha("#252525", alpha =
    0.2)),
  unranked_receptor_fill_colors = c(alpha("#FFFFFF", alpha = 0.2), alpha("#252525", alpha
    = 0.2)),
  ...
)
}
\arguments{
\item{prioritization_table}{A prioritization table as generated by \code{\link{generate_prioritization_tables}}.}

\item{top_n}{An integer indicating how many ligand-receptor pairs to show}

\item{show_rankings}{A logical indicating whether to show the ranking of the ligand-receptor pairs (default: FALSE)}

\item{show_all_datapoints}{A logical indicating whether to show all ligand-receptor pairs (default: FALSE, if true they will be grayed out)}

\item{true_color_range}{A logical indicating whether to use the default color range as determined by ggplot (TRUE, default) or set the limits to a range of 0-1 (FALSE)}

\item{use_absolute_rank}{A logical indicating to whether use the absolute prioritization rank to filter the top_n ligand-receptor pairs (default: FALSE)}

\item{size}{A string indicating which column to use for the size of the semicircles (default: "scaled_avg_exprs"; use column name without "_ligand" or "_receptor" suffix)}

\item{color}{A string indicating which column to use for the color of the semicircles (default: "scaled_p_val_adapted"; use column name without "_ligand" or "_receptor" suffix)}

\item{ligand_fill_colors}{A vector of the low and high colors to use for the ligand semicircle fill gradient (default: c("#DEEBF7", "#08306B"))}

\item{receptor_fill_colors}{A vector of the low and high colors to use for the receptor semicircle fill gradient (default: c("#FEE0D2", "#A50F15"))}

\item{unranked_ligand_fill_colors}{A vector of the low and high colors to use for the unranked ligands when show_all_datapoints is TRUE (default: c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)))}

\item{unranked_receptor_fill_colors}{A vector of the low and high colors to use for the unranked receptors when show_all_datapoints is TRUE (default: c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)))}

\item{...}{Additional arguments passed to \code{\link{ggplot2::theme}}. As there are often issues with the scales legend, it is recommended to change legend sizes and positions using this argument, i.e., \code{legend.key.height}, \code{legend.key.width}, \code{legend.title}, and \code{legend.text}.}
}
\value{
A ggplot object
}
\description{
\code{make_mushroom_plot} Make a plot in which each glyph consists of two semicircles corresponding to ligand- and receptor- information. The size of the semicircle is the percentage of cells that express the protein, while the saturation corresponds to the scaled average expression value.
}
\details{
If the values range of the column used as the "size" parameter is not between 0 and 1.001, an error will be thrown.

The sender cell types can be ordered by encoding the "sender" column as a factor. If the "sender" column is not a factor, the sender cell types will be ordered alphabetically.

By default, the top_n ligand-receptor pairs are shown despite their absolute ranking. So, if a receiver cell type has LR pairs that are only ranked from 31-40 and the top_n is set to 20, the LR pairs will be shown. If use_absolute_rank is set to TRUE, only LR pairs with absolute ranking from 1-20 will be shown.
}
\examples{
\dontrun{
# Create a prioritization table
prior_table <- generate_prioritization_tables(processed_expr_table, processed_DE_table, ligand_activities, processed_condition_markers, prioritizing_weights)
make_mushroom_plot(prior_table)


# Show only top 20, and write rankings on the plot
make_mushroom_plot(prior_table, top_n = 20, show_rankings = TRUE)

# Show all datapoints, and use true color range
make_mushroom_plot(prior_table, show_all_datapoints = TRUE, true_color_range = TRUE)


# Change the size and color columns
make_mushroom_plot(prior_table, size = "pct_expressed", color = "scaled_avg_exprs")


# For a prioritization table with multiple receiver cell types
make_mushroom_plot(prior_table_combined \%>\% filter(receiver == celltype_oi))
}

}
