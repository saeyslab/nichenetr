Inferring ligand-to-target signaling paths
================
Robin Browaeys
2019-01-17

<!-- github markdown built using 
rmarkdown::render("vignettes/ligand_target_signaling_path.Rmd", output_format = "github_document")
-->

### Infer signaling paths beween ligand(s) and target(s) of interest

To determine signaling paths between a ligand and target of interest, we
look at which transcription factors are best regulating the target genes
and are most closely downstream of the ligand (based on the weights of
the edges in the integrated ligand-signaling and gene regulatory
networks). Then, the shortest paths between these transcription factors
and the ligand of interests are determined and genes forming part in
this path are considered as important signaling mediators. Finally, we
look in our collected data source networks for all interactions between
the ligand, signaling mediators, transcription factors and target genes.
This allows to both prioritize signaling mediators and check which of
all collected data sources support the ligand-target predictions of
interest.

For this analysis, you need to define:

  - one or more ligands of interest
  - one or more target genes of interest

In this vignette, we will demonstrate how to infer signaling paths
between a CAF-ligand (CAF = cancer-associated fibroblast) of interest
and some of its top-predicted p-EMT target genes. The output of this
analysis can be easily imported into Cytoscape for exploration of the
networks.

First, we will load the necessary packages and networks to infer
signaling paths between ligand and target genes of interest.

``` r
library(nichenetr)
library(tidyverse)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
##   [3]
## names
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))
## dim
## dimnames

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
##   [1]
##   [2]
##   [3]
##   [4]
## row.names
## class
## names
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
##   [1]
##   [2]
##   [3]
##   [4]
## row.names
## class
## names
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
##   [1]
##   [2]
##   [3]
##   [4]
## row.names
## class
## names
```

As example, we will infer signaling paths between the CAF-ligand TGFB3
and its top-predicted p-EMT target genes TGFBI, LAMC2 and
TNC.

``` r
ligands_all = "TGFB3" # this can be a list of multiple ligands if required
targets_all = c("TGFBI","LAMC2","TNC")

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)
##   <environment: namespace:nichenetr>
## ligand_tf_matrix
## ligands_all
## targets_all
## top_n_regulators
## weighted_networks
## ligands_position
##   <environment: namespace:base>
## x
## do.NULL
## prefix
##   <environment: namespace:nichenetr>
## ligands_all
## targets_all
## k
## weighted_networks
## ligand_tf_matrix
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:base>
## ...
## KEEP.OUT.ATTRS
## stringsAsFactors
##   <environment: namespace:nichenetr>
## grid
## k
## weighted_networks
## ligand_tf_matrix
##   <environment: namespace:base>
## x
## rownames.force
## ...
##   <environment: namespace:base>
## x
## mode
##   <environment: namespace:base>
## x
## value
##   <environment: namespace:nichenetr>
## ligand_to_vis
## target_to_vis
## k
## weighted_networks
## ligand_tf_matrix
##   <environment: namespace:base>
## x
## ...
##   <environment: namespace:base>
## x
## value
##   <environment: namespace:dplyr>
## data
##   <environment: namespace:tibble>
## x
## ...
## validate
## .name_repair
## [1]
## [2]
## [3]
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
## validate
## ...
## .rows
## .name_repair
## rownames
##   <environment: namespace:tibble>
## .name_repair
## validate
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
## validate
## ...
## .rows
## .name_repair
##   <environment: namespace:tibble>
## x
## .rows
## .name_repair
## lengths
##   <environment: namespace:tibble>
## x
## .name_repair
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## name
## n
##   <environment: namespace:rlang>
## x
## y
##   <environment: namespace:rlang>
## x
## n
##   <environment: namespace:rlang>
## x
##   <environment: namespace:tibble>
## name
## .name_repair
##   <environment: namespace:tibble>
## name
##   <environment: namespace:tibble>
## name
##   <environment: namespace:tibble>
## name
## abort
##   <environment: namespace:tibble>
## name
## abort
##   <environment: namespace:rlang>
## x
## n
##   <environment: namespace:tibble>
## x
## [1]
## [2]
## [3]
##   <environment: namespace:tibble>
## .x
## .f
## ...
##   <environment: namespace:tibble>
## .x
## .f
## .mold
## ...
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## .x
## .f
## ...
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
## .rows
## lengths
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## .x
## .f
## ...
##   <environment: namespace:tibble>
## lengths
## .rows
##   <environment: namespace:tibble>
## x
## ...
## nrow
## class
## subclass
##   <environment: namespace:tibble>
## x
## ...
##   <environment: namespace:rlang>
## x
## n
## finite
##   <environment: namespace:tibble>
## x
## nrow
## subclass
##   <environment: namespace:pkgconfig>
## key
## fallback
##   <environment: namespace:pkgconfig>
## key
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
## names
##   <environment: namespace:dplyr>
## .data
## ...
##   <environment: namespace:tidyselect>
## .vars
## ...
## .strict
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:rlang>
## x
## ...
##   <environment: namespace:tidyselect>
## quo
## n
##   <environment: namespace:tidyselect>
## x
## n
##   <environment: namespace:purrr>
## .x
## .p
## .f
## ...
## .else
##   <environment: namespace:purrr>
## .x
## .at
## .f
## ...
##   <environment: namespace:purrr>
## nm
## .at
##   <environment: namespace:rlang>
## x
##   <environment: namespace:purrr>
## x
## sel
##   <environment: namespace:purrr>
## out
## x
##   <environment: namespace:purrr>
## .x
##   <environment: namespace:tidyselect>
## quos
## vars
##   <environment: namespace:tidyselect>
## vars
## frame
##   <environment: namespace:tidyselect>
## vars
##   <environment: namespace:tidyselect>
## vars
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
## names
##   <environment: namespace:tidyselect>
## quo
##   <environment: namespace:rlang>
## x
## default
##   <environment: namespace:purrr>
## .x
## .p
## .f
## ...
## .else
##   <environment: namespace:purrr>
## x
##   <environment: namespace:rlang>
## expr
## data
## env
##   <environment: namespace:tidyselect>
## expr
## name
## vars
##   <environment: namespace:rlang>
## .x
## ...
##   <environment: namespace:dplyr>
## df
## vars
##   <environment: namespace:dplyr>
## .data
## ...
##   <environment: namespace:dplyr>
## df
## dots
## caller_env
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
## names
##   <environment: namespace:dplyr>
## .data
## ...
##   <environment: namespace:tidyselect>
## .vars
## ...
## .include
## .exclude
## .strict
##   <environment: namespace:tidyselect>
## vars
## quos
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:dplyr>
## .x
## .p
## .f
## ...
##   <environment: namespace:dplyr>
## .x
## .p
## ...
##   <environment: namespace:dplyr>
## .x
## .f
## ...
##   <environment: namespace:dplyr>
## .x
## .f
## .mold
## ...
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:dplyr>
## .x
## .f
## ...
##   <environment: namespace:rlang>
## quo
## warn
##   <environment: namespace:dplyr>
## .x
## .f
## ...
##   <environment: namespace:tidyselect>
##   <environment: namespace:rlang>
## ...
##   <environment: namespace:rlang>
## ...
## .n_unnamed
## .ignore_empty
## .preserve_empty
## .homonyms
## .check_assign
##   <environment: namespace:tidyselect>
## x
## y
##   <environment: namespace:tidyselect>
## x
## y
##   <environment: namespace:tidyselect>
## ...
##   <environment: namespace:rlang>
## env
## data
## fn
## bind
## binder
## [1]
## [2]
## [3]
##   <environment: namespace:rlang>
## x
## fn
## env
##   <environment: namespace:tidyselect>
## quo
##   <environment: namespace:tidyselect>
## expr
##   <environment: namespace:rlang>
## x
## name
## n
## ns
##   <environment: namespace:tidyselect>
## nms
##   <environment: namespace:rlang>
## .parent
## ...
##   <environment: namespace:rlang>
## x
## parent
##   <environment: namespace:rlang>
## bottom
## top
## parent
##   <environment: namespace:rlang>
## data
##   <environment: namespace:tidyselect>
## x
## names
##   <environment: namespace:tidyselect>
## x
##   <environment: namespace:tidyselect>
## chr
## table
##   <environment: namespace:tidyselect>
## vars
## xs
##   <environment: namespace:tidyselect>
## x
## y
##   <environment: namespace:dplyr>
## .data
## ...
## .preserve
##   <environment: namespace:dplyr>
## ...
## .vectorised
##   <environment: namespace:dplyr>
## ...
## .op
## [1]
## [2]
## [3]
##   <environment: namespace:dplyr>
## df
## quo
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:dplyr>
## x
## y
## by
## copy
## suffix
## ...
## na_matches
##   <environment: namespace:dplyr>
## names
## warn_only
##   <environment: namespace:dplyr>
## by
## x
## y
##   <environment: namespace:dplyr>
## by
##   <environment: namespace:dplyr>
## by
## x
## y
##   <environment: namespace:dplyr>
## x
##   <environment: namespace:dplyr>
## na_matches
##   <environment: namespace:dplyr>
## x
## y
##   <environment: namespace:dplyr>
## x_names
## y_names
## by
## suffix
##   <environment: namespace:dplyr>
## x_names
## y_names
## by
##   <environment: namespace:dplyr>
## names
## by
##   <environment: namespace:dplyr>
## x
## y
## suffix
##   <environment: namespace:rlang>
## along
## x
##   <environment: namespace:dplyr>
## x
## y
## by_x
## by_y
## aux_x
## aux_y
## na_match
## frame
##   <environment: namespace:dplyr>
## x
## y
##   <environment: namespace:base>
## target
## current
## ...
##   <environment: namespace:base>
## x
##   <environment: namespace:tibble>
## x
## value
##   <environment: namespace:dplyr>
## out
## x
## vars
##   <environment: namespace:dplyr>
## .data
## ...
## .by_group
##   <environment: namespace:dplyr>
## df
## quosures
## frame
##   <environment: namespace:tibble>
## x
## i
## j
## drop
##   <environment: namespace:tibble>
## x
## i
##   <environment: namespace:tibble>
## x
## to
## nr
##   <environment: namespace:dplyr>
## dots
## id
##   <environment: namespace:igraph>
## d
## directed
## vertices
##   <environment: namespace:tibble>
## x
## row.names
## optional
## ...
##   <environment: namespace:base>
## x
##   <environment: namespace:igraph>
## n
## directed
##   <environment: namespace:igraph>
## percent
## message
## clean
##   <environment: namespace:igraph>
## graph
## nv
## ...
## attr
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
## edges
## ...
## attr
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
## v
## na.ok
##   <environment: namespace:nichenetr>
## ligand_oi
## signaling_df
## signaling_igraph
##   <environment: namespace:igraph>
## graph
## from
## to
## mode
## weights
## output
## predecessors
## inbound.edges
## [1]
## [2]
## [3]
##   <environment: namespace:igraph>
## arg
## choices
## several.ok
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
## P
## path
## directed
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
## es
## names
##   <environment: namespace:igraph>
## graph
## e
##   <environment: namespace:stats>
## object
## ...
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
## name
## index
##   <environment: namespace:igraph>
## vses
## graph
##   <environment: namespace:igraph>
## key
## value
## finalizer
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## x
## name
##   <environment: namespace:igraph>
## seq
##   <environment: namespace:igraph>
## ref
##   <environment: namespace:igraph>
## graph
## name
## index
##   <environment: namespace:igraph>
## seq
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
##   <environment: namespace:igraph>
## graph
## index
##   <environment: namespace:igraph>
## x
## name
##   <environment: namespace:igraph>
## seq
##   <environment: namespace:igraph>
## seq
##   <environment: namespace:igraph>
## x
## default
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
## [1]
## [2]
##   [6]
## [1]
## [2]
## [3]
##   [7]
##   [8]
##   [9]
##   [10]
##   [11]
##   [12]
##   [13]
##   [14]
##   [15]
## names
##   <environment: namespace:igraph>
## graph
## idx
## na_ok
##   <environment: namespace:igraph>
## x
## i
## na_ok
##   <environment: namespace:dplyr>
## .data
## ...
## add
## .drop
##   <environment: namespace:dplyr>
## dots
## env
## ...
## .named
## [1]
## [2]
## [3]
##   <environment: namespace:ggplot2>
## x
##   <environment: namespace:rlang>
## x
## n
##   <environment: namespace:rlang>
## quo
##   <environment: namespace:dplyr>
## .data
## vars
##   <environment: namespace:dplyr>
## quo
##   <environment: namespace:rlang>
## exprs
## width
## printer
##   <environment: namespace:rlang>
## .x
## .f
## ...
##   <environment: namespace:rlang>
## .x
## .f
## .mold
## ...
##   <environment: namespace:rlang>
## x
##   <environment: namespace:rlang>
## x
## parent
## warn
##   <environment: namespace:rlang>
## .x
## ...
##   <environment: namespace:rlang>
## x
##   <environment: namespace:rlang>
## x
## shallow
##   <environment: namespace:rlang>
## x
##   <environment: namespace:dplyr>
## data
## symbols
## drop
##   <environment: namespace:dplyr>
## .tbl
##   <environment: namespace:dplyr>
## x
## ...
##   <environment: namespace:dplyr>
## df
##   <environment: namespace:dplyr>
## .data
## ...
## .keep_all
##   <environment: namespace:dplyr>
## .data
## vars
## group_vars
## .keep_all
##   <environment: namespace:dplyr>
## df
## keep_cols
##   <environment: namespace:tibble>
## j
## x
##   <environment: namespace:tibble>
## j
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:dplyr>
## vars
## data
##   <environment: namespace:dplyr>
## df
## vars
## keep
## frame
##   <environment: namespace:dplyr>
## x
## ...

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
##   <environment: namespace:nichenetr>
## signaling_graph_list
## ligands_all
## targets_all
## sig_color
## gr_color
##   <environment: namespace:DiagrammeR>
## n
## type
## label
## ...
##   <environment: namespace:DiagrammeR>
##   <environment: namespace:stringr>
## string
## pattern
## replacement
##   <environment: namespace:stringi>
## str
## pattern
## replacement
## vectorize_all
## ...
## opts_regex
##   <environment: namespace:stringr>
## x
##   <environment: namespace:stringr>
## x
##   <environment: namespace:stringr>
## string
## pattern
## n
## simplify
##   <environment: namespace:stringi>
## str
## n
## tokens_only
## simplify
## ...
## opts_brkiter
##   <environment: namespace:stringi>
## type
## locale
## skip_word_none
## skip_word_number
## skip_word_letter
## skip_word_kana
## skip_word_ideo
## skip_line_soft
## skip_line_hard
## skip_sentence_term
## skip_sentence_sep
## ...
##   <environment: namespace:dplyr>
## ...
##   <environment: namespace:dplyr>
## dots
##   <environment: namespace:tibble>
## x
## prefix
## sep
##   <environment: namespace:tibble>
## name
## prefix
## sep
##   <environment: namespace:DiagrammeR>
## from
## to
## rel
## ...
##   <environment: namespace:DiagrammeR>
## nodes_df
## edges_df
## directed
## graph_name
## attr_theme
## write_backups
##   <environment: namespace:base>
## location
## [1]
## [2]
## [3]
##   <environment: namespace:base>
## n
## expr
## simplify
##   <environment: namespace:base>
## x
## row.names
## optional
## ...
## nm
##   <environment: namespace:base>
## x
## row.names
## optional
## ...
## nm
##   <environment: namespace:tibble>
## xs
## transform
##   <environment: namespace:rlang>
## names
## x
##   <environment: namespace:tibble>
## x
## i
##   <environment: namespace:base>
## x
## ...
## drop
##   <environment: namespace:DiagrammeR>
## graph_log
## version_id
## function_used
## time_modified
## duration
## nodes
## edges
## d_n
## d_e
##   <environment: namespace:DiagrammeR>
## start_time
##   <environment: namespace:base>
## e1
## e2
## [1]
## [2]
## [3]

# To render the graph:
# DiagrammeR::render_graph(graph_min_max, layout = "tree")
```

We will now look which of the collected data sources support the
interactions in this
network.

``` r
data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
##   <environment: namespace:nichenetr>
## signaling_graph_list
## lr_network
## sig_network
## gr_network
head(data_source_network) 
##   <environment: namespace:knitr>
## x
## ...
##   <environment: namespace:data.table>
## x
##   [1]
## [1]
## [2]
##   [2]
##   [3]
##   [4]
##   [5]
## names
##   <environment: namespace:knitr>
## x
## format
## digits
## row.names
## col.names
## align
## caption
## label
## format.args
## escape
## ...
##   <environment: namespace:tibble>
## x
## ...
## n
## width
## n_extra
##   <environment: namespace:tibble>
## ...
##   <environment: namespace:tibble>
## x
## ...
## n
## width
## n_extra
##   <environment: namespace:tibble>
## x
##   [1]
##   [2]
##   [3]
##   [4]
## names
##   <environment: namespace:tibble>
## df
## rows
## n
## star
##   <environment: namespace:tibble>
## .data
##   <environment: namespace:pillar>
## x
## has_row_id
## width
## ...
## [1]
## [2]
## [3]
## [4]
## [5]
##   <environment: 0x000001dd4acc6670>
## forget
##   <environment: namespace:base>
## con
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## .x
## .y
## .f
## ...
##   <environment: namespace:pillar>
## x
## name
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:tibble>
## .data
##   <environment: namespace:pillar>
## x
## width
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## .x
## .f
## ...
##   <environment: namespace:tibble>
## x
## ...
##   <environment: namespace:base>
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
##   <environment: namespace:tibble>
## x
## width
## ...
##   <environment: namespace:tibble>
## width
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
## right
## space
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
## width
##   <environment: namespace:tibble>
## ...
## indent
## prefix
## width
##   <environment: namespace:tibble>
## x
## width
## indent
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
##   [7]
##   [8]
## names
## class
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
##   [6]
##   [7]
##   [8]
##   [9]
##   [10]
## [1]
##   [11]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
##   [7]
##   [8]
## names
## class
##   [1]
## names
##   [1]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
## [11]
## [12]
## [13]
## [14]
## [15]
## [16]
## [17]
## [18]
## [19]
## [20]
## [21]
## [22]
## [23]
## [24]
## [25]
## [26]
## [27]
## [28]
## [29]
## [30]
## [31]
## [32]
## [33]
## [34]
## [35]
## [36]
## [37]
## [38]
## [39]
## [40]
## [41]
## [42]
## [43]
## [44]
## [45]
## [46]
## [47]
## [48]
## [49]
## [50]
## [51]
## [52]
## [53]
## [54]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
##   [7]
##   [8]
## names
## class
##   <environment: namespace:fansi>
## libname
## pkgname
##   <environment: namespace:fansi>
##   <environment: namespace:fansi>
## x
## width
## indent
## exdent
## prefix
## simplify
## initial
## warn
## term.cap
## ctl
##   <environment: namespace:pillar>
## x
## width
## ...
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:pillar>
## x
##   <environment: namespace:rlang>
## x
##   <environment: namespace:rlang>
## x
##   <environment: namespace:pillar>
## n
## has_title_row
## has_star
## ...
##   <environment: namespace:pillar>
## has_title_row
## has_star
## ...
##   <environment: namespace:pillar>
## has_title
## ...
##   <environment: namespace:pillar>
## has_star
## ...
##   <environment: namespace:pillar>
## title
## type
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## min_width
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## n
## ...
##   <environment: namespace:pillar>
## x
## ...
## width
## min_width
## class
## subclass
##   <environment: namespace:rlang>
## x
## n
##   <environment: namespace:pillar>
## capital
## shaft
## width
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## .x
## .f
## ...
##   <environment: namespace:pillar>
## .x
## .f
## .mold
## ...
##   <environment: namespace:pillar>
## x
## width
## rowid_width
##   <environment: namespace:pillar>
## width
## ncol
## rowid_width
## tier_width
##   <environment: namespace:pillar>
## .x
## .f
##   <environment: namespace:pillar>
## x
## title
## ...
##   <environment: namespace:pillar>
## x
## ...
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## width
##   <environment: namespace:pillar>
## x
## width
## align
##   <environment: namespace:pillar>
## x
## width
##   <environment: namespace:fansi>
## x
## warn
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
##   [7]
##   [8]
## names
## class
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
##   [6]
##   [7]
##   [8]
##   [9]
##   [10]
## [1]
##   [11]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
##   [7]
##   [8]
## names
## class
##   [1]
## names
##   [1]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
## [11]
## [12]
## [13]
## [14]
## [15]
## [16]
## [17]
## [18]
## [19]
## [20]
## [21]
## [22]
## [23]
## [24]
## [25]
## [26]
## [27]
## [28]
## [29]
## [30]
## [31]
## [32]
## [33]
## [34]
## [35]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
##   [7]
##   [8]
## names
## class
##   <environment: namespace:utf8>
## x
## encode
## quote
## utf8
##   <environment: namespace:utf8>
## expr
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
## [1]
## [2]
## [3]
##   <environment: namespace:utf8>
## name
## value
##   <environment: namespace:utf8>
## name
## value
##   <environment: namespace:pillar>
## width
##   <environment: namespace:pillar>
## x
## ...
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [8]
## names
## class
##   [1]
## [1]
## [2]
## [3]
## [4]
##   [2]
##   [3]
##   [4]
##   [5]
##   [6]
##   [7]
##   [8]
##   [9]
##   [10]
## [1]
##   [11]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [8]
## names
## class
##   [1]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [7]
##   [8]
## names
## class
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
##   [6]
##   [7]
##   [8]
##   [9]
##   [10]
##   [11]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
##   [4]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [7]
##   [8]
## names
## class
##   [1]
## names
##   [1]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
## [11]
## [12]
## [13]
## [14]
## [15]
## [16]
## [17]
## [18]
## [19]
## [20]
## [21]
## [22]
## [23]
## [24]
## [25]
## [26]
## [27]
## [28]
## [29]
## [30]
## [31]
## [32]
## [33]
## [34]
## [35]
## [36]
## [37]
## [38]
## [39]
## [40]
## [41]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## names
##   <environment: namespace:zeallot>
## x
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
## names
##   [1]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
## [9]
## [10]
## [11]
## [12]
## [13]
## [14]
## [15]
## [16]
## [17]
## [18]
## [19]
## [20]
## [21]
## [22]
## [23]
## [24]
## [25]
## [26]
## [27]
## [28]
## [29]
## [30]
## [31]
## [32]
## [33]
## [34]
## [35]
## [36]
## [37]
## [38]
## [39]
## [40]
## [41]
## [42]
## [43]
## [44]
## [45]
## [46]
## [47]
## [48]
## [49]
## [50]
## [51]
## [52]
## [53]
## [54]
## [55]
## [56]
## [57]
## [58]
## [59]
## [60]
## [61]
## [62]
## [63]
## [64]
## [65]
## [66]
## [67]
## [68]
## [69]
## [70]
## [71]
## [72]
## [73]
## [74]
## [75]
## [76]
## [77]
## [78]
## [79]
## [80]
## [81]
## [82]
## [83]
## [84]
## [85]
## [86]
## [87]
## [88]
## [89]
## [90]
## [91]
## [92]
## [93]
## [94]
## [95]
## [96]
## [97]
## [98]
## [99]
## [100]
## [101]
## [102]
## [103]
## [104]
## [105]
## [106]
## [107]
## [108]
## [109]
## [110]
## [111]
## [112]
## [113]
## [114]
## [115]
## [116]
## [117]
## [118]
## [119]
## [120]
## [121]
## [122]
## [123]
## [124]
## [125]
## [126]
## [127]
## [128]
## [129]
## [130]
## [131]
## [132]
## [133]
## [134]
## [135]
## [136]
## [137]
## [138]
## [139]
## [140]
## [141]
## [142]
## [143]
## [144]
## [145]
## [146]
## [147]
## [148]
## [149]
## [150]
## [151]
## [152]
## [153]
## [154]
## [155]
## [156]
## [157]
## [158]
## [159]
## [160]
## [161]
## [162]
## [163]
## [164]
## [165]
## [166]
## [167]
## [168]
## [169]
## [170]
## [171]
## [172]
## [173]
## [174]
## [175]
## [176]
## [177]
## [178]
## [179]
## [180]
## [181]
## [182]
## [183]
## [184]
## [185]
## [186]
## [187]
## [188]
## [189]
## [190]
## [191]
## [192]
## [193]
## [194]
## [195]
## [196]
## [197]
## [198]
## [199]
## [200]
## [201]
## [202]
## [203]
## [204]
## [205]
## [206]
## [207]
## [208]
## [209]
## [210]
## [211]
## [212]
## [213]
## [214]
## [215]
## [216]
## [217]
## [218]
## [219]
## [220]
## [221]
## [222]
## [223]
## [224]
## [225]
## [226]
## [227]
## [228]
## [229]
## [230]
## [231]
## [232]
## [233]
## [234]
## [235]
## [236]
## [237]
## [238]
## [239]
## [240]
## [241]
## [242]
## [243]
## [244]
## [245]
## [246]
## [247]
## [248]
## [249]
## [250]
## [251]
## [252]
## [253]
## [254]
## [255]
## [256]
## [257]
## [258]
## [259]
## [260]
## [261]
## [262]
## [263]
## [264]
## [265]
## [266]
## [267]
## [268]
## [269]
## [270]
## [271]
## [272]
## [273]
## [274]
## [275]
## [276]
## [277]
## [278]
## [279]
## [280]
## [281]
## [282]
## [283]
## [284]
## [285]
## [286]
## [287]
## [288]
## [289]
## [290]
## [291]
## [292]
## [293]
## [294]
## [295]
## [296]
## [297]
## [298]
## [299]
## [300]
## [301]
## [302]
## [303]
## [304]
## [305]
## [306]
## [307]
## [308]
## [309]
## [310]
## [311]
## [312]
## [313]
## [314]
## [315]
## [316]
## [317]
## [318]
## [319]
## [320]
## [321]
## [322]
## [323]
## [324]
## [325]
## [326]
## [327]
## [328]
## [329]
## [330]
## [331]
## [332]
## [333]
## [334]
## [335]
## [336]
## [337]
## [338]
## [339]
## [340]
## [341]
## [342]
## [343]
## [344]
## [345]
## [346]
## [347]
## [348]
## [349]
## [350]
## [351]
## [352]
## [353]
## [354]
## [355]
## [356]
## [357]
## [358]
## [359]
## [360]
## [361]
## [362]
## [363]
## [364]
## [365]
## [366]
## [367]
## [368]
## [369]
## [370]
## [371]
## [372]
## [373]
## [374]
## [375]
## [376]
## [377]
## [378]
## [379]
## [380]
## [381]
## [382]
## [383]
## [384]
## [385]
## [386]
## [387]
## [388]
## [389]
## [390]
## [391]
## [392]
## [393]
## [394]
## [395]
## [396]
## [397]
## [398]
## [399]
## [400]
## [401]
## [402]
## [403]
## [404]
## [405]
## [406]
## [407]
## [408]
## [409]
## [410]
## [411]
## [412]
## [413]
## [414]
## [415]
## [416]
## [417]
## [418]
## [419]
## [420]
## [421]
## [422]
## [423]
## [424]
## [425]
## [426]
## [427]
## [428]
## [429]
## [430]
## [431]
## [432]
## [433]
## [434]
## [435]
## [436]
## [437]
## [438]
## [439]
## [440]
## [441]
## [442]
## [443]
## [444]
## [445]
## [446]
## [447]
## [448]
## [449]
## [450]
## [451]
## [452]
## [453]
## [454]
## [455]
## [456]
## [457]
## [458]
## [459]
## [460]
## [461]
## [462]
## [463]
## [464]
## [465]
## [466]
## [467]
## [468]
## [469]
## [470]
## [471]
## [472]
## [473]
## [474]
## [475]
## [476]
## [477]
## [478]
## [479]
## [480]
## [481]
## [482]
## [483]
## [484]
## [485]
## [486]
## [487]
## [488]
## [489]
## [490]
## [491]
## [492]
## [493]
## [494]
## [495]
## [496]
## [497]
## [498]
## [499]
## [500]
## [501]
## [502]
## [503]
## [504]
## [505]
## [506]
## [507]
## [508]
## [509]
## [510]
## [511]
## [512]
## [513]
## [514]
## [515]
## [516]
## [517]
## [518]
## [519]
## [520]
## [521]
## [522]
## [523]
## [524]
## [525]
## [526]
## [527]
## [528]
## [529]
## [530]
## [531]
## [532]
## [533]
## [534]
## [535]
## [536]
## [537]
## [538]
## [539]
## [540]
## [541]
## [542]
## [543]
## [544]
## [545]
## [546]
## [547]
## [548]
## [549]
## [550]
## [551]
## [552]
## [553]
## [554]
## [555]
## [556]
## [557]
## [558]
## [559]
## [560]
## [561]
## [562]
## [563]
## [564]
## [565]
## [566]
## [567]
## [568]
## [569]
## [570]
## [571]
## [572]
## [573]
## [574]
## [575]
## [576]
## [577]
## [578]
## [579]
## [580]
## [581]
## [582]
## [583]
## [584]
## [585]
## [586]
## [587]
## [588]
## [589]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## names
##   <environment: namespace:vctrs>
## x
## outer_name
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## op
## x
## y
##   <environment: namespace:vctrs>
## op
## x
## y
##   <environment: namespace:vctrs>
## op
## x
## y
##   <environment: namespace:vctrs>
## op
## x
## y
##   <environment: namespace:vctrs>
## op
## x
## y
##   <environment: namespace:vctrs>
## op
## x
## y
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## fun
## x
## ...
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
## to
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
## y
##   <environment: namespace:vctrs>
## x
##   [1]
##   [2]
##   [3]
##   [4]
##   [5]
## names
##   [1]
##   [2]
## [1]
## [2]
## [3]
## [4]
##   [3]
## [1]
## [2]
## [3]
##   [4]
## [1]
##   [5]
##   [6]
## [1]
## [2]
## [3]
## [4]
## [5]
## [6]
## [7]
## [8]
##   [7]
## [1]
## [2]
## [3]
## [4]
## [5]
##   [8]
## names
## class
##   <environment: namespace:vctrs>
## libname
## pkgname
##   <environment: namespace:vctrs>
## generic
## class
## method
## [1]
## [2]
## [3]
##   <environment: namespace:vctrs>
## x
## ...
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
## levels
## ...
##   <environment: namespace:vctrs>
## x
## levels
## ...
##   <environment: namespace:vctrs>
## x
## units
## ...
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:vctrs>
## x
##   <environment: namespace:rlang>
## .x
## .to
## ...
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
## [1]
## [2]
## [3]
## [4]
## [5]
##   <environment: 0x000001dd48ed60f8>
## x
##   <environment: namespace:pillar>
## pillars
## tier_widths
##   <environment: namespace:pillar>
## .x
## .f
## ...
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## col_df
## tier_widths
## fixed_tier_df
##   <environment: namespace:pillar>
## widths
## tier_widths
##   <environment: namespace:pillar>
## tier_df
##   <environment: namespace:pillar>
## x
## ...
## min_width
##   <environment: namespace:utf8>
## x
## width
## quote
## justify
## escapes
## display
## utf8
##   <environment: namespace:utf8>
## name
## value
## nonnegative
##   <environment: namespace:utf8>
## name
## value
## nonnegative
##   <environment: namespace:utf8>
## name
## value
##   <environment: namespace:utf8>
## name
## value
## choices
##   <environment: namespace:utf8>
## name
## value
##   <environment: namespace:utf8>
## name
## value
## utf8
##   <environment: namespace:utf8>
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## ...
## class
##   <environment: namespace:pillar>
## x
## ...
## min_width
##   <environment: namespace:pillar>
## formatted
## ...
## width
## align
## min_width
## na_indent
##   <environment: namespace:pillar>
## col_widths_df
## tier_widths
## [1]
## [2]
## [3]
##   <environment: namespace:pillar>
## col_widths
## max_widths
## width
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:fansi>
## x
## start
## stop
## warn
## term.cap
## ctl
##   <environment: namespace:fansi>
## x
## start
## stop
## type
## round
## tabs.as.spaces
## tab.stops
## warn
## term.cap
## ctl
##   <environment: namespace:fansi>
## x
## start
## stop
## type.int
## round
## tabs.as.spaces
## tab.stops
## warn
## term.cap.int
## round.start
## round.stop
## x.len
## ctl.int
##   <environment: namespace:fansi>
## x
##   <environment: namespace:pillar>
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## use_brackets_if_no_color
##   <environment: namespace:pillar>
## x
##   <environment: namespace:pillar>
## x
## width
## align
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
## [1]
## [2]
## [3]
## [4]
## [5]
##   <environment: 0x000001dd486d9388>
## x
##   <environment: namespace:pillar>
## x
## width
## ...
##   <environment: namespace:pillar>
## x
##   <environment: namespace:rlang>
## from
## to
##   <environment: namespace:rlang>
## ...
##   <environment: namespace:pillar>
## x
## colonnade
## extra_cols
##   <environment: namespace:tibble>
## x
##   <environment: namespace:pillar>
## x
## ...
##   <environment: namespace:pillar>
## x
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
## squeezed_colonnade
## [1]
## [2]
## [3]
##   <environment: namespace:tibble>
## x
##   <environment: namespace:tibble>
## x
## extra_cols
##   <environment: namespace:pillar>
## x
## ...
## n
##   <environment: namespace:pillar>
## .x
## .f
##   <environment: namespace:pillar>
## .x
## .y
## .f
## ...
##   <environment: namespace:pillar>
## x
## title
## # A tibble: 6 x 5
##   from  to    source            database       layer     
##   <chr> <chr> <chr>             <chr>          <chr>     
## 1 SMAD1 TGFBI regnetwork_source regnetwork     regulatory
## 2 SMAD1 TGFBI Remap_5           Remap          regulatory
## 3 SMAD2 LAMC2 harmonizome_CHEA  harmonizome_gr regulatory
## 4 SMAD2 TGFBI harmonizome_CHEA  harmonizome_gr regulatory
## 5 SMAD2 TNC   harmonizome_CHEA  harmonizome_gr regulatory
## 6 SMAD2 TNC   regnetwork_source regnetwork     regulatory
```

Export the following to e.g. Cytoscape for exploration of the networks

``` r
output_path = ""
write_output = FALSE # change to TRUE for writing output

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"), active_signaling_network$gr %>% mutate(layer = "regulatory")) %>% write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
}

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl = bind_rows(
  tibble(gene = ligands_all, annotation = "ligand"),
  tibble(gene = targets_all, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_all,ligands_all)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
##   <environment: namespace:tibble>
## x
## length
##   <environment: namespace:dplyr>
## x
## y
## ...
##   <environment: namespace:dplyr>
## x
## y
## ...
non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
bind_rows(specific_annotation_tbl,non_specific_annotation_tbl) %>% write_tsv(paste0(output_path,"annotation_table.txt"))
}
```
