make_ligand_activity_target_lfc_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, prioritization_tbl_ligand_target, DE_table_receiver_processed_targets, lfc_cutoff, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()

  # ligand lfc
  ordered_ligands = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score) %>% distinct() %>% group_by(ligand) %>% summarise(ligand_score = max(ligand_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, ligand_score) %>% distinct()) %>% arrange(sender, ligand_score)
  ordered_ligands = ordered_ligands %>% mutate(ligand_ordered = factor(ligand, ordered = T, levels = ordered_ligands$ligand)) %>% distinct(ligand, ligand_ordered, niche) %>% rename(niche_prior = niche)
  plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligands)
  p1 = plot_data %>%
    ggplot(aes(sender, ligand_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches")
  max_lfc = abs(plot_data$ligand_score) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p_ligand_lfc = p1 + custom_scale_fill + ylab("Prioritized Ligands") + xlab("Ligand LogFC in Sender")

  # Target expression
  targets_oi = prioritization_tbl_ligand_target %>% filter(target_score >= lfc_cutoff) %>% filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% pull(target) %>% unique()

  ordered_targets = prioritization_tbl_ligand_target %>% filter(target %in% targets_oi) %>% select(niche, receiver, target, target_score) %>% distinct()  %>% arrange(receiver, -target_score)
  ordered_targets = ordered_targets %>% mutate(target_ordered = factor(target, ordered = T, levels = ordered_targets$target)) %>% distinct(target, target_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = DE_table_receiver_processed_targets %>% inner_join(ordered_targets) %>% filter(receiver %in% (prioritization_tbl_ligand_target$receiver %>% unique()))

  # p1 = plot_data %>% mutate(receiver = factor(receiver, levels = c("MoMac2","MoMac1","KCs"))) %>%
  p1 = plot_data %>%
    # ggplot(aes(target_ordered, receiver , color = target_expression_scaled_myeloid, size = target_fraction )) +
    ggplot(aes(target_ordered, receiver , fill = target_score)) +
    # geom_point() +
    geom_tile(color = "black") +
    # facet_grid(receiver~niche_prior, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x =  element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_blank()
    ) + labs(fill = "Target:\nmin LFC vs\nother niches") + xlab("Target LogFC in Receiver")
  max_exprs = abs(plot_data$target_score ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.41, 0.4850, 0.5, 0.5150, 0.6, 0.7, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_targets = p1 + custom_scale_fill


  # Ligand-Target heatmap
  active_ligand_target_links_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::select(ligand, target, ligand_target_weight ) %>% dplyr::rename(weight = ligand_target_weight )

  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.33
  }

  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)

  order_ligands_ = ordered_ligands$ligand_ordered %>% levels()
  order_targets_ = ordered_targets$target_ordered %>% levels()

  order_ligands = order_ligands_ %>% make.names()
  order_targets = order_targets_ %>% make.names()

  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

  if( length(setdiff(order_ligands, colnames(active_ligand_target_links))) > 0){
    removed_ligands = setdiff(order_ligands, colnames(active_ligand_target_links))
    new_lt_tibble = removed_ligands %>% lapply(function(ligand_oi){
      tibble(ligand = ligand_oi, target = order_targets, weight = 0)
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(ligand, weight)

    active_ligand_target_links = new_lt_tibble %>% select(-target) %>% data.frame() %>% as.matrix(ncol = length(removed_ligands)) %>% cbind(active_ligand_target_links)
  }

  if( length(setdiff(order_targets, rownames(active_ligand_target_links))) > 0){
    removed_targets = setdiff(order_targets, rownames(active_ligand_target_links))
    new_lt_tibble = removed_targets %>% lapply(function(target_oi){
      tibble(target = target_oi, ligand = order_ligands, weight = 0)
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(target, weight)

    active_ligand_target_links = new_lt_tibble %>% select(-ligand) %>% data.frame() %>% as.matrix(ncol = length(removed_targets)) %>% t() %>% rbind(active_ligand_target_links)
  }

  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {

    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }

  p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")



  # Ligand-Activity-Scaled

  order_receivers = prioritization_tbl_ligand_target %>% pull(receiver) %>% levels()

  # ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi)  %>% dplyr::select(ligand, niche, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(niche, activity_normalized)
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity_normalized)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Scaled Ligand activity", color = "purple",legend_position = "top", x_axis_position = "top", legend_title = "Scaled\nLigand\nActivity") + theme(legend.text = element_text(size = 9))

  # limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE))
  limits = c(-max(abs(vis_ligand_pearson), na.rm = TRUE), max(abs(vis_ligand_pearson), na.rm = TRUE))
  # print(limits)
  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.50, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill

  # Ligand-Activity
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson = p_ligand_pearson + custom_scale_fill


  # Combine the plots
  n_groups = ncol(vis_ligand_pearson)
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_senders = prioritization_tbl_ligand_receptor_filtered$sender %>% unique() %>% length()

  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc )), ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))

  if(is.null(heights)){
    heights = c(n_ligands, n_groups + 0.5)
  }
  if(is.null(widths)){
    widths = c(n_senders + 0.5, n_groups, n_groups, n_targets)
  }

  if(plot_legend == FALSE){
    design <- "SAaB
               ###C"
    combined_plot = patchwork::wrap_plots(S = p_ligand_lfc  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))

  } else {
    design <- "SAaB
               L##C"

    combined_plot = patchwork::wrap_plots(S = p_ligand_lfc  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }

}
make_ligand_zonation_activity_target_exprs_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()

  # ligand expression
  ordered_ligands = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score, ligand_score_zonation) %>% distinct() %>% group_by(ligand) %>% summarise(ligand_score = max(ligand_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, ligand_score, ligand_score_zonation) %>% distinct()) %>% arrange(sender, ligand_score)
  ordered_ligands = ordered_ligands %>% mutate(ligand_ordered = factor(ligand, ordered = T, levels = ordered_ligands$ligand)) %>% distinct(ligand, ligand_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = exprs_tbl_ligand %>% inner_join(prioritization_tbl_ligand_receptor %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score, ligand_score_zonation))   %>% inner_join(ordered_ligands) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver))



  p1 = plot_data  %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligands")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_ligands = p1 + custom_scale_fill
  p_ligands


  # ligand zonation
  p1 = plot_data  %>% filter(niche == "KC") %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Periportal-vs-Pericentral LFC") + xlab("Zonation LFC") + ylab("Prioritized Ligands")
  max_lfc = abs(plot_data$ligand_score_zonation ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_ligands_zonation = p1 + custom_scale_fill
  p_ligands_zonation

  # Target expression
  targets_oi = prioritization_tbl_ligand_target %>% filter(target_score >= lfc_cutoff) %>% filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% pull(target) %>% unique()

  ordered_targets = prioritization_tbl_ligand_target %>% filter(target %in% targets_oi) %>% select(niche, receiver, target, target_score) %>% distinct()  %>% arrange(receiver, -target_score)
  # ordered_targets = ordered_targets %>% select(-niche) %>% distinct() %>% mutate(niche = receiver)
  ordered_targets = ordered_targets %>% mutate(target_ordered = factor(target, ordered = T, levels = ordered_targets$target)) %>% distinct(target, target_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = exprs_tbl_target %>% inner_join(ordered_targets) %>% filter(receiver %in% (prioritization_tbl_ligand_target$receiver %>% unique()))
  plot_data = plot_data %>% group_by(target) %>% mutate(target_expression_scaled_myeloid = nichenetr::scaling_zscore(target_expression))

  # p1 = plot_data %>% mutate(receiver = factor(receiver, levels = c("MoMac2","MoMac1","KCs"))) %>%
  p1 = plot_data %>%
    # ggplot(aes(target_ordered, receiver , color = target_expression_scaled_myeloid, size = target_fraction )) +
    ggplot(aes(target_ordered, receiver , fill = target_expression_scaled_myeloid)) +
    # geom_point() +
    geom_tile(color = "black") +
    # facet_grid(receiver~niche_prior, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x =  element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_blank()
    ) + labs(fill = "Scaled Exprs Target") + xlab("Target Expression")
  max_exprs = abs(plot_data$target_expression_scaled_myeloid ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_targets = p1 + custom_scale_fill


  # Ligand-Target heatmap
  active_ligand_target_links_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::select(ligand, target, ligand_target_weight ) %>% dplyr::rename(weight = ligand_target_weight )

  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.33
  }

  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)

  order_ligands_ = ordered_ligands$ligand_ordered %>% levels()
  order_targets_ = ordered_targets$target_ordered %>% levels()

  order_ligands = order_ligands_ %>% make.names()
  order_targets = order_targets_ %>% make.names()

  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

  if( length(setdiff(order_ligands, colnames(active_ligand_target_links))) > 0){
    removed_ligands = setdiff(order_ligands, colnames(active_ligand_target_links))
    new_lt_tibble = removed_ligands %>% lapply(function(ligand_oi){
      tibble(ligand = ligand_oi, target = order_targets, weight = 0)
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(ligand, weight)

    active_ligand_target_links = new_lt_tibble %>% select(-target) %>% data.frame() %>% as.matrix(ncol = length(removed_ligands)) %>% cbind(active_ligand_target_links)
  }

  if( length(setdiff(order_targets, rownames(active_ligand_target_links))) > 0){
    removed_targets = setdiff(order_targets, rownames(active_ligand_target_links))
    new_lt_tibble = removed_targets %>% lapply(function(target_oi){
      tibble(target = target_oi, ligand = order_ligands, weight = 0)
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(target, weight)

    active_ligand_target_links = new_lt_tibble %>% select(-ligand) %>% data.frame() %>% as.matrix(ncol = length(removed_targets)) %>% t() %>% rbind(active_ligand_target_links)
  }

  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {

    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }

  p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")

  order_receivers = prioritization_tbl_ligand_target %>% pull(receiver) %>% levels()
  # Ligand-Activity-Scaled
  # ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi)  %>% dplyr::select(ligand, niche, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(niche, activity_normalized)
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity_normalized)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Scaled Ligand activity", color = "purple",legend_position = "top", x_axis_position = "top", legend_title = "Scaled\nLigand\nActivity") + theme(legend.text = element_text(size = 9))

  # limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE))
  limits = c(-max(abs(vis_ligand_pearson), na.rm = TRUE), max(abs(vis_ligand_pearson), na.rm = TRUE))

  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.50, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill

  # Ligand-Activity
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson = p_ligand_pearson + custom_scale_fill


  # Combine the plots
  n_groups = ncol(vis_ligand_pearson)
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_senders = prioritization_tbl_ligand_receptor_filtered$sender %>% unique() %>% length()

  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligands)), ggpubr::as_ggplot(ggpubr::get_legend(p_ligands_zonation)), ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))



  if(is.null(heights)){
    heights = c(n_ligands, n_groups + 0.5)
  }
  if(is.null(widths)){
    widths = c(n_senders + 0.5, 3, n_groups, n_groups, n_targets)
  }

  if(plot_legend == FALSE){
    design <- "SZAaB
               ####C"
    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          Z = p_ligands_zonation + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))

  } else {
    design <- "SZAaB
               L###C"

    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          Z = p_ligands_zonation + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }

}

make_ligand_activity_target_exprs_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor_filtered, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()

  # ligand expression
  ordered_ligands = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand %in% best_upstream_ligands) %>% select(niche, sender, ligand, ligand_score) %>% distinct() %>% group_by(ligand) %>% summarise(ligand_score = max(ligand_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, ligand_score) %>% distinct()) %>% arrange(sender, ligand_score)
  ordered_ligands = ordered_ligands %>% mutate(ligand_ordered = factor(ligand, ordered = T, levels = ordered_ligands$ligand)) %>% distinct(ligand, ligand_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligands) %>% filter(sender %in% (prioritization_tbl_ligand_receptor_filtered$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver, niche))

  p1 = plot_data  %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligands")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_ligands = p1 + custom_scale_fill
  p_ligands




  # Target expression
  targets_oi = prioritization_tbl_ligand_target %>% filter(target_score >= lfc_cutoff) %>% filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% pull(target) %>% unique()

  ordered_targets = prioritization_tbl_ligand_target %>% filter(target %in% targets_oi) %>% select(niche, receiver, target, target_score) %>% distinct()  %>% arrange(receiver, -target_score)

  # if duplicated: eg MoMac1-vs-MoMac1_CV
  ordered_targets = ordered_targets %>% select(-niche) %>% distinct() %>% mutate(niche = receiver)
  ordered_targets = ordered_targets %>% mutate(target_ordered = factor(target, ordered = T, levels = ordered_targets$target)) %>% distinct(target, target_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = exprs_tbl_target %>% inner_join(ordered_targets) %>% filter(receiver %in% (prioritization_tbl_ligand_target$receiver %>% unique()))
  plot_data = plot_data %>% group_by(target) %>% mutate(target_expression_scaled_myeloid = nichenetr::scaling_zscore(target_expression))

  # p1 = plot_data %>% mutate(receiver = factor(receiver, levels = c("MoMac2","MoMac1","KCs"))) %>%
  p1 = plot_data %>%

    # ggplot(aes(target_ordered, receiver , color = target_expression_scaled_myeloid, size = target_fraction )) +
    ggplot(aes(target_ordered, receiver , fill = target_expression_scaled_myeloid)) +
    # geom_point() +
    geom_tile(color = "black") +
    # facet_grid(receiver~niche_prior, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x =  element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_blank()
    ) + labs(fill = "Scaled Exprs Target") + xlab("Target Expression")
  max_exprs = abs(plot_data$target_expression_scaled_myeloid ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_targets = p1 + custom_scale_fill
  p_targets

  # Ligand-Target heatmap
  active_ligand_target_links_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% dplyr::select(ligand, target, ligand_target_weight ) %>% dplyr::rename(weight = ligand_target_weight )

  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.33
  }

  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)

  order_ligands_ = ordered_ligands$ligand_ordered %>% levels()
  order_targets_ = ordered_targets$target_ordered %>% levels()

  order_ligands = order_ligands_ %>% make.names()
  order_targets = order_targets_ %>% make.names()

  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

  if( length(setdiff(order_ligands, colnames(active_ligand_target_links))) > 0){
    removed_ligands = setdiff(order_ligands, colnames(active_ligand_target_links))
    new_lt_tibble = removed_ligands %>% lapply(function(ligand_oi){
      tibble(ligand = ligand_oi, target = order_targets, weight = 0)
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(ligand, weight)

    active_ligand_target_links = new_lt_tibble %>% select(-target) %>% data.frame() %>% as.matrix(ncol = length(removed_ligands)) %>% cbind(active_ligand_target_links)
  }

  if( length(setdiff(order_targets, rownames(active_ligand_target_links))) > 0){
    removed_targets = setdiff(order_targets, rownames(active_ligand_target_links))
    new_lt_tibble = removed_targets %>% lapply(function(target_oi){
      tibble(target = target_oi, ligand = order_ligands, weight = 0)
    }) %>% bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(target, weight)

    active_ligand_target_links = new_lt_tibble %>% select(-ligand) %>% data.frame() %>% as.matrix(ncol = length(removed_targets)) %>% t() %>% rbind(active_ligand_target_links)
  }

  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {

    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }

  p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  p_ligand_target_network

  if(is.null( prioritization_tbl_ligand_target %>% pull(receiver) %>% levels())){
    order_receivers = prioritization_tbl_ligand_target %>% pull(receiver) %>% unique()
  } else{
    order_receivers = prioritization_tbl_ligand_target %>% pull(receiver) %>% levels()
  }

  # Ligand-Activity-Scaled
  # ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi)  %>% dplyr::select(ligand, niche, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(niche, activity_normalized)
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity_normalized)
  # print(ligand_pearson_df)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers  %>% make.names()] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Scaled Ligand activity", color = "purple",legend_position = "top", x_axis_position = "top", legend_title = "Scaled\nLigand\nActivity") + theme(legend.text = element_text(size = 9))

  # limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE))
  limits = c(-max(abs(vis_ligand_pearson), na.rm = TRUE), max(abs(vis_ligand_pearson), na.rm = TRUE))
  # print(limits)

  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.50, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill
  p_ligand_pearson_scaled

  # Ligand-Activity
  ligand_pearson_df = prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands)  %>% dplyr::select(ligand, receiver, activity) %>% dplyr::distinct() %>% tidyr::spread(receiver, activity)
  ligand_pearson_matrix = ligand_pearson_df %>% dplyr::select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(ligand_pearson_df$ligand)
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands %>% generics::intersect(rownames(ligand_pearson_matrix)), order_receivers %>% make.names()] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Ligand\nActivity") + theme(legend.text = element_text(size = 9))
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE)))
  p_ligand_pearson = p_ligand_pearson + custom_scale_fill
  p_ligand_pearson

  # Combine the plots
  n_groups = ncol(vis_ligand_pearson)
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_senders = prioritization_tbl_ligand_receptor_filtered$sender %>% unique() %>% length()

  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligands)), ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))

  if(is.null(heights)){
    heights = c(n_ligands, n_groups + 0.5)
  }
  if(is.null(widths)){
    widths = c(n_senders + 0.5, n_groups, n_groups, n_targets)
  }

  if(plot_legend == FALSE){
    design <- "SAaB
               ###C"
    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))

  } else {
    design <- "SAaB
               L##C"

    combined_plot = patchwork::wrap_plots(S = p_ligands  + theme(legend.position = "none", axis.ticks = element_blank()),
                                          A = p_ligand_pearson_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
                                          a = p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none"),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }

}

make_ligand_receptor_lfc_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = TRUE, heights = NULL, widths = NULL){

  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()

  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score)
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()

  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score)
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligand_receptors)
  p_lig_lfc = plot_data %>%
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in Sender")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p_lig_lfc = p_lig_lfc + custom_scale_fill
  p_lig_lfc

  p_rec_lfc = plot_data %>%
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches")  + xlab("Receptor LFC\n in Receiver")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p_rec_lfc = p_rec_lfc + custom_scale_fill
  p_rec_lfc

  design = "A#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, B = p_rec_lfc, nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}


make_ligand_receptor_exprs_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, exprs_tbl_ligand, exprs_tbl_receiver, plot_legend = TRUE, heights = NULL, widths = NULL){

  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()

  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score)
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()

  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score)
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligand_receptors) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver))

  p1 = plot_data  %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligand-Receptor pairs")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_lig_lfc = p1 + custom_scale_fill

  plot_data = exprs_tbl_receptor %>% inner_join(ordered_ligand_receptors) %>% filter(receiver %in% (prioritization_tbl_ligand_receptor$receiver %>% unique()))
  plot_data = plot_data %>% group_by(receptor) %>% mutate(receptor_expression_scaled_sender = nichenetr::scaling_zscore(receptor_expression)) %>% inner_join(prioritization_tbl_ligand_receptor %>% distinct(sender, receiver))

  p1 = plot_data  %>%
    # ggplot(aes(receptor_ordered, sender , color = receptor_expression_scaled_myeloid, size = receptor_fraction )) +
    ggplot(aes(receiver,ligand_receptor_ordered , fill = receptor_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs receptor") + xlab("Receptor Expression")
  max_exprs = abs(plot_data$receptor_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_rec_lfc = p1 + custom_scale_fill
  p_rec_lfc

  design = "A#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, B = p_rec_lfc, nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}
make_ligand_receptor_exprs_zonation_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, exprs_tbl_ligand, exprs_tbl_receiver, plot_legend = TRUE, heights = NULL, widths = NULL){

  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()

  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score)
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()

  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score)
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligand_receptors) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% distinct(sender, receiver))

  ## ligand figure
  p1 = plot_data  %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligand-Receptor pairs")
  max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_lig_lfc = p1 + custom_scale_fill

  plot_data = exprs_tbl_ligand %>% inner_join(ordered_ligand_receptors) %>% filter(sender %in% (prioritization_tbl_ligand_receptor$sender %>% unique()))
  plot_data = plot_data %>% group_by(ligand) %>% mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% inner_join(prioritization_tbl_ligand_receptor %>% distinct(sender, receiver, ligand, ligand_score_zonation))


  # ligand zonation
  p1 = plot_data  %>% filter(niche == "KC") %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Sender\nPeriportal-vs-Pericentral\nLigand LFC") + xlab("Sender Zonation\nLFC ligand")
  max_lfc = abs(plot_data$ligand_score_zonation ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_ligands_zonation = p1 + custom_scale_fill
  # p_ligands_zonation

  ## receptor fugure
  plot_data = exprs_tbl_receptor %>% inner_join(ordered_ligand_receptors) %>% filter(receiver %in% (prioritization_tbl_ligand_receptor$receiver %>% unique()))
  plot_data = plot_data %>% group_by(receptor) %>% mutate(receptor_expression_scaled_sender = nichenetr::scaling_zscore(receptor_expression)) %>% inner_join(prioritization_tbl_ligand_receptor %>% distinct(sender, receiver))

  p1 = plot_data  %>%
    # ggplot(aes(receptor_ordered, sender , color = receptor_expression_scaled_myeloid, size = receptor_fraction )) +
    ggplot(aes(receiver,ligand_receptor_ordered , fill = receptor_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs receptor") + xlab("Receptor Expression")
  max_exprs = abs(plot_data$receptor_expression_scaled_sender ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_rec_lfc = p1 + custom_scale_fill
  # p_rec_lfc

  design = "AZ#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, Z = p_ligands_zonation + ylab(""), B = p_rec_lfc + ylab(""), nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 3, 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}


make_ligand_receptor_lfc_zonation_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_receptor_filtered, plot_legend = TRUE, heights = NULL, widths = NULL){

  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()

  ordered_ligand_receptors = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct() %>% group_by(ligand_receptor) %>% summarise(prioritization_score = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score)
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor_filtered %>% filter(ligand_receptor %in% filtered_ligand_receptors) %>% select(niche, sender, ligand, prioritization_score) %>% distinct() %>% group_by(ligand) %>% summarise(prioritization_score_ligand = max(prioritization_score)) %>% inner_join(prioritization_tbl_ligand_receptor_filtered %>% select(niche, sender, ligand, prioritization_score) %>% distinct()) %>% arrange(sender, prioritization_score_ligand) %>% distinct()

  ordered_ligand_receptors = ordered_ligand_receptors %>% inner_join(ordered_ligand_receptors_max_ligand_score) %>% arrange(sender, prioritization_score_ligand, prioritization_score)
  ordered_ligand_receptors = ordered_ligand_receptors %>% mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, levels = ordered_ligand_receptors$ligand_receptor)) %>% distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% rename(niche_prior = niche)

  plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligand_receptors)
  p_lig_lfc = plot_data %>%
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in sender")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p_lig_lfc = p_lig_lfc + custom_scale_fill
  # p_lig_lfc

  # ligand zonation

  senders_zonated = plot_data %>% filter(ligand_score_zonation != 0) %>% pull(sender) %>% unique()

  p1 = plot_data %>% filter(sender %in% senders_zonated) %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Sender\nRegion-oi-vs-Other\nLigand LFC") + xlab("Sender Spatial\nLFC ligand")
  max_lfc = abs(plot_data$ligand_score_zonation ) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_ligands_zonation = p1 + custom_scale_fill


  p_rec_lfc = plot_data %>%
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches") + xlab("Receptor LFC\n in Receiver")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p_rec_lfc = p_rec_lfc + custom_scale_fill
  # p_rec_lfc

  design = "AZ#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, Z = p_ligands_zonation + ylab(""), B = p_rec_lfc + ylab(""), nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 3, 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
}
############################## ############################## ############################## ##############################
############################## smaller functions############################## ##############################
############################## ############################## ############################## ##############################

ligand_lfc_plot = function(plot_data, max_lfc){
  p_lig_lfc = plot_data %>%
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand:\nmin LFC vs\nother niches")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in sender")

  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  p_lig_lfc = p_lig_lfc + custom_scale_fill

  return(p_lig_lfc)

}

ligand_lfc_zonation_plot = function(plot_data, max_lfc){
  # ligand zonation
  p1 = plot_data   %>% filter(niche == "KC") %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_receptor_ordered , fill = ligand_score_zonation)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Sender\nPeriportal-vs-Pericentral\nLigand LFC") + xlab("Sender Zonation\nLFC ligand")
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.3, 0.466, 0.5, 0.533, 0.6, 1),  limits = c(-1*max_lfc, max_lfc))
  # custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

  p_ligands_zonation = p1 + custom_scale_fill

  return(p_ligands_zonation)

}

receptor_lfc_plot = function(plot_data, max_lfc){
  p_rec_lfc = plot_data %>%
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    # facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches") + xlab("Receptor LFC\n in Receiver")
  # max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))

  p_rec_lfc = p_rec_lfc + custom_scale_fill
  # p_rec_lfc
  return(p_rec_lfc)
}
normalized_activity_plot = function(plot_data, max_activity, min_activity){
  p_lig_lfc = plot_data %>%
    ggplot(aes(receiver, ligand_receptor_ordered, fill = activity_normalized)) +
    geom_tile(color = "black") +
    # facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Normalized Ligand activity\nin Receiver")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Normalized Ligand activity\nin Receiver")

  limits = c(-max_normalized_activity, max_activity)
  custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),values = c(0, 0.49, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  limits = limits)
  p_activity = p_lig_lfc +custom_scale_fill

  return(p_activity)

}
activity_plot = function(plot_data, max_activity, min_activity){
  p_lig_lfc = plot_data %>%
    ggplot(aes(receiver, ligand_receptor_ordered, fill = activity)) +
    geom_tile(color = "black") +
    # facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Ligand activity\nin Receiver")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand activity\nin Receiver")

  limits = c(min_activity, max_activity)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  limits = limits) # non-scaled
  p_activity = p_lig_lfc +custom_scale_fill

  return(p_activity)

}
prioritization_score_plot = function(plot_data){
  p_score = plot_data %>%
    ggplot(aes(scoretype   , ligand_receptor_ordered, fill = score)) +
    geom_tile(color = "black") +
    # facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "LR Prioritization Scores")  + ylab("Prioritized Ligand-Receptor pairs") + xlab("Prioritization Scores\nLR pair")

  limits = c(0, 1.01)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "YlGn")),values = c(0, 0.50, 0.625, 0.75, 0.825, 0.885, 0.945, 1),  limits = limits) # non-scaled
  p_score = p_score +custom_scale_fill

  return(p_score)

}

make_circos_lr= function(prioritized_tbl_oi, colors_sender, colors_receiver, cutoff, scale, transparency = NULL, circos_type, border = TRUE){

  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")

  # Link each cell type to a color
  grid_col_ligand = colors_sender
  # names(grid_col_ligand) = prioritized_tbl_oi$sender %>% unique() %>% sort()

  grid_col_receptor = colors_receiver
  # names(grid_col_receptor) = prioritized_tbl_oi$receiver %>% unique() %>% sort()

  grid_col_tbl_ligand = tibble::tibble(sender = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
  grid_col_tbl_receptor = tibble::tibble(receiver = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

  # Make the plot for condition of interest - title of the plot
  circos_links_oi = prioritized_tbl_oi

  # deal with duplicated sector names
  # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
  # only do it with duplicated ones!
  circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score)

  #circos_links = circos_links_oi %>% dplyr::rename(weight = prioritization_score) %>% dplyr::mutate(ligand = paste(sender, ligand, sep = "_"), receptor = paste(receptor, receiver, sep = "_"))

  df = circos_links %>% mutate(ligand_receptor_sender_receiver = paste0(sender, receiver, ligand_receptor))

  ligand.uni = unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i = df[df$ligand == ligand.uni[i], ]
    sender.uni = unique(df.i$sender)
    for (j in 1:length(sender.uni)) {
      df.i.j = df.i[df.i$sender == sender.uni[j], ]
      df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$ligand_receptor_sender_receiver %in% df.i.j$ligand_receptor_sender_receiver] = df.i.j$ligand
    }
  }
  receptor.uni = unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i = df[df$receptor == receptor.uni[i], ]
    receiver.uni = unique(df.i$receiver)
    for (j in 1:length(receiver.uni)) {
      df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
      df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$ligand_receptor_sender_receiver %in% df.i.j$ligand_receptor_sender_receiver] = df.i.j$receptor
    }
  }

  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))

  # print(intersecting_ligands_receptors)

  while(length(intersecting_ligands_receptors) > 0){
    df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
    df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
    df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(receptor, " ", sep = ""))
    df = dplyr::bind_rows(df_unique, df_duplicated)
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  }

  circos_links = df

  # Link ligands/Receptors to the colors of senders/receivers
  circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
  links_circle = circos_links %>% dplyr::distinct(ligand,receptor, weight)
  ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
  receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
  grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
  grid_col =c(grid_ligand_color,grid_receptor_color)

  # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
  # transparency = circos_links %>% mutate(old_weight = weight) %>% mutate_cond(old_weight < cutoff, weight = 0) %>% mutate_cond(old_weight >= cutoff, weight = 2)  %>% mutate_cond(old_weight > cutoff+0.15, weight = 3) %>% mutate_cond(old_weight > cutoff+0.2, weight = 3.5)  %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
  if(is.null(transparency)) {
    # transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    # transparency = circos_links %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency
    transparency = circos_links %>% dplyr::mutate(transparency = 1-weight) %>% .$transparency

  }

  # Define order of the ligands and receptors and the gaps
  ligand_order = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
    # ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(weight) %>% dplyr::distinct(ligand)

  }) %>% unlist()


  receptor_order = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    circos_links_n = circos_links %>% dplyr::filter(receiver == receiver_oi) %>% group_by(receptor) %>% count() %>% ungroup()
    receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>% inner_join(circos_links_n) %>% dplyr::arrange(-n, ligand) %>% dplyr::distinct(receptor)
    # receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(weight) %>% dplyr::distinct(receptor)

  }) %>% unlist()

  # receptor_order_last = c("CSF1R","NOTCH2")
  # receptor_order_first = c("SLC40A1","BMPR1A","BMPR2","ACVRL1","LRP6")
  # receptor_order = c(receptor_order_first, setdiff(receptor_order, c(receptor_order_last,receptor_order_first)), receptor_order_last)


  order = c(ligand_order,receptor_order)
  # print(length(order))

  width_same_cell_same_ligand_type = 1
  width_different_cell = 5
  width_ligand_receptor = 15
  width_same_cell_same_receptor_type = 1

  sender_gaps = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  sender_gaps = sender_gaps[-length(sender_gaps)]

  receiver_gaps = prioritized_tbl_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  receiver_gaps = receiver_gaps[-length(receiver_gaps)]

  gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)

  # print(length(gaps))
  # print(length(union(circos_links$ligand, circos_links$receptor) %>% unique()))
  if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
  }



  # links_circle$weight[links_circle$weight == 0] = 0.01
  circos.clear()
  circos.par(gap.degree = gaps)

  if(circos_type == "arrow"){
    chordDiagram(links_circle,
                 directional = 1,
                 order=order,
                 link.sort = TRUE,
                 link.decreasing = FALSE,
                 grid.col = grid_col,
                 transparency = transparency,
                 diffHeight = 0.0075,
                 direction.type = c("diffHeight", "arrows"),
                 link.visible = links_circle$weight >= cutoff,
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.25),
                 # grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                 link.arr.col = arr.col, link.arr.length = 0.4, link.arr.lwd = 5, link.arr.width = 0.2,
                 reduce = 0
                 , scale = scale ### TRUE: width of the sectors does not depend on the nr of links
    )
  } else {
    if(border == TRUE) {
      chordDiagram(links_circle,
                   directional = 1,
                   order=order,
                   link.sort = TRUE,
                   link.decreasing = FALSE,
                   grid.col = grid_col,
                   transparency = transparency,
                   diffHeight = 0.0075,
                   direction.type = c("diffHeight", "arrows"),
                   link.visible = links_circle$weight >= cutoff,
                   annotationTrack = "grid",
                   preAllocateTracks = list(track.height = 0.25),
                   grid.border = "gray35", link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1, link.border="gray35",
                   # link.arr.col = arr.col, link.arr.length = 0.4, link.arr.lwd = 5, link.arr.width = 0.2,
                   reduce = 0
                   , scale = scale ### TRUE: width of the sectors does not depend on the nr of links
      )
    } else {
      chordDiagram(links_circle,
                   directional = 1,
                   order=order,
                   link.sort = TRUE,
                   link.decreasing = FALSE,
                   grid.col = grid_col,
                   transparency = transparency,
                   diffHeight = 0.0075,
                   direction.type = c("diffHeight", "arrows"),
                   link.visible = links_circle$weight >= cutoff,
                   annotationTrack = "grid",
                   preAllocateTracks = list(track.height = 0.25),
                   link.arr.length = 0.05, link.arr.type = "big.arrow",
                   # link.arr.col = arr.col, link.arr.length = 0.4, link.arr.lwd = 5, link.arr.width = 0.2,
                   reduce = 0
                   , scale = scale ### TRUE: width of the sectors does not depend on the nr of links
      )
    }

  }




  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
  }, bg.border = NA) #

  p_circos = recordPlot()

  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(grid_col_receptor, grid_col_ligand)
  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = grid_col_receptor[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))

  legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = grid_col_ligand[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))

  p_legend = grDevices::recordPlot()


  return(list(p_circos = p_circos, p_legend = p_legend))

}


