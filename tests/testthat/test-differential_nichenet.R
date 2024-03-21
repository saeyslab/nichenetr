# context("Differential NicheNet")
# test_that("Differential NicheNet pipeline works", {
#   options(timeout = 3600)
#   ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#   ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
#
#   lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#   lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor)
#
#   seurat_object_lite = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
#   seurat_objs = list(Seurat::UpdateSeuratObject(seurat_object_lite))
#
#   # Test for v5
#   if (grepl("^5", packageVersion("Seurat"))){
#     seurat_objs[[2]] <- CreateSeuratObject(counts = Seurat::GetAssayData(seurat_object_lite, layer = "counts"),
#                                            meta.data = seurat_object_lite@meta.data) %>% SetIdent(value = .$celltype) %>%
#       NormalizeData()
#   }
#
#   seurat_object_lite@meta.data$celltype_aggregate = paste(seurat_object_lite@meta.data$celltype, seurat_object_lite@meta.data$aggregate,sep = "_") # user adaptation required on own dataset
#
#   celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
#   seurat_obj = SetIdent(seurat_object_lite, value = seurat_object_lite[[celltype_id, drop=TRUE]])
#
#   niches = list(
#     "LCMV_niche" = list(
#       "sender" = c("CD8 T_LCMV", "Mono_LCMV"),
#       "receiver" = c("CD8 T_LCMV")),
#     "SS_niche" = list(
#       "sender" = c("CD8 T_SS",  "Mono_SS"),
#       "receiver" = c("CD8 T_SS"))
#   )
#
#   assay_oi = "RNA" # other possibilities: RNA,...
#
#   DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
#   DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE
#
#   expression_pct = 0.10
#   DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
#   DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
#
#   specificity_score_LR_pairs = "min_lfc"
#   DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
#
#   include_spatial_info_sender = TRUE # if not spatial info to include: put this to false # user adaptation required on own dataset
#   include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
#   spatial_info = tibble(celltype_region_oi = "Mono_LCMV", celltype_other_region = "CD8 T_LCMV", niche =  "LCMV_niche", celltype_type = "sender") # user adaptation required on own dataset
#   specificity_score_spatial = "lfc"
#
#   if(include_spatial_info_sender == TRUE){
#     sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"), assay_oi = assay_oi)
#     sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
#
#     # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
#     sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
#     sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
#
#     sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
#
#   } else {
#     # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
#     sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
#     sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
#
#   }
#   if(include_spatial_info_receiver == TRUE){
#     receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"), assay_oi = assay_oi)
#     receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
#
#     # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
#     receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
#     receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
#
#     receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
#
#   } else {
#     # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
#     receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
#     receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
#   }
#
#   lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff.
#   specificity_score_targets = "min_lfc"
#
#   DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi)
#   DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
#
#   background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
#   geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
#   geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
#
#   # Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
#   # If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
#   geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
#   geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
#
#   length(geneset_niche1)
#   length(geneset_niche2)
#
#   top_n_target = 250
#
#   niche_geneset_list = list(
#     "LCMV_niche" = list(
#       "receiver" = niches[[1]]$receiver,
#       "geneset" = geneset_niche1,
#       "background" = background),
#     "SS_niche" = list(
#       "receiver" = niches[[2]]$receiver,
#       "geneset" = geneset_niche2 ,
#       "background" = background)
#   )
#
#   ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
#
#
#   features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
#
#   dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
#   exprs_tbl = dotplot$data %>% as_tibble()
#   exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
#     mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
#
#   exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction)
#   exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
#   exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
#
#   exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
#
#   exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
#
#   exprs_sender_receiver = lr_network %>%
#     inner_join(exprs_tbl_ligand, by = c("ligand")) %>%
#     inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
#
#   ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction) %>% distinct() %>% ungroup()
#
#   prioritizing_weights = c("scaled_ligand_score" = 5,
#                            "scaled_ligand_expression_scaled" = 1,
#                            "ligand_fraction" = 1,
#                            "scaled_ligand_score_spatial" = 2,
#                            "scaled_receptor_score" = 0.5,
#                            "scaled_receptor_expression_scaled" = 0.5,
#                            "receptor_fraction" = 1,
#                            "ligand_scaled_receptor_expression_fraction" = 1,
#                            "scaled_receptor_score_spatial" = 0,
#                            "scaled_activity" = 0,
#                            "scaled_activity_normalized" = 1)
#
#   output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
#                 ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
#   prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
#
#   prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
#   prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
#
#   prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
#   prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[2]]$receiver) %>% head(10)
#
#   expect_type(prioritization_tables,"list")
#
#   top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
#   top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
#
#   ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50
#
#
#   receiver_oi = "CD8 T_LCMV"
#
#   filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% pull(ligand) %>% unique()
#
#   prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup()
#
#   lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
#   lfc_plot
#
#   lfc_plot_spatial = make_ligand_receptor_lfc_spatial_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, ligand_spatial = include_spatial_info_sender, receptor_spatial = include_spatial_info_receiver, plot_legend = FALSE, heights = NULL, widths = NULL)
#   lfc_plot_spatial
#
#   exprs_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
#   exprs_plot$combined_plot
#
#   colors_sender = c("blue","red") %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
#   colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())
#
#   circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
#   expect_type(prioritized_tbl_oi,"list")
#
# })
