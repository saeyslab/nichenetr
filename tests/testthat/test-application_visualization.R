context("Visualization functions for application")
test_that("ligand-target signaling paths can be visualized", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","TGFB1",c("IL4","IL13"))
  ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.1, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
  all_ligands = c("TGFB1")
  all_targets = c("CCND1")

  k = 2

  ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,k,weighted_networks)
  data_source_info_network = infer_supporting_datasources(ligand_target_signaling_list, lr_network, sig_network, gr_network)
  graph = diagrammer_format_signaling_graph(ligand_target_signaling_list, all_ligands,all_targets)

  expect_type(ligand_target_signaling_list,"list")
  expect_type(data_source_info_network,"list")
  expect_type(graph,"list")
  expect_equal(length(ligand_target_signaling_list),2)

  all_receptors = "TGFBR1"
  ligand_target_signaling_list = get_ligand_signaling_path_with_receptor(ligand_tf_matrix,all_ligands,all_receptors, all_targets,k,weighted_networks)
  expect_type(ligand_target_signaling_list,"list")

})
test_that("heatmaps can be shown", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","TGFB1",c("IL4","IL13"))
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)

  p = make_heatmap_ggplot(ligand_target_matrix[1:50,] %>% t(), y_name = "ligand", x_name = "target")
  expect_equal(class(p)[2],"ggplot")

  p = make_threecolor_heatmap_ggplot(ligand_target_matrix[1:50,] %>% t(), y_name = "ligand", x_name = "target")
  expect_equal(class(p)[2],"ggplot")

  ligand_target_matrix_vis_genedirection = ligand_target_matrix %>% apply(1,scaling_modified_zscore) %>% .[,1:50]
  ligand_target_matrix_vis_genedirection[ligand_target_matrix_vis_genedirection < 2] = 0
  ligand_target_matrix_vis_genedirection[ligand_target_matrix_vis_genedirection != 0] = 1
  #'
  ligand_target_matrix_vis_liganddirection = ligand_target_matrix %>% apply(2,scaling_modified_zscore) %>% .[1:50, ] %>% t()
  ligand_target_matrix_vis_liganddirection[ligand_target_matrix_vis_liganddirection < 2] = 0
  ligand_target_matrix_vis_liganddirection[ligand_target_matrix_vis_liganddirection != 0] = 2
  #'
  bidirectional_ligand_target_matrix_vis = ligand_target_matrix_vis_genedirection + ligand_target_matrix_vis_liganddirection
  bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 0] = "none"
  bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 1] = "top-ligand"
  bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 2] = "top-target"
  bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 3] = "top"

  p = make_heatmap_bidir_lt_ggplot(bidirectional_ligand_target_matrix_vis, y_name = "ligand", x_name = "target")
  expect_equal(class(p)[2],"ggplot")


})
