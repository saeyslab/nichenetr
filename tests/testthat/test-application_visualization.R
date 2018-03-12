context("Visualization functions for application")
test_that("ligand-target signaling paths can be visualized", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","BMP2",c("IL4","IL13"))
  ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
  all_ligands = c("BMP2")
  all_targets = c("HEY1")
  k = 2
  ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,k,weighted_networks)
  data_source_info_network = infer_supporting_datasources(ligand_target_signaling_list, lr_network, sig_network, gr_network)
  graph = diagrammer_format_signaling_graph(ligand_target_signaling_list, all_ligands,all_targets)

  expect_type(ligand_target_signaling_list,"list")
  expect_type(data_source_info_network,"list")
  expect_type(graph,"list")
  expect_equal(length(ligand_target_signaling_list),2)
})
