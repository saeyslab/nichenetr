context("Target gene prediction functions for application")
test_that("Gene-cluster assignments can be converted to settings", {

  genes_clusters = c("TGFB1" = 1,"TGFB2" = 1,"TGFB3" = 2)
  cluster_settings = lapply(seq(length(unique(genes_clusters))), convert_cluster_to_settings, cluster_vector = genes_clusters, setting_name = "example", setting_from = "BMP2")

  expect_type(cluster_settings,"list")
  expect_type(cluster_settings[[1]]$response,"logical")
  expect_type(cluster_settings[[1]]$from,"character")
  expect_type(cluster_settings[[1]]$name,"character")
  expect_type(cluster_settings[[2]]$response,"logical")

})

context("Ligand activity prediction for application")
test_that("Ligand activities can be predicted well", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","IL4")
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
  potential_ligands = c("TNF","IL4")
  geneset = c("VCAM1","SOCS3", "IRF1","LTBR","STAT4")
  background_expressed_genes = c("VCAM1","SOCS3","IRF1","ICAM1","ID1","ID2","ID3","LTBR",'STAT4')
  ligand_activities = predict_ligand_activities(geneset = geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

  expect_type(ligand_activities,"list")
  expect_type(ligand_activities$test_ligand,"character")
  expect_type(ligand_activities$pearson,"double")

  active_ligand_target_links_df = potential_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows()
  expect_type(active_ligand_target_links_df,"list")

  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.001)
  expect_type(active_ligand_target_links,"double")

})
