context("Model construction functions")
#
test_that("Construct weighted networks: correct sum", {
  lrn_toy = tibble(from = "A", to = "B", source = "toy_lr")
  sign_toy = tibble(from = c("B","B"), to = c("C","C"), source = c("toy_sig1","toy_sig2"))
  grn_toy = tibble(from = "C", to = "D", source = "toy_grn")
  source_weights_df_toy = tibble(source = c("toy_lr","toy_sig1","toy_sig2","toy_grn"), weight = c(0.5,1,0.5,1))

  expected_wn = tibble(from = c("A","B","C"),
                       to = c("B","C","D"),
                       weight =c(0.5,1.5,1))

  expect_equal(construct_weighted_networks(lrn_toy, sign_toy, grn_toy, source_weights_df_toy) %>% bind_rows(), expected_wn)
})
test_that("Construct weighted networks: input weights not higher than 1", {
  lrn_toy = tibble(from = "A", to = "B", source = "toy_lr")
  sign_toy = tibble(from = c("B","B"), to = c("C","C"), source = c("toy_sig1","toy_sig2"))
  grn_toy = tibble(from = "C", to = "D", source = "toy_grn")
  source_weights_df_toy = tibble(source = c("toy_lr","toy_sig1","toy_sig2","toy_grn"), weight = c(0.5,2,0.5,1))

  expect_error(construct_weighted_networks(lrn_toy, sign_toy, grn_toy, source_weights_df_toy) %>% bind_rows())
})
test_that("Construct weighted networks: n_output networks 2 or 3", {
  lrn_toy = tibble(from = "A", to = "B", source = "toy_lr")
  sign_toy = tibble(from = c("B","B"), to = c("C","C"), source = c("toy_sig1","toy_sig2"))
  grn_toy = tibble(from = "C", to = "D", source = "toy_grn")
  source_weights_df_toy = tibble(source = c("toy_lr","toy_sig1","toy_sig2","toy_grn"), weight = c(0.5,1,0.5,1))

  expect_error(construct_weighted_networks(lrn_toy, sign_toy, grn_toy, source_weights_df_toy,1) %>% bind_rows())
  expect_error(construct_weighted_networks(lrn_toy, sign_toy, grn_toy, source_weights_df_toy,4) %>% bind_rows())

})
test_that("Add new data sources: common input base network", {
  lr_toy = tibble(from = "A", to = "B", source = "toy")
  new_lr_network = add_new_datasource(lr_toy, lr_network,1,source_weights_df)
  expect_equal(new_lr_network$network, bind_rows(lr_network, lr_toy))
  expect_equal(new_lr_network$source_weights_df %>% filter(source == "toy") %>% .$weight, 1)

})
test_that("Add new data sources: NULL input base network", {
  lr_toy = tibble(from = "A", to = "B", source = "toy")
  output = add_new_datasource(new_source = lr_toy, network = NULL,new_weight = 1,source_weights_df = source_weights_df)
  expect_equal(output$network, lr_toy)
  expect_equal(output$source_weights_df %>% filter(source == "toy") %>% .$weight, 1)

})
test_that("Add new data sources: weight higher than 0", {
  lr_toy = tibble(from = "A", to = "B", source = "toy")
  expect_error(add_new_datasource(new_source = lr_toy, network = NULL,new_weight = 1.2, source_weights_df = source_weights_df))
})
test_that("Correct application hub correction: gr_hub", {
  lrn_toy = tibble(from = "A", to = "B", source = "toy_lr")
  sign_toy = tibble(from = c("B","B"), to = c("C","C"), source = c("toy_sig1","toy_sig2"))
  grn_toy = tibble(from = c("C","E"), to = c("D","D"), source = "toy_grn")

  source_weights_df_toy = tibble(source = c("toy_lr","toy_sig1","toy_sig2","toy_grn"), weight = c(0.5,1,0.5,1))
  w = construct_weighted_networks(lrn_toy, sign_toy, grn_toy, source_weights_df_toy)
  expected_wn = tibble(from = c("A","B","C","E"),
                       to = c("B","C","D","D"),
                       weight =c(0.5,1.5,0.5,0.5))
  expect_equal(apply_hub_corrections(w,0,1) %>% bind_rows(), expected_wn)
})
test_that("Correct application hub correction: lr_sig_hub", {
  lrn_toy = tibble(from = "A", to = "B", source = "toy_lr")
  sign_toy = tibble(from = c("B","B","E"), to = c("C","C","C"), source = c("toy_sig1","toy_sig2","toy_sig2"))
  grn_toy = tibble(from = c("C"), to = c("D"), source = "toy_grn")

  source_weights_df_toy = tibble(source = c("toy_lr","toy_sig1","toy_sig2","toy_grn"), weight = c(0.5,1,1,1))
  w = construct_weighted_networks(lrn_toy, sign_toy, grn_toy, source_weights_df_toy)
  expected_wn2 = tibble(from = c("A","B","E","C"),
                        to = c("B","C","C","D"),
                        weight =c(0.5,1,0.5,1))
  expect_equal(apply_hub_corrections(w,1,0) %>% bind_rows(), expected_wn2)
})
test_that("Construct ligand_to_tf_matrix: no error", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("BMP2",c("IL4","IL13"))
  expect_type(construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5),"double")
  expect_type(construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "SPL", damping_factor = 0.5),"double")
  expect_type(construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "direct", damping_factor = 0.5),"double")

})
test_that("Construct ligand_to_target_matrix: no error", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("BMP2",c("IL4","IL13"))
  expect_type(construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE),"double")
  expect_type(construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "SPL", damping_factor = 0.5, secondary_targets = FALSE),"double")
  expect_type(construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "direct", damping_factor = 0.5, secondary_targets = FALSE),"double")
  expect_type(construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, secondary_targets = TRUE),"double")
  expect_type(construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, secondary_targets = TRUE, remove_direct_links = "ligand"),"double")
  expect_type(construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, secondary_targets = TRUE, remove_direct_links = "ligand-receptor"),"double")
})
test_that("Construct tf_to_target_matrix: no error", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  expect_type(construct_tf_target_matrix(weighted_networks, standalone_output = TRUE),"S4")
  expect_type(construct_tf_target_matrix(weighted_networks, standalone_output = FALSE),"S4")
})
test_that("Correct PPR-ligand-target matrices for topolgy: no error", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("BMP2",c("IL4","IL13"))
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
  expect_type(correct_topology_ppr(ligand_target_matrix,weighted_networks),"double")
})
test_that("Convert probabilistic ligand-target to discrete: no error", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("BMP2",c("IL4","IL13"))
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
  expect_equal(dim(make_discrete_ligand_target_matrix(ligand_target_matrix, error_rate = 0.1, cutoff_method = "distribution", ligands_position = "cols")),dim(ligand_target_matrix))
})


