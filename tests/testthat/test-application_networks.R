context("Network extraction functions for application")
test_that("active ligand-receptor, signaling and gene regulatory networks can be constructed", {
  expression_vector_sender = rnorm(n = 10000, mean = 6, sd = 3)
  expression_vector_receiver = rnorm(n = 10000, mean = 6, sd = 3)
  names(expression_vector_sender) = sample(x = geneinfo_human$symbol,size = 10000,replace = FALSE)
  names(expression_vector_receiver) = sample(x = geneinfo_human$symbol,size = 10000,replace = FALSE)

  sender_receiver_lr_network = get_active_ligand_receptor_network(expression_vector_sender,expression_vector_receiver,lr_network,expression_cutoff_sender = 0, expression_cutoff_receiver = 4)
  receiver_sig_network = get_active_signaling_network(expression_vector_receiver,sig_network,expression_cutoff_receiver = 4)
  receiver_gr_network = get_active_regulatory_network(expression_vector_receiver,gr_network,expression_cutoff_receiver = 4)

  weighted_networks_receiver = construct_weighted_networks(sender_receiver_lr_network, receiver_sig_network, receiver_gr_network, source_weights_df)

  expect_type(sender_receiver_lr_network,"list")
  expect_gt(nrow(sender_receiver_lr_network),0)

  expect_type(receiver_sig_network,"list")
  expect_gt(nrow(receiver_sig_network),0)

  expect_type(receiver_gr_network,"list")
  expect_gt(nrow(receiver_gr_network),0)

  expect_type(weighted_networks_receiver,"list")
})
test_that("active ligand-target matrix and network can be constructed", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  setting = lapply(expression_settings_validation[1:3],convert_expression_settings_evaluation)
  ligands = extract_ligands_from_settings(setting)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands,ligands_as_cols = TRUE)

  # # temporary debugging solution:
  # ligands = colnames(ligand_target_matrix)
  # targets = rownames(ligand_target_matrix)
  # error_rate = 0.1
  # cutoff_method = "quantile"
  # fdr_method = "global"
  # ligands_position = "cols"
  # list_targets = lapply(ligands,get_target_genes_ligand_oi,ligand_target_matrix,cutoff_method = cutoff_method, fdr_method = fdr_method,error_rate = error_rate)
  # print(length(list_targets))
  # print(length(list_targets[[1]]))
  # names(list_targets) = ligands
  # ligand_target_matrix_discrete = list_targets %>% bind_rows() %>% as.matrix()
  #
  # # this means: ligand-target matrix discrete gets 3 rows and n targets columns??? weird
  #
  # print(dim(ligand_target_matrix_discrete))
  # expect_equal(dim(ligand_target_matrix_discrete), dim(ligand_target_matrix))
  #
  # rownames(ligand_target_matrix_discrete) = targets %>% make.names() # maybe change this again
  # expect_equal(dim(ligand_target_matrix_discrete), dim(ligand_target_matrix))

  ligand_target_matrix_discrete = make_discrete_ligand_target_matrix(ligand_target_matrix, cutoff_method = "quantile")
  active_lt = get_active_ligand_target_matrix(setting[[1]] %>% .$response, ligand_target_matrix_discrete)
  active_lt_df = get_active_ligand_target_df(setting[[1]] %>% .$response, ligand_target_matrix_discrete)

  expect_type(active_lt_df,"list")
  expect_gt(nrow(active_lt_df),0)

  expect_type(active_lt,"logical")
  active_lt = get_active_ligand_target_matrix(setting[[1]] %>% .$response, ligand_target_matrix)
  expect_type(active_lt,"double")

})

