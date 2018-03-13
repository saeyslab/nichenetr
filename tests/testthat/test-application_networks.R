context("Network extraction functions for application")
test_that("active ligand-receptor, signaling and gene regulatory networks can be constructed", {
  library(Biobase)
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
  mousesymbol2humansymbol = mapper(geneinfo_human,"symbol","symbol_mouse")
  expression_vector_sender = Exprs_lsec[,Exprs_lsec$celltype == "LSEC_12h"] %>% apply(1,mean)
  expression_vector_receiver = Exprs_mono_kc[,Exprs_mono_kc$celltype == "BM_mono"] %>% apply(1,mean)
  names(expression_vector_sender) = names(expression_vector_sender) %>% mousesymbol2humansymbol[.]
  names(expression_vector_receiver) = names(expression_vector_receiver) %>% mousesymbol2humansymbol[.]
  expression_vector_sender = expression_vector_sender %>% .[!is.na(names(.))]
  expression_vector_receiver = expression_vector_receiver %>% .[!is.na(names(.))]

  lsec_mono_lr_network = get_active_ligand_receptor_network(expression_vector_sender,expression_vector_receiver,lr_network,expression_cutoff_sender = 0, expression_cutoff_receiver = 4)
  mono_sig_network = get_active_signaling_network(expression_vector_receiver,sig_network,expression_cutoff_receiver = 4)
  mono_gr_network = get_active_regulatory_network(expression_vector_receiver,gr_network,expression_cutoff_receiver = 4)

  weighted_networks_monocyte = construct_weighted_networks(lsec_mono_lr_network, mono_sig_network, mono_gr_network, source_weights_df)

  expect_type(lsec_mono_lr_network,"list")
  expect_gt(nrow(lsec_mono_lr_network),0)

  expect_type(mono_sig_network,"list")
  expect_gt(nrow(mono_sig_network),0)

  expect_type(mono_gr_network,"list")
  expect_gt(nrow(mono_gr_network),0)

  expect_type(weighted_networks_monocyte,"list")
})
test_that("active ligand-target matrix and network can be constructed", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  setting = lapply(expression_settings_validation[1:2],convert_expression_settings_evaluation)
  ligands = extract_ligands_from_settings(setting)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_target_matrix_discrete = make_discrete_ligand_target_matrix(ligand_target_matrix)
  active_lt = get_active_ligand_target_matrix(setting[[1]] %>% .$response, ligand_target_matrix_discrete)
  active_lt_df = get_active_ligand_target_df(setting[[1]] %>% .$response, ligand_target_matrix_discrete)

  expect_type(active_lt_df,"list")
  expect_gt(nrow(active_lt_df),0)

  expect_type(active_lt,"logical")
  active_lt = get_active_ligand_target_matrix(setting[[1]] %>% .$response, ligand_target_matrix)
  expect_type(active_lt,"double")

})

