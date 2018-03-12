context("Network extraction functions for application")
test_that("active ligand-receptor network can be constructed", {
  library(Biobase)
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

