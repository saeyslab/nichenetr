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


