context("Model randomization functions")
#
test_that("Construct randomized networks", {
  random_lr = randomize_complete_network_source_specific(lr_network)
  expect_equal(tibble::is.tibble(random_lr), TRUE)
  expect_gte(nrow(random_lr),1)

  datasource_lr = lr_network$source[1]
  lr_randomized_source = randomize_datasource_network(datasource_lr, lr_network)
  expect_equal(tibble::is.tibble(lr_randomized_source), TRUE)
  expect_gte(nrow(lr_randomized_source),1)

  random_lr = randomize_network(lr_network)
  expect_equal(tibble::is.tibble(random_lr), TRUE)
  expect_gte(nrow(random_lr),1)

})

