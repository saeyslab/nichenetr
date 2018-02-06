context("Construct ligand-target matrix functions")

test_that("Construct ligand-target matrix", {
  time <- runif(1000)
  groups <- cut(time, breaks = 4, labels = F)
  expect_equal(evaluate_trajectory(time, groups), 1)
})
