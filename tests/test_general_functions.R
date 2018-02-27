library(dplyr)
library(testthat)

test_check("nichenetr",filter = "use_ligand_to_target")
test_check("nichenetr",filter = "calculate_popularity_bias")
test_check("nichenetr",filter = "random_networks")
