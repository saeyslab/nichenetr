library(dplyr)
library(nichenetr)
library(testthat)

test_check("nichenetr",filter = "test_use_ligand_to_target")
test_check("nichenetr",filter = "test_calculate_popularity_bias")
