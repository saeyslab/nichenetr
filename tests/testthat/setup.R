options(timeout = 3600)

lr_network_mouse = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
lr_network_mouse = lr_network_mouse %>% distinct(from, to)

seurat_obj_test = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
seurat_obj_test = Seurat::UpdateSeuratObject(seurat_obj_test)

if (as.numeric(substr(packageVersion("Seurat"), 1, 1)) < 5){
    seurat_obj_test = alias_to_symbol_seurat(seurat_obj_test, "mouse")
} else if (grepl("^5", packageVersion("Seurat")) & grepl("^5", seurat_obj_test@version)){
    expect_warning(alias_to_symbol_seurat(seurat_obj_test, "mouse"))
}

ligand_target_matrix_mouse = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks_mouse = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

lr_network_human = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network_human = lr_network_human %>% distinct(from, to)
gr_network_human = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))
signaling_network_human = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))