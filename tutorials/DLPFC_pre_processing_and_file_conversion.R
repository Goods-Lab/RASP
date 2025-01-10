#Load necessary libraries

library(Seurat)
library(SeuratDisk)
library("spatialLIBD")

# Download the filted feature matrix, sample 151673 from http://spatial.libd.org/spatialLIBD/
# re-name the raw count matrix: filtered_feature_bc_matrix.h5


# Load the count matrix 

counts = Read10X_h5(".../filtered_feature_bc_matrix.h5")

#Initialize seurat object 
object=CreateSeuratObject(counts = seurat_obj,assay = "RNA",min.cells = 20,min.features = 20)

#Run SCTransform on the object 
object <- SCTransform(object,assay = "RNA",
                      verbose = TRUE,
                      variable.features.n = NULL,
                      variable.features.rv.th = 1.3,return.only.var.genes = FALSE)

#Convert the object to .h5ad object for use in Python! object[["RNA"]] <- as(object = object[["RNA"]],Class = "Assay")
SaveH5Seurat(object,filename = ".../SCT.h5Seurat",overwrite = T)
Convert(".../SCT.h5Seurat",overwrite = T,assay = "SCT",dest = 'h5ad')