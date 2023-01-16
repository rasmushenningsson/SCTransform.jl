library(sctransform)
library(Seurat)

filename = "../../SingleCell10x/data/500_PBMC_3p_LT_Chromium_X_50genes/filtered_feature_bc_matrix.h5"

raw.data = Read10X_h5(filename)

clip.range = c(-sqrt(ncol(raw.data)/30), sqrt(ncol(raw.data)/30)) # default in Seurat (but not in sctransform package)

transformed = vst(raw.data, theta_regularization="log_theta", n_genes=NULL, res_clip_range=clip.range)


gm = as.numeric(transformed[["genes_log_gmean_step1"]])
write.table(gm, "../data/500_PBMC_50genes/genemean.csv", row.names=FALSE, col.names=FALSE)


model_pars = as.data.frame(transformed$model_pars)
colnames(model_pars) = c("theta", "beta0", "beta1")
write.csv(model_pars, "../data/500_PBMC_50genes/model_pars.csv", quote=FALSE)

# ------------------------------

f_ind = setdiff(1:50, c(7,15,19,22,25,42,43,47))
raw.data2 = raw.data[f_ind,]
transformed2 = vst(raw.data2, theta_regularization="log_theta", n_genes=NULL, res_clip_range=clip.range)

gm2 = as.numeric(transformed2[["genes_log_gmean_step1"]])
write.table(gm2, "../data/500_PBMC_50genes/genemean2.csv", row.names=FALSE, col.names=FALSE)

message("bandwidth: ", bw.SJ(gm2))

model_pars = as.data.frame(transformed2$model_pars)
colnames(model_pars) = c("theta", "beta0", "beta1")
write.csv(model_pars, "../data/500_PBMC_50genes/model_pars2.csv", quote=FALSE)

model_pars_fit = as.data.frame(transformed2$model_pars_fit)
colnames(model_pars_fit) = c("theta", "beta0", "beta1")
write.csv(model_pars_fit, "../data/500_PBMC_50genes/model_pars_fit2.csv", quote=FALSE)


out_c_ind = seq(from=1, to=ncol(transformed2$y), by=10)
write.table(signif(transformed2$y[,out_c_ind],digits=6), "../data/500_PBMC_50genes/transformed2_every_10th_cell.csv", sep=',', quote=FALSE, col.names=FALSE, row.names=FALSE)

# Check that calling sctransform through Seurat gives the same results
# seurat_object2 = CreateSeuratObject(counts = raw.data2)
# seurat_object2 <- Seurat::SCTransform(seurat_object2, theta_regularization="log_theta", n_genes=NULL, do.center=FALSE)
# write.table(signif(seurat_object2@assays[["SCT"]]@scale.data[,out_c_ind],digits=6), "../data/500_PBMC_50genes/transformed2_every_10th_cell_b.csv", sep=',', quote=FALSE, col.names=FALSE, row.names=FALSE)
