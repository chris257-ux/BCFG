rm(list=ls())
gc()
getwd()
set.seed(256)
library(SingleR)
library(Seurat)
library(celldex)
dataset="EndometrialCarcinoma"
output_dir="singleR"
if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
DefaultAssay(scRNA)
scRNA<-seurat_combined
head(scRNA@meta.data)
resolution_cut<-"SCT_snn_res.0.5"
ref_label="hpca_bpe"
getwd()
ref_path="/home/chenliqun/project/processScRNAData/Data/Solid tumor/Bladder urothelial carcinoma/code_Bladder urothelial carcinoma/ref/"
if(ref_label=="hpca"){
print(paste("ref dataset",ref_label))
load(paste0(ref_path,"HumanPrimaryCellAtlas_hpca.se_human.RData"))
refd=hpca.se
}else if(ref_label=="bpe"){
print(paste("ref dataset",ref_label))
load(paste0(ref_path,"BlueprintEncode_bpe.se_human.RData"))
refd=bpe.se
}else if(ref_label=="hpca_bpe"){
print(paste("ref dataset",ref_label))
load(paste0(ref_path,"HumanPrimaryCellAtlas_hpca.se_human.RData"))
load(paste0(ref_path,"BlueprintEncode_bpe.se_human.RData"))
ref1=hpca.se
ref2=bpe.se
}else if(ref_label=="hpca_nh"){
print(paste("ref dataset",ref_label))
load(paste0(ref_path,"HumanPrimaryCellAtlas_hpca.se_human.RData"))
load(paste0(ref_path,"NovershternHematopoieticData.RData"))
ref1=hpca.se
ref2=nh.se
}else if(ref_label=="hpca_dic"){
print(paste("ref dataset",ref_label))
load(paste0(ref_path,"HumanPrimaryCellAtlas_hpca.se_human.RData"))
load(paste0(ref_path,"DatabaseImmuneCellExpressionData.RData"))
ref1=hpca.se
ref2=dic.se
}else if(ref_label=="hpca_mi"){
print(paste("ref dataset",ref_label))
load(paste0(ref_path,"HumanPrimaryCellAtlas_hpca.se_human.RData"))
load(paste0(ref_path,"MonacoImmuneData_monaco.se_human.RData"))
ref1=hpca.se
ref2=mi.se
}else{
print("please choose ref dataset")
}
getwd()
norm_count<-SCT_matrix
dim(norm_count)
norm_count[1:4,1:4]
if(!dir.exists("singler")){dir.create("singler",recursive=TRUE)}
if(!dir.exists(paste0("singler/","tsne"))){dir.create(paste0("singler/","tsne"),recursive=TRUE)}
if(!dir.exists(paste0("singler/","umap"))){dir.create(paste0("singler/","umap"),recursive=TRUE)}
if(ref_label=="hpca"|ref_label=="bpe"){
print("start cell type annotation with one dataset")
table(refd@colData@listData$label.main)
scRNA.pred.cell<-SingleR(test=norm_count,ref=refd,labels=refd$label.main)
table(scRNA.pred.cell$labels)
rownames(scRNA.pred.cell)[1:5]
scRNA<-AddMetaData(scRNA,scRNA.pred.cell$labels,col.name="SingleR_cell")
summary(is.na(scRNA.pred.cell$labels))
write.table(unique(scRNA.pred.cell$labels),file=paste0("singler/",dataset,"_anno_cell_",ref_label,".txt"),row.names=FALSE,quote=FALSE)
scRNA.pred.clu<-SingleR(test=norm_count,ref=refd,labels=refd$label.main,clusters=scRNA@meta.data[[resolution_cut]])
table(scRNA.pred.clu$labels)
clu_scrna<-data.frame(ClusterID=scRNA@meta.data[[resolution_cut]])
celltype=data.frame(ClusterID=rownames(scRNA.pred.clu),celltype=scRNA.pred.clu$labels,stringsAsFactors=F)
matched_index=match(clu_scrna$ClusterID,celltype$ClusterID)
annotated_clusters<-cbind(clu_scrna,celltype[matched_index,])
scRNA<-AddMetaData(scRNA,annotated_clusters$celltype,col.name="SingleR_cluster")
summary(is.na(scRNA@meta.data[["SingleR_cluster"]]))
write.table(unique(scRNA.pred.clu$labels),file=paste0("singler/",dataset,"_anno_clu_",ref_label,".txt"),row.names=FALSE,quote=FALSE)
saveRDS(scRNA,file=paste0("singler/",dataset,"_anno_",ref_label,".rds"))
plot1<-DimPlot(scRNA,reduction="tsne",group.by=c(resolution_cut,"SingleR_cell"),label=TRUE)
plot2<-DimPlot(scRNA,reduction="umap",group.by=c(resolution_cut,"SingleR_cell"),label=TRUE)
plot3<-DimPlot(scRNA,reduction="tsne",group.by=c(resolution_cut,"SingleR_cluster"),label=TRUE)
plot4<-DimPlot(scRNA,reduction="umap",group.by=c(resolution_cut,"SingleR_cluster"),label=TRUE)
ggsave(paste0("singler/","tsne/",dataset,"_anno_",ref_label,"_cell_tsne_pac.pdf"),plot=plot1,width=20,height=8)
ggsave(paste0("singler/","umap/",dataset,"_anno_",ref_label,"_cell_umap_pac.pdf"),plot=plot2,width=20,height=8)
ggsave(paste0("singler/","tsne/",dataset,"_anno_",ref_label,"_cluster_tsne_pac.pdf"),plot=plot3,width=20,height=8)
ggsave(paste0("singler/","umap/",dataset,"_anno_",ref_label,"_cluster_umap_pac.pdf"),plot=plot4,width=20,height=8)
}
if(ref_label!="hpca"&ref_label!="bpe"){
print("start cell type annotation with two reference datasets")
table(ref1@colData@listData$label.main)
table(ref2@colData@listData$label.main)
scRNA.pred.cell<-SingleR(test=norm_count,ref=list(refa=ref1,refb=ref2),labels=list(ref1$label.main,ref2$label.main))
table(scRNA.pred.cell$labels)
rownames(scRNA.pred.cell)[1:5]
scRNA<-AddMetaData(scRNA,scRNA.pred.cell$labels,col.name="SingleR_cell")
summary(is.na(scRNA.pred.cell$labels))
write.table(unique(scRNA.pred.cell$labels),file=paste0("singler/",dataset,"_anno_cell_",ref_label,".txt"),row.names=FALSE,quote=FALSE)
scRNA.pred.clu<-SingleR(test=norm_count,ref=list(refa=ref1,refb=ref2),labels=list(ref1$label.main,ref2$label.main),clusters=scRNA@meta.data[[resolution_cut]])
table(scRNA.pred.clu$labels)
clu_scrna<-data.frame(ClusterID=scRNA@meta.data[[resolution_cut]])
celltype=data.frame(ClusterID=rownames(scRNA.pred.clu),celltype=scRNA.pred.clu$labels,stringsAsFactors=F)
matched_index=match(clu_scrna$ClusterID,celltype$ClusterID)
annotated_clusters<-cbind(clu_scrna,celltype[matched_index,])
scRNA<-AddMetaData(scRNA,annotated_clusters$celltype,col.name="SingleR_cluster")
summary(is.na(scRNA@meta.data[["SingleR_cluster"]]))
write.table(unique(scRNA.pred.clu$labels),file=paste0("singler/",dataset,"_anno_clu_",ref_label,".txt"),row.names=FALSE,quote=FALSE)
saveRDS(scRNA,file=paste0("singler/",dataset,"_anno_",ref_label,".rds"))
plot1<-DimPlot(scRNA,reduction="tsne",group.by=c(resolution_cut,"SingleR_cell"),label=TRUE)
plot2<-DimPlot(scRNA,reduction="umap",group.by=c(resolution_cut,"SingleR_cell"),label=TRUE)
plot3<-DimPlot(scRNA,reduction="tsne",group.by=c(resolution_cut,"SingleR_cluster"),label=TRUE)
plot4<-DimPlot(scRNA,reduction="umap",group.by=c(resolution_cut,"SingleR_cluster"),label=TRUE)
ggsave(paste0("singler/","tsne/",dataset,"_anno_",ref_label,"_cell_tsne_pac.pdf"),plot=plot1,width=20,height=8)
ggsave(paste0("singler/","umap/",dataset,"_anno_",ref_label,"_cell_umap_pac.pdf"),plot=plot2,width=20,height=8)
ggsave(paste0("singler/","tsne/",dataset,"_anno_",ref_label,"_cluster_tsne_pac.pdf"),plot=plot3,width=20,height=8)
ggsave(paste0("singler/","umap/",dataset,"_anno_",ref_label,"_cluster_umap_pac.pdf"),plot=plot4,width=20,height=8)
}