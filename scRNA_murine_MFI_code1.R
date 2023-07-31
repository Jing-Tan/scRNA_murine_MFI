###Scripts of data progressing for single-cell transcriptome analyse on the maternal-fetal interface###

#############################################################################################
#########///////=========Seurat version=3.2.0========\\\\\\\#############
#############################################################################################

### Detail information see online vignettes of Seurat
### https://satijalab.org/seurat/vignettes.html

library(dplyr)
library(Seurat)
library(patchwork)


### read inputdata 
Seurat_Object <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
colnames(x = Seurat_Object) <- paste('Seurat_Object', colnames(x = Seurat_Object), sep = '_')
#### Construct Seurat object
Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = "Seurat_Object", min.cells = 3, min.features = 200)
Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Seurat_Object, pattern = "^mt-")
VlnPlot(Seurat_Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(Seurat_Object, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Seurat_Object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
### QC and selecting cells for further analysis
Seurat_Object <- subset(Seurat_Object, subset =  nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 200 & percent.mt < 8)
VlnPlot(Seurat_Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Seurat_Object <- NormalizeData(Seurat_Object, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssay(Seurat_Object,assay = "RNA")
Seurat_Object <- FindVariableFeatures(Seurat_Object, selection.method = "vst", nfeatures = 2000) ## 2000 features in default

Seurat_Object$group<-"group"
all.genes <- rownames(Seurat_Object)
Seurat_Object <- ScaleData(Seurat_Object, features = all.genes)
Seurat_Object@meta.data$time <- "timepoint"
### remove douplets cells with Doublefinder
### Detail information see online vignettes https://github.com/CCCofficial/DoubletFinder
Seurat_Object <- RunPCA(Seurat_Object)
Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:10)
Seurat_Object <- JackStraw(Seurat_Object, num.replicate = 100)
Seurat_Object <- ScoreJackStraw(Seurat_Object, dims = 1:20)
JackStrawPlot(Seurat_Object, dims = 1:20, xmax = 0.1, ymax = 0.9)
ElbowPlot(Seurat_Object)
### pK Identification (no ground-truth)
sweep.res.list_Seurat_Object <- paramSweep_v3(Seurat_Object, PCs = 1:20, sct = FALSE)
head(sweep.res.list_Seurat_Object)
sweep.stats_Seurat_Object <- summarizeSweep(sweep.res.list_Seurat_Object, GT = FALSE)
bcmvn_Seurat_Object <- find.pK(sweep.stats_Seurat_Object)
### Homotypic Doublet Proportion Estimate 

Seurat_Object <- FindNeighbors(Seurat_Object,reduction ="pca",dims = 1:20 )
Seurat_Object <- FindClusters(Seurat_Object,resolution = 0.5)
annotationscon_Seurat_Object <- Seurat_Object@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_Seurat_Object) 
nExp_poi <- round(0.08*(length(Seurat_Object@active.ident))) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_Seurat_Object <- doubletFinder_v3(Seurat_Object, PCs = 1:20, pN = 0.25, pK = 0.08)
head(seu_Seurat_Object@meta.data)
table(seu_Seurat_Object$DF.classifications_0.25_0.08)

seu_Seurat_Object@meta.data$cellfilter <- seu_Seurat_Object@meta.data$DF.classifications_0.25_0.08
seu_Seurat_Object@meta.data <-seu_Seurat_Object@meta.data[,-10]

### Extract the singlets
head(seu_Seurat_Object@meta.data)
Seurat_Object.singlet <- subset(seu_Seurat_Object, subset = cellfilter == 'Singlet')
saveRDS(Seurat_Object.singlet, file = "~/SampleName_singlet.rds")

### Perform integration
MFI.anchors <- FindIntegrationAnchors(object.list = list(SampleName_singlet1, SampleName_singlet2,~,SampleName_singlet8), dims = 1:20)
Data.integrated <- IntegrateData(anchorset = MFI.anchors)
Data.integrated <- ScaleData(Data.integrated, verbose = FALSE)
Data.integrated <- RunPCA(Data.integrated, npcs = 30, verbose = FALSE)
Data.integrated <- SCTransform(Data.integrated) %>% RunPCA()
Data.integrated <- RunHarmony(Data.integrated, group.by.vars = "time",
                              assay.use = "SCT", max.iter.harmony = 10) 
pc.num=1:10
Data.integrated <- RunUMAP(Data.integrated, reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.25) 
DimPlot(Data.integrated, reduction = "umap", label = TRUE, repel = TRUE)
saveRDS(Data.integrated, file = "~/Data.integrated.rds")

### Find Markers and assign cell type identity to clusters
new.cluster.ids <- c("Neut","MNP","B","Plas","Baso","NK1","NK2","T",
                     "Trop","Endo1","Endo2","Endo3","Peri","Mechy","Stro1","Stro2","Stro3")
names(new.cluster.ids) <- levels(Data.integrated)
Data.integrated <- RenameIdents(Data.integrated, new.cluster.ids)
Data.integrated@meta.data$celltype <- Data.integrated@active.ident
subcell.markers <- FindAllMarkers(Data.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MFI_clusters.markers <- FindAllMarkers(Data.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(Data.integrated, file = "~/Data.integrated_name.rds")

### Subset single cell type and recluster
cell_type <- subset(Data.integrated, idents = c("cell_type"))
saveRDS(Data.integrated, file = "~/cell_type.rds")
DefaultAssay(cell_type) <- "Data.integrated"
cell_type <- ScaleData(cell_type, verbose = FALSE)
cell_type <- RunPCA(cell_type, npcs = 30, verbose = FALSE)
cell_type <- SCTransform(cell_type) %>% RunPCA()
cell_type <- RunHarmony(cell_type, group.by.vars = "time",
                        assay.use = "SCT", max.iter.harmony = 10) 
pc.num=1:10
cell_type <- RunUMAP(cell_type, reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.25) 
DimPlot(cell_type, reduction = "umap", label = TRUE, repel = TRUE)
saveRDS(cell_type, file = "~/cell_type.rds")

### mapping the T and MNP subtypes to the initial cell clusters
new.cluster.ids <- c("Neut","SC_Mo","SC_M","Mo","CM_Mo","M","cDC","pDC",
                     "B","Plas","Baso","NK1","NK2","CD4_TN", "CD8_TN","MAIT_rgT","CD_TEFF","CD4-TP", "Th2","Treg",
                     "Trop1","Trop2","Endo1","Endo2","Endo3","Peri","Mechy","Stro1","Stro2","Stro3")
names(new.cluster.ids) <- levels(Data.integrated)
Data.integrated <- RenameIdents(Data.integrated, new.cluster.ids)
Data.integrated@meta.data$cellsubtype <- Data.integrated@active.ident
saveRDS(Data.integrated, file = "~/Data.integrated_name.rds")

#############################################################################################
#########///////=========SCENIC (version 1.2.4)========\\\\\\\#############
#############################################################################################
#####load the dataset 
library(SCENIC)
subcell <- subset(Data.integrated_name,idents = c("Baso"),invert = T)
subcell <- subset(subcell,downsample = 200)
setwd("~/SCENIC")
dir.create("int")

##Initialize SCENIC settings
library(SCENIC)
cellInfo <- data.frame(subcell@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="time")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "orig"
colnames(cellInfo)[which(colnames(cellInfo)=="Rename")] <- "celltype"
cellInfo <- cellInfo[,c("sample","orig","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(subcell, file="subcell.rds")

exprMat <- as.matrix(subcell@assays$RNA@counts)
mydbDIR <- "~/cisTarget_databases/cisTarget_databases"
mydbs <- c("mm9-500bp-upstream-7species.mc9nr.feather",
           "mm9-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="mgi", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "SCENIC of MFI cell")
# Save to use at a later time~
saveRDS(scenicOptions, "int/scenicOptions.rds")
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]

### The following steps are executed according to SCENIC pipeline with default parameters. 
### Please refer to the url below: 
### https://github.com/aertslab/SCENIC
### https://cloud.tencent.com/developer/article/1692240


#############################################################################################
#########///////=========seurat (version 3.2.0)Gene set scoring ========\\\\\\\#############
#############################################################################################
#####load the dataset 
MDSClist <- read.csv("~/MDSCmarker_geneset.csv",header = T,row.names = 1)
Immunecells <- subset(Data.integrated_name,idents = c("immunecells"))
gene <- MDSClist$gene
immunemdsc <- AddModuleScore(
  object = MAC2,
  features = gene,
  ctrl = 100, #Ä¬ÈÏÖµÊÇ100
  name = 'MDSC_Features')
immunemdsc@meta.data[["MDSC.Features1"]]
VlnPlot(immunemdsc,pt.size = 0.1,features = 'MDSC.Features1',cols = col)
### Detail information see online vignettes https://www.rdocumentation.org/packages/Seurat/versions/4.3.0.1/topics/AddModuleScore











#############################################################################################
#########///////=========Stereo-seq analysis ========\\\\\\\#############
#############################################################################################
#####load the dataset 

library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)
library(ggsci)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ggpubr)

setwd('~/Stereo-seq/pip')
infile <- 'SampleName.gem'
bs <- 30
pro <- tail(unlist(strsplit(infile,"/")),1)
pro <- gsub(".gem|.txt|.tsv|.gz|_filterd","",pro)
############################## 1. bin data  ##############################
dat <- fread(file = infile)
if(length(grep("MIDCounts|MIDCount",colnames(dat))>0)){
  colnames(dat) <- gsub("MIDCounts|MIDCount","UMICount",colnames(dat))}
out <- as.data.frame(dat)

dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
dat$y <- trunc((dat$y - min(dat$y)) / bs + 1)

out <- cbind(dat$y,dat$x,out)
colnames(out)[1:2] <- c(paste0("bin",bs,".y"),paste0("bin",bs,".x"))
fwrite(out,paste0(pro,"_bin",bs,"_information.txt"),col.names=T,row.names=F,sep="\t",quote=F)


dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
bin.coor <- dat[, sum(V1), by = .(x, y)]

#out <- as.data.frame(cbind(paste0("BIN.",unique(dat$bin_ID)),bin.coor$y,bin.coor$x))
#colnames(out) <- c(paste0("BIN.",bs),paste0("bin",bs,".y"),paste0("bin",bs,".x"))
#rownames(out) <- out[,1]
#fwrite(out,paste0(pro,"_bin",bs,"_position.txt"),col.names=T,row.names=F,sep="\t",quote=F)

##
geneID <- seq(length(unique(dat$geneID)))
hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
gen <- hash.G[dat$geneID, 'values']


##
bin_ID <- unique(dat$bin_ID)
hash.B <- data.frame(row.names = sprintf('%d', bin_ID), values = bin_ID)
bin <- hash.B[sprintf('%d', dat$bin_ID), 'values']


##
cnt <- dat$V1


##
rm(dat)
gc()


##
tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))

tissue_positions_list <- data.frame(row.names = paste('BIN', rownames(hash.B), sep = '.'),
                                    tissue = 1, 
                                    row = bin.coor$y, col = bin.coor$x,
                                    imagerow = bin.coor$y, imagecol = bin.coor$x)

scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))


##
mat <- sparseMatrix(i = gen, j = bin, x = cnt)

rownames(mat) <- rownames(hash.G)
colnames(mat) <- paste('BIN', sprintf('%d', seq(max(hash.B[, 'values']))), sep = '.')

seurat_spatialObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial',min.cells=5, min.features=5)


##
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) 
{
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  
  spot.radius <- unnormalized.radius / max(dim(image))
  
  return(new(Class = 'VisiumV1', 
             image = image, 
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                          fiducial = scale.factors$fiducial_diameter_fullres, 
                                          hires = scale.factors$tissue_hires_scalef, 
                                          lowres = scale.factors$tissue_lowres_scalef), 
             coordinates = tissue.positions, 
             spot.radius = spot.radius))
}

spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                  scale.factors = fromJSON(scalefactors_json), 
                                  tissue.positions = tissue_positions_list)
##
spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'

seurat_spatialObj[['slice1']] <- spatialObj



##############################  3. Spatial Data QC  ##############################
### for mouse
seurat_spatialObj[["percent.mt"]] <- PercentageFeatureSet(seurat_spatialObj, pattern = "^mt-") 
Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3+1.5*(Q3-Q1))
lower <- as.numeric(Q1-1.5*(Q3-Q1))
#raw Quality
VlnPlot(seurat_spatialObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol=3,pt.size = 0) +theme(axis.text.x=element_text(angle=20,size=0),axis.title.x=element_text(angle=20,size=10))+labs(x=paste0("nGene: ",dim(seurat_spatialObj)[1],"; ","nBIN: ",dim(seurat_spatialObj)[2]))
plot(density(seurat_spatialObj$nFeature_Spatial))+abline(v=c(500,650,800),col="grey",lwd=2,lty=6)
SpatialFeaturePlot(seurat_spatialObj, features="nFeature_Spatial", stroke=0) + theme(legend.position = "right") + scale_y_reverse()
SpatialFeaturePlot(seurat_spatialObj, features="nCount_Spatial", stroke=0) + theme(legend.position = "right") + scale_y_reverse()
#QC
pdf(paste0("figures/",pro,"_bin",bs,"_QC.Feature.pdf"))
#temp <- seurat_spatialObj
seurat_spatialObj <- subset(seurat_spatialObj,subset = nFeature_Spatial > 20 &  percent.mt < 20)
VlnPlot(seurat_spatialObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol=3,pt.size = 0) +theme(axis.text.x=element_text(angle=20,size=0),axis.title.x=element_text(angle=20,size=10))+labs(x=paste0("nGene: ",dim(seurat_spatialObj)[1],"; ","nBIN: ",dim(seurat_spatialObj)[2]))
plot(density(seurat_spatialObj$nFeature_Spatial))+abline(v=c(500,650,800),col="grey",lwd=2,lty=6)
SpatialFeaturePlot(seurat_spatialObj, features="nFeature_Spatial", stroke=0) + theme(legend.position = "right") + scale_y_reverse()
SpatialFeaturePlot(seurat_spatialObj, features="nCount_Spatial", stroke=0) + theme(legend.position = "right") 
#seurat_spatialObj <- temp
#saveRDS(seurat_spatialObj,file=paste0(pro,"_bin",bs,"_QC.rds"))

##############################  4. Data SCTransform  Clusters ##############################

seurat_spatialObj1 <- SCTransform(seurat_spatialObj, assay = "Spatial", verbose = FALSE, variable.features.n=2000)
seurat_spatialObj1 <- RunPCA(seurat_spatialObj1,assay = "SCT",verbose = F)
seurat_spatialObj1 <- FindNeighbors(seurat_spatialObj1, reduction = "pca", dims = 1:20)
seurat_spatialObj1 <- RunUMAP(seurat_spatialObj1, reduction = "pca", dims = 1:20)
seurat_spatialObj1 <- FindClusters(seurat_spatialObj1, verbose = FALSE,resolution = 1.0)
saveRDS(seurat_spatialObj1, file = "~/SampleName/seurat_spatialObj1.rds")


anchors <- FindTransferAnchors(reference = Data.integrated_name, query = seurat_spatialObj1, normalization.method = "SCT", 
                               k.anchor = 8,features = featrelist)

predictions.assay <- TransferData(anchorset = anchors, refdata = Data.integrated_name$Rename, prediction.assay = TRUE,
                                  weight.reduction = seurat_spatialObj1[["pca"]], dims = 1:30)
seurat_spatialObj1[["predictions"]] <- predictions.assay
DefaultAssay(seurat_spatialObj1) <- "predictions"
SpatialDimPlot(seurat_spatialObj1, label = TRUE, label.size=0, pt.size.factor =1.5,stroke=0.1,crop = FALSE, label.box = F,
               cols=c('Immune2' = '#008000','Immune1' = '#FFB6C1','Trop' = '#4169E1','Endo' = '#FFFF00','Stro' = '#F5F5F5'))#4169E1







#########///////=========nichenetr version=1.1.1========\\\\\\\#############
### Detail information see online vignettes of nichenetr
### https://github.com/saeyslab/nichenetr

#### Library packages ####
library(readr)
library(dplyr)
library(tidyr)
library(nichenetr)
library(Seurat) ## R version <= 4.3.0
library(ggplot2)

#### Load raw data ####
setwd('~/RawData')
options(timeout = 10000)

lr_network_human <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network_mouse <- lr_network_human %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go")) %>%
  dplyr::rename(ligand = from, receptor = to) %>% distinct() %>% 
  mutate(ligand = convert_human_to_mouse_symbols(ligand),
         receptor = convert_human_to_mouse_symbols(receptor)) %>% drop_na()%>%
  filter(bonafide == 'TRUE')

ligand_target_matrix_human <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
colnames(ligand_target_matrix_human) = ligand_target_matrix_human %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix_human) = ligand_target_matrix_human %>% rownames() %>% convert_human_to_mouse_symbols() 
ligand_target_matrix_mouse = ligand_target_matrix_human %>% .[!is.na(rownames(ligand_target_matrix_human)), !is.na(colnames(ligand_target_matrix_human))]
dim(ligand_target_matrix_mouse)

##### Differential expression analysis of 28 cellsubsets
AllCellMap <- readRDS("~/Data.integrated_name.rds")
library(Seurat)
levels(AllCellMap)
Cluster_28 <- c("PMN-MDSC","M-MDSC","Mac-MDSC","Mo","CM-Mo","Mac","cDC","pDC","B","Plas" ,
                "NK1","NK2","CD4-TN", "CD8-TN", "CD8-TEFF","MAIT/¦Ã¦ÄT","CD4-TP","Th2","Treg",
                "Trop","Endo1", "Endo2","Endo3","Peri", "Mechy","Stro1","Stro2","Stro3")
AllCellMap_28 <- subset(AllCellMap,idents = Cluster_28)
levels(AllCellMap_28)
DimPlot(AllCellMap_28,label = T,label.size = 3)+NoLegend()
AllCellMap_28 <- AllCellMap_28


AllCellMap_28.markers <- FindAllMarkers(AllCellMap_28,logfc.threshold = 0,min.pct = 0.1)
write.csv(AllCellMap_28.markers,"~/AllCellMap_28.markers.csv")
# AllCellMap_28.markers <- read.csv("~/R/NicheNet/NewResult/AllCellMap_28.markers.csv",
#                                   header = T,row.names = 1)

## MFI signature genes
gene127  <- read.csv("~/127 common signature genes of MFI-MDSCs.csv")

intersect(AllCellMap_28.markers$gene,gene127)

### DE_sender_receiver ######
## 1.5 foldchange receptors
levels(AllCellMap_28)
Non_MDSC = c("Mo","CM-Mo","Mac","cDC","pDC","B","Plas","NK1","NK2",     
             "CD4-TN","CD8-TN","CD8-TEFF","MAIT/¦Ã¦ÄT","CD4-TP","Th2" , 
             "Treg","Trop","Endo1","Endo2","Endo3","Peri",    
             "Mechy","Stro1","Stro2","Stro3")
MDSC = c("PMN-MDSC","M-MDSC","Mac-MDSC")

PMN-MDSC.markers = FindMarkers(AllCellMap_28,ident.1= MDSC[1],ident.2 = Non_MDSC)
M-MDSC.markers = FindMarkers(AllCellMap_28,ident.1= MDSC[2],ident.2 = Non_MDSC)
Mac-MDSC.markers = FindMarkers(AllCellMap_28,ident.1= MDSC[3],ident.2 = Non_MDSC)
PMN-MDSC.markers = PMN-MDSC.markers %>% mutate(gene = row.names(PMN-MDSC.markers),cluster = 'PMN-MDSC' )
M-MDSC.markers = M-MDSC.markers %>% mutate(gene = row.names(M-MDSC.markers), cluster = 'M-MDSC')
Mac-MDSC.markers = Mac-MDSC.markers %>% mutate(gene = row.names(Mac-MDSC.markers),cluster = 'Mac-MDSC')
MDSC.dif.26.markers <- rbind(rbind(PMN-MDSC.markers,M-MDSC.markers),Mac-MDSC.markers)
write.csv(MDSC.dif.26.markers,'~/MDSC.dif.26.markers.csv')

## DE_receiver_processed_receptors
DE_receiver_processed_receptors =  MDSC.dif.26.markers %>% 
  filter(receiver %in% c("PMN-MDSC","M-MDSC","Mac-MDSC")) %>%
  filter(avg_log2FC >= 0.5849625,gene %in% lr_network_mouse$receptor) %>%
  mutate(significant = p_val_adj <= 0.05, present = pct.1 >= 0.05) %>% 
  mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% 
  mutate(diff = (pct.1/pct.2)) %>% 
  mutate(score = diff*avg_log2FC) %>% arrange(-score) %>% 
  group_by(gene, receiver) %>% 
  summarise(mean_avg_log2FC = mean(avg_log2FC), min_avg_log2FC = min(avg_log2FC), mean_significant = mean(significant), mean_present = mean(present), mean_score = mean(score), min_score = min(score)) %>% 
  arrange(-min_avg_log2FC) %>% 
  rename(receptor = gene, receptor_score = min_avg_log2FC,
         receptor_significant = mean_significant,  receptor_present = mean_present)  %>% 
  ungroup() %>% filter(receptor %in% lr_network_mouse$receptor) %>% 
  mutate(scaled_receptor_score = scale_quantile_adapted(receptor_score))
receptor = unique(DE_receiver_processed_receptors$receptor)

## DE_receiver_processed_targets
DE_receiver_processed_targets = AllCellMap_28.markers %>% 
  rename(receiver = cluster) %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= 0.1) %>% 
  mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% 
  mutate(score = diff*avg_log2FC) %>% arrange(-score) %>% group_by(gene, receiver) %>% 
  summarise(mean_avg_log2FC = mean(avg_log2FC), min_avg_log2FC = min(avg_log2FC), mean_significant = mean(significant), 
            mean_present = mean(present), mean_score = mean(score), min_score = min(score)) %>% arrange(-min_avg_log2FC) %>%
  rename(target = gene, target_score = min_avg_log2FC,target_significant = mean_significant,target_present = mean_present) %>% 
  filter(target %in% gene127) %>% select(receiver, target, target_score, mean_avg_log2FC,target_significant, target_present) %>% 
  arrange(-target_score)
mean(DE_receiver_processed_targets$mean_avg_log2FC)
setdiff(gene127,unique(DE_receiver_processed_targets$target))

## DE_sender_processed_ligands
DE_sender = AllCellMap_28.markers %>%
  filter(gene %in% lr_network_mouse$ligand) %>% 
  rename(sender = cluster) %>%
  mutate(significant = p_val_adj <= 0.05, present = pct.1 >= 0.1) %>% 
  mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% 
  mutate(diff = (pct.1/pct.2)) %>% 
  mutate(score = diff*avg_log2FC) %>% arrange(-score) %>% 
  group_by(gene, sender) %>% 
  summarise(mean_avg_log2FC = mean(avg_log2FC),min_avg_log2FC = min(avg_log2FC), 
            mean_significant = mean(significant), mean_present = mean(present), 
            mean_score = mean(score), min_score = min(score)) %>% 
  arrange(-min_avg_log2FC)
DE_sender_processed_ligands = DE_sender %>% 
  rename(ligand = gene, ligand_score = min_avg_log2FC, 
         ligand_significant = mean_significant, 
         ligand_present = mean_present) %>% ungroup() %>% 
  filter(ligand %in% lr_network_mouse$ligand) %>% 
  dplyr::mutate(scaled_ligand_score = scale_quantile_adapted(ligand_score))

### DE_sender_receiver
DE_sender_receiver = lr_network_mouse %>% 
  inner_join(DE_sender_processed_ligands, by = c("ligand")) %>%
  inner_join(DE_receiver_processed_receptors, by = c("receptor")) %>% 
  select(sender, receiver, ligand, receptor, ligand_score, 
         ligand_significant, ligand_present, receptor_score, 
         receptor_significant, receptor_present, 
         scaled_ligand_score, scaled_receptor_score) %>%
  mutate(avg_score_ligand_receptor = 0.5*(ligand_score + receptor_score)) %>% 
  arrange(-avg_score_ligand_receptor) %>% ungroup()%>% 
  mutate(scaled_avg_score_ligand_receptor = scale_quantile_adapted(avg_score_ligand_receptor))
ligand = unique(DE_sender_receiver$ligand)
receptor = unique(DE_sender_receiver$receptor)

####### Ligand Activity ########
ligand_target_matrix_mouse <- as.matrix(ligand_target_matrix_mouse)
background_expressed_genes = unique(MDSC.dif.26.markers$gene)
geneset_oi = gene127
ligand_target_matrix = ligand_target_matrix_mouse[rownames(ligand_target_matrix_mouse) %in% background_expressed_genes, ]
ligands = intersect(colnames(ligand_target_matrix),ligand)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = ligands)
ligand_activities = ligand_activities  %>% dplyr::rename(ligand = test_ligand, activity = pearson) %>% dplyr::select(-aupr, -auroc) %>% filter(!is.na(activity))
best_upstream_ligands = ligand_activities %>%  arrange(-activity) %>% pull(ligand)

Ligand_All <- AllCellMap_28.markers %>% filter(gene %in% lr_network_mouse$ligand) %>% distinct()
best_upstream_ligands = intersect(best_upstream_ligands,Ligand_All$gene)

active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% 
  bind_rows() %>% na.omit()

unique(active_ligand_target_links_df$target)
unique(active_ligand_target_links_df$ligand)
intersect(active_ligand_target_links_df$ligand,Ligand_All$gene)

LR = lr_network_mouse %>%
  filter(ligand %in% intersect(active_ligand_target_links_df$ligand,
                               Ligand_All$gene),
         receptor %in% unique(DE_receiver_processed_receptors$receptor))%>% 
  select(ligand,receptor)%>% distinct()

ligand_activities = ligand_activities %>% 
  dplyr::inner_join(active_ligand_target_links_df, by = c("ligand")) %>% 
  dplyr::mutate(activity_normalized = nichenetr::scaling_zscore(activity))
scaled_ligand_activities_targets = ligand_activities %>% ungroup() %>% 
  distinct(ligand, activity, activity_normalized) %>% 
  mutate(scaled_activity_normalized = scale_quantile_adapted(activity_normalized), 
         scaled_activity = scale_quantile_adapted(activity))
ligand_activities_targets = ligand_activities  %>% 
  mutate(ligand_target_weight = weight) %>%
  inner_join(scaled_ligand_activities_targets) %>% ungroup()
ligand_TF_matrix <- ligand_activities %>%
  filter(target %in% gene127) %>%
  select(ligand,target,weight) %>%
  pivot_wider(names_from = 'target',values_from = 'weight')
ligand_TF_matrix <- ligand_TF_matrix %>% as.data.frame()
row.names(ligand_TF_matrix) <- ligand_TF_matrix$ligand
ligand_TF_matrix <- ligand_TF_matrix %>% as.data.frame()
ligand_TF_matrix <- ligand_TF_matrix[,-1]
ligand_TF_matrix[is.na(ligand_TF_matrix)] = 0
ligand_TF_matrix <- as.matrix(ligand_TF_matrix)
write.csv(ligand_TF_matrix,'~/ligand_TF_matrix_13.csv')

library(pheatmap)
pheatmap(ligand_TF_matrix,cluster_cols = F,cluster_rows = F,
         color =colorRampPalette(c( 'blue',"white","red"))(100),
         fontsize = 5,cellwidth = 4.5,cellheight = 4.5)

#### The expression of ligands, receptors and targets ####
ligand = lr_network_mouse %>%
  filter(receptor %in% unique(DE_sender_receiver$receptor),
         ligand %in% unique(active_ligand_target_links_df$ligand)) %>% 
  pull(ligand) %>% unique()
receptor = unique(DE_sender_receiver$receptor)

p <- DotPlot(AllCellMap_28,features = union(union(ligand,receptor),gene127) ,
             cols = c("#00BFFF", "#FF00FF"))
exprs_tbl <- p$data %>% as_tibble() %>% 
  rename(celltype = id,gene = features.plot,expression = avg.exp,
         expression_scaled =avg.exp.scaled,fraction = pct.exp) %>% 
  mutate(fraction = fraction/100) %>% as_tibble() %>% 
  select(celltype,gene,expression,expression_scaled,fraction) %>%
  distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene)) %>%
  filter(!is.na(fraction))   

exprs_tbl_receptor <- exprs_tbl %>% 
  filter(gene %in% receptor, celltype %in% c("PMN-MDSC","M-MDSC","Mac-MDSC") ) %>% 
  rename(receiver = celltype, receptor = gene,  receptor_expression = expression, 
         receptor_expression_scaled =expression_scaled,receptor_fraction = fraction) %>%
  mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>%
  mutate(receptor_fraction_adapted = receptor_fraction) %>%
  mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
unique(exprs_tbl_receptor$receptor)
exprs_tbl_target <- exprs_tbl %>% 
  filter(gene %in% unique(ligand_activities_targets$target),
         celltype %in% c("PMN-MDSC","M-MDSC","Mac-MDSC")) %>% 
  rename(receiver = celltype, target = gene, 
         target_expression = expression, 
         target_expression_scaled =expression_scaled,
         target_fraction = fraction ) 

exprs_tbl_ligand <- exprs_tbl %>% filter(gene %in% ligand) %>% 
  rename(sender = celltype,ligand = gene, ligand_expression = expression, 
         ligand_expression_scaled = expression_scaled,ligand_fraction = fraction ) %>%
  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>%
  mutate(ligand_fraction_adapted = ligand_fraction) %>%
  mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
unique(exprs_tbl_ligand$ligand)
exprs_sender_receiver <- lr_network_mouse %>% 
  inner_join(exprs_tbl_ligand, by=c('ligand')) %>% 
  inner_join(exprs_tbl_receptor,by=c('receptor')) %>%
  inner_join(DE_sender_receiver)
ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>% group_by(ligand,receiver) %>%
  mutate(rank_receptor_expression = dense_rank(receptor_expression),rank_receptor_fraction = dense_rank(receptor_fraction)) %>%
  mutate(ligand_scaled_receptor_expression_fraction = 0.5*((rank_receptor_fraction/max(rank_receptor_fraction))+((rank_receptor_expression/max(rank_receptor_expression))))) %>%
  distinct(ligand,receptor,receiver,ligand_scaled_receptor_expression_fraction) %>% distinct() %>% ungroup()

combined_information = DE_sender_receiver %>%
  inner_join(ligand_scaled_receptor_expression_fraction_df, by = c("receiver", "ligand", "receptor")) %>%
  inner_join(ligand_activities_targets) %>%
  left_join(DE_receiver_processed_targets) %>% ## if ligand has no target genes --> it gets NA as target value --> these ligands should not be removed --> therefore left_join instead of inner_join
  inner_join(exprs_tbl_ligand) %>%
  inner_join(exprs_tbl_receptor)  %>%
  left_join(exprs_tbl_target) %>% 
  mutate(ligand_receptor = paste(ligand, receptor, sep = "--"))
length(unique(combined_information$ligand))
combined_information = combined_information %>% 
  select(receiver, sender, ligand_receptor, ligand, receptor, target, ligand_score, ligand_significant, 
         ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction,receptor_score, 
         receptor_significant, receptor_present, receptor_expression,receptor_expression_scaled, receptor_fraction,
         ligand_scaled_receptor_expression_fraction,avg_score_ligand_receptor, 
         target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, 
         ligand_target_weight,activity, activity_normalized,scaled_ligand_score, scaled_ligand_expression_scaled, 
         scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
         scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,scaled_activity,scaled_activity_normalized) %>% distinct()
prioritizing_weights <- c("scaled_ligand_score"=5,"scaled_ligand_expression_scaled"=1,"ligand_fraction"=1,
                          "scaled_receptor_score"=1,"scaled_receptor_expression_scaled"=0.5,"receptor_fraction"=1,
                          "scaled_activity"=0,"scaled_activity_normalized"=1,"ligand_scaled_receptor_expression_fraction"=1)
combined_information_prioritized <- combined_information %>%
  dplyr::mutate(prioritization_score =
                  ((prioritizing_weights["scaled_ligand_score"] *scaled_ligand_score) +
                     (prioritizing_weights["scaled_ligand_expression_scaled"] * scaled_ligand_expression_scaled) +
                     (prioritizing_weights["scaled_receptor_score"] * scaled_receptor_score) +
                     (prioritizing_weights["scaled_receptor_expression_scaled"] * scaled_receptor_expression_scaled) +
                     (prioritizing_weights["ligand_scaled_receptor_expression_fraction"] * ligand_scaled_receptor_expression_fraction) +
                     (prioritizing_weights["scaled_activity"] * scaled_activity) +
                     (prioritizing_weights["scaled_activity_normalized"] * scaled_activity_normalized) +
                     (prioritizing_weights["ligand_fraction"] * scaled_ligand_fraction_adapted ) +
                     (prioritizing_weights["receptor_fraction"] * scaled_receptor_fraction_adapted  )
                  )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)
write.csv(combined_information_prioritized,'~/combined_information_prioritized_16.csv')

prioritization_tbl_ligand_target = combined_information_prioritized %>% 
  select(receiver, sender, ligand_receptor, ligand, receptor, target, target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, 
         ligand_target_weight,activity, activity_normalized, scaled_activity, scaled_activity_normalized,prioritization_score) %>% distinct() %>%
  arrange(sender,receiver)

top_ligand_receptor_niche_df <- prioritization_tbl_ligand_target %>% select(sender,receiver, ligand, receptor,prioritization_score)
top_ligand_receptor_niche_df$sender <- factor(top_ligand_receptor_niche_df$sender,levels = levels(AllCellMap_28))
top_ligand_receptor_niche_df$receiver <- factor(top_ligand_receptor_niche_df$receiver,levels = c("PMN-MDSC","M-MDSC","Mac-MDSC"))
top_ligand_receptor_niche_df <- top_ligand_receptor_niche_df %>%
  arrange(sender,receiver) %>% mutate(SenderReceiver = paste(sender,receiver,sep = ' -> ')) %>%
  mutate(ligand_receptor = paste(ligand, receptor,sep='-')) %>% distinct() %>%
  group_by(SenderReceiver,ligand_receptor) %>% top_n(1,prioritization_score) %>% ungroup() %>%
  select(SenderReceiver,ligand_receptor,prioritization_score,receptor)%>%distinct()


###### Doing plot ####
level_SR <- c("PMN-MDSC -> PMN-MDSC", "PMN-MDSC -> M-MDSC",  "PMN-MDSC -> Mac-MDSC"   ,   
              "M-MDSC -> PMN-MDSC",   "M-MDSC -> M-MDSC",    "M-MDSC -> Mac-MDSC",
              "Mac-MDSC -> PMN-MDSC",  "Mac-MDSC -> M-MDSC",   "Mac-MDSC -> Mac-MDSC",
              "Mo -> PMN-MDSC",       "Mo -> M-MDSC" ,       "Mo -> Mac-MDSC" ,      
              "CM-Mo -> PMN-MDSC" ,   "CM-Mo -> M-MDSC",     "CM-Mo -> Mac-MDSC",     
              "M -> PMN-MDSC"  ,      "M -> M-MDSC"  ,       "M -> Mac-MDSC" ,     
              "cDC -> PMN-MDSC" ,     "cDC -> M-MDSC" ,      "cDC -> Mac-MDSC",
              "pDC -> PMN-MDSC" ,     "pDC -> M-MDSC" ,      "pDC -> Mac-MDSC",      
              "B -> PMN-MDSC"  ,      "B -> M-MDSC" ,        "B -> Mac-MDSC" ,
              "Plas -> PMN-MDSC" ,    "Plas -> M-MDSC",      "Plas -> Mac-MDSC",
              "NK1 -> PMN-MDSC" ,     "NK1 -> M-MDSC",       "NK1 -> Mac-MDSC" ,
              "NK2 -> PMN-MDSC" ,     "NK2 -> M-MDSC",       "NK2 -> Mac-MDSC",      
              "CD4-TN -> PMN-MDSC",   "CD4-TN -> M-MDSC",    "CD4-TN -> Mac-MDSC",
              "CD8-TN -> PMN-MDSC" ,  "CD8-TN -> M-MDSC",    "CD8-TN -> Mac-MDSC",
              "CD8-TEFF -> PMN-MDSC", "CD8-TEFF -> M-MDSC",  "CD8-TEFF -> Mac-MDSC",
              "MAIT/¦Ã¦ÄT -> PMN-MDSC", "MAIT/¦Ã¦ÄT -> M-MDSC",  "MAIT/¦Ã¦ÄT -> Mac-MDSC", 
              "CD4-TP -> PMN-MDSC",   "CD4-TP -> M-MDSC",    "CD4-TP -> Mac-MDSC",
              "Th2 -> PMN-MDSC" ,     "Th2 -> M-MDSC",       "Th2 -> Mac-MDSC",
              "Treg -> PMN-MDSC",     "Treg -> M-MDSC",      "Treg -> Mac-MDSC",
              "Trop -> PMN-MDSC",     "Trop -> M-MDSC",      "Trop -> Mac-MDSC",    
              "Endo1 -> PMN-MDSC",    "Endo1 -> M-MDSC",     "Endo1 -> Mac-MDSC",    
              "Endo2 -> PMN-MDSC",    "Endo2 -> M-MDSC",     "Endo2 -> Mac-MDSC",    
              "Endo3 -> PMN-MDSC" ,   "Endo3 -> M-MDSC",     "Endo3 -> Mac-MDSC",     
              "Peri -> PMN-MDSC"  ,   "Peri -> M-MDSC",      "Peri -> Mac-MDSC",     
              "Mechy -> PMN-MDSC" ,   "Mechy -> M-MDSC",     "Mechy -> Mac-MDSC",     
              "Stro1 -> PMN-MDSC" ,   "Stro1 -> M-MDSC",     "Stro1 -> Mac-MDSC",     
              "Stro2 -> PMN-MDSC" ,   "Stro2 -> M-MDSC",     "Stro2 -> Mac-MDSC",     
              "Stro3 -> PMN-MDSC" ,   "Stro3 -> M-MDSC",     "Stro3 -> Mac-MDSC")

setdiff(unique(top_ligand_receptor_niche_df$SenderReceiver),level_SR)
top_ligand_receptor_niche_df$SenderReceiver <- factor(top_ligand_receptor_niche_df$SenderReceiver,levels = level_SR)

top_ligand_receptor_niche_df_mat <- top_ligand_receptor_niche_df %>% 
  pivot_wider(names_from = 'SenderReceiver',values_from = 'prioritization_score') %>% as.data.frame()
row.names(top_ligand_receptor_niche_df_mat) <- top_ligand_receptor_niche_df_mat$ligand_receptor
top_ligand_receptor_niche_df_mat <- as.data.frame(top_ligand_receptor_niche_df_mat)

top_ligand_receptor_niche_df_mat.1 <- top_ligand_receptor_niche_df_mat[,-c(1,2)]
top_ligand_receptor_niche_df_mat.1[is.na(top_ligand_receptor_niche_df_mat.1)] = 0
top_ligand_receptor_niche_df_mat.1 <- as.matrix(top_ligand_receptor_niche_df_mat.1)
write.csv(top_ligand_receptor_niche_df_mat.1, '~/top_ligand_receptor_niche_df_mat_127.csv')

#### Choose needed ligand-receptors acording to biological knowledge 
library(pheatmap)
## All
pheatmap(top_ligand_receptor_niche_df_mat.1, cluster_cols = F,
         cluster_rows = F,color =colorRampPalette(c("#4682B4", "white","red"))(100),
         fontsize = 5,  cellwidth = 4.5,   cellheight = 4.5)
## common ligand-receptors
common_LR_127.gene <- c(
  'Ccl7-Ccr1','Ccl12-Ccr1','Ccl8-Ccr1','Ccl5-Ccr1','Ccl2-Ccr1',
  'Ccl4-Ccr1','Ccl3-Ccr1','Cxcl12-Cxcr4','Cxcl1-Cxcr2',
  'Cxcl2-Cxcr2','Podxl2-Sell','Cd34-Sell','Cfh-Sell',
  'Selp-Selplg','Icam1-Itgal','Icam4-Itgal','F11r-Itgal',
  'Icam2-Itgal','Icam4-Itgam','Jam3-Itgam','Icam2-Itgam',
  'Vcam1-Itgb2','Icam4-Itgb2','Icam2-Itgb2','Jam3-Itgb2',
  'C3-Itgb2','Icam4-Itgb2','Icam1-Itgb2','Apoe-Sorl1',
  'Plau-Plaur','Spp1-Cd44','Fn1-Cd44','Col1a1-Cd44',
  'Lama1-Cd44','Has2-Cd44','Csf2-Csf2ra','Csf2-Csf2rb',
  'Tnf-Tnfrsf1b','Lta-Tnfrsf1b','Il17a-Il17ra',
  'Il1b-Il1r2','Il1a-Il1r2','Rps19-C5ar1','App-Fpr2',
  'Hebp1-Fpr2','Anxa1-Fpr2')#Cd80-Cd274
intersect = intersect(common_LR_127.gene,row.names(top_ligand_receptor_niche_df_mat.1))
top_ligand_receptor_niche_df_mat.com <- top_ligand_receptor_niche_df_mat.1[intersect,]
pheatmap(top_ligand_receptor_niche_df_mat.com, cluster_cols = F,cluster_rows = F,border=F,
         color =colorRampPalette(c("#4682B4", "white","red"))(100),fontsize = 6,cellwidth = 5.5,cellheight =5.5)

## Mac-MDSC specific ligand-receptors
Mac_MDSC_LR_127.gene <- c('Vtn-Itgb5','Adam9-Itgb5','Spp1-Itgb5','Sema4d-Plxnb2',
                          'C3-C3ar1', 'Tgfb1-Tgfbr1','Tgfb2-Tgfbr1','Tgfb3-Tgfbr1',
                          'Gas6-Axl','Pros1-Axl','Spn-Siglec1','A2m-Lrp1','Plat-Lrp1',
                          'Serping1-Lrp1','Serpine1-Lrp1','Serpine2-Lrp1','Plau-Lrp1',
                          'C3-Lrp1','Lrpap1-Lrp1')
intersect = intersect(Mac-MDSC_LR_127.gene,row.names(top_ligand_receptor_niche_df_mat.1))
top_ligand_receptor_niche_df_mat.mac.spe <- top_ligand_receptor_niche_df_mat.1[intersect,]
pheatmap(top_ligand_receptor_niche_df_mat.com, cluster_cols = F,cluster_rows = F,border=F,
         color =colorRampPalette(c("#4682B4", "white","red"))(100),fontsize = 6,cellwidth = 5.5,cellheight =5.5)

## Mac-MDSC and M-MDSC shared ligand-receptors
Mac_M_MDSC_LR_127.gene <- c('Ccl7-Ccr2','Ccl8-Ccr2','Ccl12-Ccr2','Ccl24-Ccr2','Ccl2-Ccr2',
                            'Ccl2-Ccr5','Ccl3-Ccr5','Ccl4-Ccr5','Ccl5-Ccr5','Ccl7-Ccr5','Ccl8-Ccr5',
                            'Ccl12-Ccr5','Cx3cl1-Cx3cr1','Col5a3-Sdc3','Csf1-Csf1r','Ifng-Ifngr2')
intersect = intersect(Mac-M_MDSC_LR_127.gene,row.names(top_ligand_receptor_niche_df_mat.1))
top_ligand_receptor_niche_df_mat.mac.m.com <- top_ligand_receptor_niche_df_mat.1[intersect,]
pheatmap(top_ligand_receptor_niche_df_mat.mac.m.com, cluster_cols = F,cluster_rows = F,border=F,
         color =colorRampPalette(c("#4682B4", "white","red"))(100),fontsize = 6,cellwidth = 5.5,cellheight =5.5)





