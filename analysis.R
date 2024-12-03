Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(ggplot2)

scRNA_harmony2<-readRDS('E:\\Lab\\GLM_CKW\\results\\scRNA_harmony2.RDS')
setwd('E:\\Lab\\GLM_CKW\\results')
H1<-Read10X('E:\\Lab\\GLM_CKW\\2024_analysis_reports_zhouyunyi__zhouyunyi__2024096_10XRNA_mouse_1053_report\\6CloupeFile\\H1\\filtered_feature_bc_matrix')
V1<-Read10X('E:\\Lab\\GLM_CKW\\2024_analysis_reports_zhouyunyi__zhouyunyi__2024096_10XRNA_mouse_1053_report\\6CloupeFile\\V1\\filtered_feature_bc_matrix')



scRNA.H1 = CreateSeuratObject(H1 ,project="HTH-01-015",min.cells = 3, min.features = 200)
scRNA.V1 = CreateSeuratObject(V1 ,project="Vehicle",min.cells = 3, min.features = 200)
scRNA_harmony <- merge(x=scRNA.H1, y=scRNA.V1)


scRNA.H1[["percent.mt"]] <- PercentageFeatureSet(scRNA.H1, pattern = "^mt-")
scRNA.V1[["percent.mt"]] <- PercentageFeatureSet(scRNA.V1, pattern = "^mt-")


VlnPlot(scRNA.H1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt" ), 
        
        pt.size = 0.01, #若不需要显示点，可以设置pt.size = 0
        ncol = 3) 


VlnPlot(scRNA.V1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt" ), 
        
        pt.size = 0.01, #若不需要显示点，可以设置pt.size = 0
        ncol = 3) 

scRNA2 <- subset(scRNA1, subset = nFeature_RNA > 500 & 
                   percent.mt < 20 &   nCount_RNA > 1000)


scRNA_harmony2<-scRNA_harmony
scRNA_harmony_data<-GetAssayData(scRNA_harmony2,slot = 'counts')
mito.genes<-grep(pattern="^mt-", x=rownames(x=scRNA_harmony_data),value=TRUE)
rib.genes<-grep(pattern="^Rp", x=rownames(scRNA_harmony_data), value=TRUE)
percent.mito<-Matrix::colSums(scRNA_harmony_data[mito.genes, ])/Matrix::colSums(scRNA_harmony_data)
x=setdiff(rownames(x =scRNA_harmony_data),rib.genes)
x=setdiff(x,mito.genes)
scRNA_harmony_data=scRNA_harmony_data[x,]


scRNA_harmony<- CreateSeuratObject(counts = scRNA_harmony_data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = scRNA_harmony@var.genes, npcs = 50, verbose = FALSE)

scRNA_harmony@meta.data<-scRNA_harmony2@meta.data
colnames(scRNA_harmony@meta.data)[1]<-'group'

unique(scRNA_harmony$group)
scRNA_harmony <- scRNA_harmony %>% 
  RunHarmony("group", plot_convergence = TRUE)
ElbowPlot(scRNA_harmony,ndims=50)
scRNA_harmony <- scRNA_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  RunTSNE(reduction = "harmony", dims = 1:40)%>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()


UMAPPlot(scRNA_harmony,group.by="seurat_clusters",label=T, pt.size = 0.5)+
  UMAPPlot(scRNA_harmony,group.by="group",label=T, pt.size = 0.5)
TSNEPlot(scRNA_harmony,group.by="group",label=F,cols = cluster_colors)



UMAPPlot(scRNA_harmony,group.by="seurat_clusters",label=T, pt.size = 0.5)+
  DotPlot(scRNA_harmony, features = c('Ptprc','Il7r','Cd3e','Cd3d','Cd4','Foxp3','Cd8a','Tcf7','Lef1',#T
                                      'Ncr1',#NK
                                      'Ms4a1','Cd19','Cd79a','Cd79b',#B
                                      'C1qa','C1qb','C1qc',#Macro
                                      'Clec9a',#cDC1
                                      'Clec10a',#cDC2
                                      'S100a9','S100a8','G0s2','Csf3r',#Neutrophils
                                      'Mzb1','Derl3','Igkc','Iglc2',#Plasma cells
                                      'Tm4sf1','Sele',#Endothelial
                                      'Dcn','Col1a2','Col1a1','Pecam1',#Stromal cells
                                      'Tpsab1','Tpsb2','Cpa3','Ms4a2',#Mast cells
                                      'Ctsk'),group.by = 'seurat_clusters')+scale_color_gradientn(colors = c("blue", "yellow", "#E95C59"))+
  RotatedAxis()# DC2


saveRDS(scRNA_harmony,'E:\\Lab\\GLM_CKW\\results\\scRNA_harmony.RDS')
scRNA_harmony<-readRDS(scRNA_harmony,'E:\\Lab\\GLM_CKW\\results\\scRNA_harmony.RDS')

scRNA_harmony2<-subset(scRNA_harmony,idents = c('0','1','2','4','6','7','8','9','10','11','12',
                                                '13','14','15','16','17','18','19','20'))
scRNA_harmony3<-scRNA_harmony2



scRNA_harmony2<- CreateSeuratObject(counts = scRNA_harmony_data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = scRNA_harmony2@var.genes, npcs = 50, verbose = FALSE)

scRNA_harmony2@meta.data<-scRNA_harmony3@meta.data
colnames(scRNA_harmony2@meta.data)[1]<-'group'

unique(scRNA_harmony2$group)
scRNA_harmony2 <- scRNA_harmony2 %>% 
  RunHarmony("group", plot_convergence = TRUE)
ElbowPlot(scRNA_harmony2,ndims=50)
scRNA_harmony2 <- scRNA_harmony2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  RunTSNE(reduction = "harmony", dims = 1:40)%>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()




UMAPPlot(scRNA_harmony2,group.by="celltype3",label=T, pt.size = 0.5)
DotPlot(scRNA_harmony2, features = c('Ptprc','Il7r','Cd3e','Cd3d','Cd4','Foxp3','Cd8a',#T
                                     'Ncr1',#NK
                                     'Ms4a1','Cd19','Cd79a','Cd79b',#B
                                     'Siglech',#pDC
                                     'Ly6c2','Nr4a1','Itgal',#Monocyto
                                     'C1qa','C1qb','C1qc',#Macro
                                     'Clec10a',#cDC2
                                     'Clec9a','Ccl22',#cDC1
                                     'S100a9','S100a8','G0s2','Csf3r'#Neutrophils
),group.by = 'celltype3')+scale_color_gradientn(colors = c('#5353ec', '#E4C755', '#E95C59'))+
  RotatedAxis()

UMAPPlot(scRNA_harmony2,group.by="seurat_clusters",label=T, pt.size = 0.5)+
  DotPlot(scRNA_harmony2, features = c('Ptprc','Il7r','Cd3e','Cd3d','Cd4','Foxp3','Cd8a',#T
                                       'Ncr1',#NK
                                       'Ms4a1','Cd19','Cd79a','Cd79b',#B
                                       'C1qa','C1qb','C1qc','Mafb','Maf','Ccl12','Mgl2','Vegfa',#Macro
                                       'Clec9a',#cDC1
                                       'Clec10a',#cDC2
                                       'S100a9','S100a8','G0s2','Csf3r',#Neutrophils
                                       'Mzb1','Derl3','Igkc','Iglc2',#Plasma cells
                                       'Dcn','Col1a2','Col1a1','Pecam1',#Stromal cells
                                       'Tpsab1','Tpsb2','Cpa3',#Mast cells
                                       'Csf1r','Ccr2',
                                       'Il1b','Wfdc17','Arg2','Ly6g',
                                       'Siglech'#pDC
  ),group.by = 'celltype')+scale_color_gradientn(colors = c("blue", "yellow", "#E95C59"))+
  RotatedAxis()# DC2



scRNA_harmony2$cellnames<-rownames(scRNA_harmony2@meta.data)
scRNA_harmony2$celltype3<-scRNA_harmony2$celltype
scRNA_harmony2$celltype3 <- factor(scRNA_harmony2$celltype3, levels = c("Gamma delta T cell","Conventional CD4","Tregs","CD8 T cells",
                                                                        "NK cells","B cells","pDC",'Monocyto',"Macrophages",
                                                                        "cDC2","cDC1","Neutrophils"))




DotPlot(scRNA_harmony2,features = c('Ly6c2','Nr4a1','Itgal'))+RotatedAxis()
UMAPPlot(scRNA_harmony2, group.by='seurat_clusters',label=T)+
  FeaturePlot(scRNA_harmony2,features = c('Ly6c2','Nr4a1','Itgal'))#+FeaturePlot(scRNA_harmony2,features = 'Cd8a')
VlnPlot(scRNA_harmony2,features = 'Lamp3')

scRNA_harmony2$celltype<-'cell'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='0')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='1')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='2')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='3')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='4')]<-'CD8 T cells'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='5')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='6')]<-'CD8 T cells'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='7')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='8')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='9')]<-'CD8 T cells'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='10')]<-'Neutrophils'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='11')]<-'cDC1'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='12')]<-'NK cells'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='13')]<-'Monocyto'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='14')]<-'Tregs'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='15')]<-'cDC1'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='16')]<-'cDC2'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='17')]<-'Macrophages'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='18')]<-'Neutrophils'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='19')]<-'B cells'
scRNA_harmony2$celltype[which(scRNA_harmony2$seurat_clusters=='20')]<-'Gamma delta T cell'

scRNA_harmony2$celltype[scRNA_harmony2$cellnames%in%cCD4]<-'Conventional CD4'


scRNA_harmony2$celltype[scRNA_harmony2$cellnames%in%Bcell]<-'B cells'
scRNA_harmony2$celltype[scRNA_harmony2$cellnames%in%pDC]<-'pDC'



TSNEPlot(scRNA_harmony2,group.by='seurat_clusters',label=T,pt.size=0.5)+
  TSNEPlot(scRNA_harmony2,group.by='celltype',label=T,pt.size=0.5)

UMAPPlot(scRNA_harmony2,group.by='seurat_clusters',label=T,pt.size=0.5)+
  UMAPPlot(scRNA_harmony2,group.by='celltype',label=T,pt.size=0.5)


Tcell<-subset()


table(scRNA_harmony2$celltype,scRNA_harmony2$group)
scRNA_harmony2_markers<-FindAllMarkers(scRNA_harmony2, logfc.threshold = 0.5,
                                       test.use = "roc", 
                                       return.thresh = 0.25, 
                                       min.pct = 0.3, only.pos = T)
write.csv(scRNA_harmony2_markers,'E:\\Lab\\GLM_CKW\\results\\scRNA_harmony2_markers.csv')





Tcell<-subset(scRNA_harmony2,idents = c('4','6'))
Tcell2<-Tcell

Tcell_data<-GetAssayData(Tcell,slot = 'counts')
Tcell<- CreateSeuratObject(counts = Tcell_data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = Tcell@var.genes, npcs = 50, verbose = FALSE)

Tcell@meta.data<-Tcell2@meta.data
colnames(Tcell@meta.data)[1]<-'group'

unique(Tcell$group)
Tcell <- Tcell %>% 
  RunHarmony("group", plot_convergence = TRUE)
ElbowPlot(Tcell,ndims=50)
Tcell <- Tcell %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  RunTSNE(reduction = "harmony", dims = 1:40)%>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 2) %>% 
  identity()


saveRDS(scRNA_harmony2,'E:\\Lab\\GLM_CKW\\results\\scRNA_harmony2.RDS')
saveRDS(Tcell,'E:\\Lab\\GLM_CKW\\results\\Tcell.RDS')
scRNA_harmony2<-readRDS('E:\\Lab\\GLM_CKW\\results\\scRNA_harmony2.RDS')
Tcell<-readRDS('E:\\Lab\\GLM_CKW\\results\\Tcell.RDS')
saveRDS(Plasma,'E:\\Lab\\GLM_CKW\\results\\Plasma.RDS')



UMAPPlot(Tcell ,label=T,group.by='seurat_clusters')+
  FeaturePlot(Tcell,features = 'Cd4')+FeaturePlot(Tcell,features = 'Cd8a')
Tcell$cellnames<-rownames(Tcell@meta.data)
cCD4<-Tcell$cellnames[which(Tcell$seurat_clusters=='8')]
VlnPlot(Tcell,features = 'Cd4')
DotPlot(Tcell, features = c('Cd4'),group.by = 'seurat_clusters')+scale_color_gradientn(colors = c("blue", "yellow", "#E95C59"))+
  RotatedAxis()# DC2



UMAPPlot(scRNA_harmony2,group.by="seurat_clusters",label=T, pt.size = 0.5)+
  DotPlot(scRNA_harmony2, features = c('Ly6g'),group.by = 'seurat_clusters')+scale_color_gradientn(colors = c("blue", "yellow", "#E95C59"))+
  RotatedAxis()# DC2



table(scRNA_harmony2$celltype2,scRNA_harmony2$group)







FeaturePlot(Plasma,features = c('Siglech','Ms4a1','Cd19','Cd79a','Cd79b'))
UMAPPlot(Plasma,split.by='group')

Idents(scRNA_harmony2)<-scRNA_harmony2$seurat_clusters
Plasma<-subset(scRNA_harmony2,idents = c('19'))
Plasma2<-Plasma

Plasma_data<-GetAssayData(Plasma,slot = 'counts')
Plasma<- CreateSeuratObject(counts = Plasma_data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = Plasma@var.genes, npcs = 20, verbose = FALSE)

Plasma@meta.data<-Plasma2@meta.data
colnames(Plasma@meta.data)[1]<-'group'

unique(Plasma$group)
Plasma <- Plasma %>% 
  RunHarmony("group", plot_convergence = TRUE)
ElbowPlot(Plasma,ndims=50)
Plasma <- Plasma %>% 
  
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
)


Plasma<-RunUMAP(Plasma,reduction = "harmony", dims = 1:10) #%>% 
#RunTSNE(reduction = "harmony", dims = 1:10, perplexity = 20)



UMAPPlot(Plasma,label=T)






Bcell<-Plasma$cellnames[which(Plasma$seurat_clusters=='0')]
pDC<-Plasma$cellnames[which(Plasma$seurat_clusters=='1')]






















library(reshape2)
#table(x)
pB2_df <- table(scRNA_harmony2$celltype3,scRNA_harmony2$group) %>% melt()
#orig.ident为按照样本的分组进行展示
colnames(pB2_df) <- c("Cluster","Group","Number")

#pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster)
pB2_df=na.omit( pB2_df)   #去除NA值
sample_color <- color
pB2 <- ggplot(data = pB2_df, aes(x = Group, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB2







pB3 <- ggplot(data = pB2_df, aes(x =Number, y = Cluster , fill = Group)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB3


sample_color <- c("Group1" = "#FF0000", "Group2" = "#0000FF")
pB3 <- ggplot(data = pB2_df, aes(x = Number, y = Cluster, fill = Group)) +
  geom_bar(stat = "identity", width=0.8, position="fill") +
  scale_fill_manual(values = sample_color) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Ratio") +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB3



cluster_colors<-c('#5353ec' , '#E95C59')
sample_color <- c("Control" = "#1f77b4", "Motolimod" = "#ff7f0e")
sample_color <- c("Control" = "#5353ec", "Motolimod" = "#E95C59")
pB3 <- ggplot(data = pB2_df, aes(x = Number, y = Cluster, fill = Group)) +
  geom_bar(stat = "identity", width=0.8, position="fill") +
  scale_fill_manual(values = sample_color) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Ratio") +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB3























Idents(scRNA_harmony2)<-scRNA_harmony2$group
unique(scRNA_harmony2$group)
deg=FindMarkers(scRNA_harmony8,ident.1 = c("HTH-01-015"),
                ident.2 = c("Vehicle"),      #FindMarkers是差异分析
                group.by = "group",logfc.threshold =0.5)
write.csv(deg,'E:\\Lab\\PRR\\实验数据\\Motolimod_Control_scRNA\\转录组\\deg.csv')

setwd('E:\\Lab\\GLM_CKW\\results\\DEG')
type<-as.character(unique(scRNA_harmony2$celltype3))
r.deg=data.frame()
for (i in 1:12) {
  Idents(scRNA_harmony2)="celltype3"
  deg=FindMarkers(scRNA_harmony2,ident.1 = "HTH-01-015",ident.2 = "Vehicle",      #FindMarkers是差异分析
                  group.by = "group",subset.ident =type[i]   )    #subset.ident =type[i]表示做提取
  
  write.csv(deg,file = paste0( type[i],'deg.csv') )
  deg$gene=rownames(deg)   ##为了防止行名出现重复，创建新的列记录细胞类型
  deg$celltype=type[i]
  deg$unm=i-1
  r.deg=rbind(deg,r.deg)
  
}







r.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0 , 'Up', 'Down'))
dim(r.deg)

r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
r.deg$unm %<>% as.vector(.) %>% as.numeric(.)

##自定义显示想要展示的基因名
##这里挑选log2FC为top5的基因进行展示

top_up_label <- r.deg %>% 
  subset(., threshold%in%"Up") %>% 
  group_by(unm) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  as.data.frame()

top_down_label <- r.deg %>% 
  subset(., threshold %in% "Down") %>% 
  group_by(unm) %>% 
  top_n(n = -5, wt = avg_log2FC) %>% 
  as.data.frame()
library(magrittr)

top_label <- rbind(top_up_label,top_down_label)
top_label$thr_signi %<>% 
  factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))
background_position <- r.deg %>%
  dplyr::group_by(unm) %>%
  dplyr::summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  as.data.frame()

## `summarise()` ungrouping output (override with `.groups` argument)
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4

### 准备绘制中间区域cluster彩色bar所需数据
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cluster_bar_position$unm %<>% 
  factor(., levels = c(0:max(as.vector(.))))

## 设置填充颜色
cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black")
cols_cluster <- c("0" = "#35978f",
                  "1" = "#8dd3c7",
                  "2" = "#ffffb3",
                  "3" = "#bebada",
                  "4" = "#fb8072",
                  "5" = "#80b1d3",
                  "6" = "#fdb462",
                  "7" = "#b3de69","8" = "#b3de89","9"='#F3B1A0',
                  "10"='#D6E7A3',
                  "11"='#57C3F3')
library(ggplot2)
library(ggrepel)
p= ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "#525252", alpha = 0.1) + ###添加灰色背景色
  geom_jitter(data = r.deg, aes(x =unm, y = avg_log2FC, colour = thr_signi),
              size = 0.5,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi) +
  scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                     breaks = seq(0, max(r.deg$unm), 1),
                     label = seq(0, max(r.deg$unm),1)) + #修改坐标轴显示刻度
  # 根据top_label标注基因名
  geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                             ymax = 0.4, fill = unm), color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "Cluster", y = "average log2FC") +
  theme_bw()


plot1 <- p + theme(panel.grid.minor = element_blank(), ##去除网格线
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour = 'black', size = 14, vjust = 65), #调整x轴坐标,vjust的值按照最终结果稍加调整
                   panel.border = element_blank(), ## 去掉坐标轴
                   axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                   axis.line.y = element_line(colour = "black")) #添加y轴坐标轴

plot1

unique(r.deg$celltype)









CD8_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/CD8 T cellsdeg.csv')
CD8_D<-CD8_DEG$X

Macro_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/Macrophagesdeg.csv')
Macro_D<-Macro_DEG$X

cDC1_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/cDC1deg.csv')
cDC1_D<-cDC1_DEG$X

cDC2_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/cDC2deg.csv')
cDC2_D<-cDC2_DEG$X

Neutro_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/Neutrophilsdeg.csv')
Neutro_D<-Neutro_DEG$X

cCD4_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/Conventional CD4deg.csv')
cCD4_D<-cCD4_DEG$X

Tregs_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/Tregsdeg.csv')
Tregs_D<-Tregs_DEG$X

B_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/B cellsdeg.csv')
B_D<-B_DEG$X

NK_DEG<-read.csv('E:\\Lab\\GLM_CKW\\results\\DEG/NK cellsdeg.csv')
NK_D<-NK_DEG$X


# install.packages('AnnotationHub')
# BiocManager::install('org.Mm.eg.db')
library(clusterProfiler)
library(AnnotationHub)
library(dbplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
gsym.id <- bitr(CD8_D[1:50], #基因名
                fromType = "SYMBOL", #从gene symbol
                toType = "ENTREZID", #提取ENTREZ ID
                OrgDb = "org.Mm.eg.db") #相应物种的包，人
egohr <- enrichGO(gene = gsym.id$ENTREZID,
                  #小鼠用这行
                  OrgDb = org.Mm.eg.db,
                  #人类用这行
                  #OrgDb = org.Hs.eg.db,
                  #非模式生物用这行，例如玉米
                  #OrgDb = maize.db,
                  ont = "BP", #或MF或CC
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,readable = T) 
barplot(egohr)
egohrresult<-egohr@result
# egohrresult<-as.matrix(egohrresult)
# barplot(egohrresult)
write.csv(egohrresult,'E:\\Lab\\GLM_CKW\\results\\DEG\\GO\\CD8_D.csv')
write.table(egohrresult,"E:/Lab/PRR/实验数据/article/Fig2_2/DEG/CD8_DEG.txt",sep = "\t",col.names = NA)

library(dplyr)
selected_pathways <- egohr@result %>%
  filter(grepl("apoptosis", Description, ignore.case = TRUE))
barplot(selected_pathways)











kegg <- enrichKEGG(gsym.id$ENTREZID, 
                   organism="mmu", #hsa,mmu
                   pvalueCutoff=0.05, 
                   #pAdjustMethod="BH",
                   keyType="kegg") #pvaluecutoff 是pvalue的阈值，显著富集性要<0.01
keggresult<-kegg@result








Idents(scRNA_harmony2)<-scRNA_harmony2$celltype3

Macrophage<-subset(scRNA_harmony2 , ident=c('Macrophages'))
Macrophage2<-Macrophage

Macrophage_data<-GetAssayData(Macrophage,slot = 'counts')
Macrophage<- CreateSeuratObject(counts = Macrophage_data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = Macrophage@var.genes, npcs = 50, verbose = FALSE)

Macrophage@meta.data<-Macrophage2@meta.data
colnames(Macrophage@meta.data)[1]<-'group'

unique(Macrophage$group)
Macrophage <- Macrophage %>% 
  RunHarmony("group", plot_convergence = TRUE)
ElbowPlot(Macrophage,ndims=50)
Macrophage <- Macrophage %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  RunTSNE(reduction = "harmony", dims = 1:40)%>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 1) %>% 
  identity()


UMAPPlot(Macrophage,label=T)

saveRDS(Macrophage,'E:\\Lab\\GLM_CKW\\results/Macrophages/Macrophage.RDS')
Macrophage<-readRDS('E:\\Lab\\GLM_CKW\\results/Macrophages/Macrophage.RDS')

Macrophage_markers<-FindAllMarkers(Macrophage, logfc.threshold = 0.5,
                                   test.use = "roc", 
                                   return.thresh = 0.25, 
                                   min.pct = 0.3, only.pos = T)
write.csv(Macrophage_markers,'E:\\Lab\\GLM_CKW\\results/Macrophages/Macrophage_markers.csv')



DotPlot(Macrophage,features = c('Ccl8'))+RotatedAxis()
UMAPPlot(Macrophage, group.by='seurat_clusters',label=T)+
  FeaturePlot(Macrophage,features = c('Ccl8'))
VlnPlot(Macrophage,features = 'Npy')

Macrophage$subcelltype<-'Macrophages'
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='0')]<-'Ifitm1+ Macrophages'#Ifitm1
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='1')]<-'Plac8+ Macrophages'#Plac8
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='2')]<-'Pf4+ Macrophages'#Pf4
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='3')]<-'Ccl8+ Macrophages'#Ccl8
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='4')]<-'Gpnmb+ Macrophages'#Mmp13,Gpnmb
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='5')]<-'Ccr2+ Macrophages'#Ccr2
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='6')]<-'Nos2+ Macrophages'
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='7')]<-'Birc5+ Macrophages'#Stmn1,Birc5
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='8')]<-'Pf4+ Macrophages'
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='9')]<-'Isg15+ Macrophages'#Isg15,Ifit2,Ifit3
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='10')]<-'Mcm6+ Macrophages'#Mcm6,Mcm3,Mcm7
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='11')]<-'NK'
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='12')]<-'Gdf15+ Macrophages'#Gdf15,M2 (tumor-supportive) macrophages may upregulate growth differentiation factor 15 (GDF15)
Macrophage$subcelltype[which(Macrophage$seurat_clusters=='13')]<-'Pf4+ Macrophages'


NK<-Macrophage$cellnames[which(Macrophage$seurat_clusters=='11')]

table(Macrophage$subcelltype,Macrophage$group)


UMAPPlot(Macrophage,group.by="seurat_clusters",label=T, pt.size = 0.5)+
  DotPlot(Tcell, features = c('C1qa','C1qb','C1qc','Ifitm1',
                              'Plac8',
                              'Pf4',
                              'Ccl8',
                              'Mmp13','Gpnmb',
                              'Ccr2',
                              'Nos2',
                              'Stmn1','Birc5',
                              'Isg15','Ifit2','Ifit3',
                              'Mcm6','Mcm3','Mcm7',
                              'Gdf15'),group.by = 'seurat_clusters')+scale_color_gradientn(colors = c("blue", "yellow", "#E95C59"))+
  RotatedAxis()




UMAPPlot(Macrophage,group.by="seurat_clusters",label=T, pt.size = 0.5)+
  DotPlot(Macrophage, features = c('Ptprc','Il7r','Cd3e','Cd3d','Cd4','Foxp3','Cd8a',#T
                                   'Ncr1',#NK
                                   'Ms4a1','Cd19','Cd79a','Cd79b',#B
                                   'Siglech',#pDC
                                   'Ly6c2','Nr4a1','Itgal',#Monocyto
                                   'C1qa','C1qb','C1qc',#Macro
                                   'Clec10a',#cDC2
                                   'Clec9a','Ccl22',#cDC1
                                   'S100a9','S100a8','G0s2','Csf3r'#Neutrophils
  ),group.by = 'seurat_clusters')+scale_color_gradientn(colors = c("blue", "yellow", "#E95C59"))+
  RotatedAxis()












my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4',
               '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', 
               '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', 
               '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E',
               '#68A180', '#3A6963', '#968175')
unique(Macrophage$subcelltype)
Idents(Macrophage)<-Macrophage$subcelltype
Macrophage2<-subset(Macrophage,ident=c("Ifitm1+ Macrophages" ,"Nos2+ Macrophages" ,  "Plac8+ Macrophages" ,  "Ccl8+ Macrophages" ,  "Gpnmb+ Macrophages" ,
                                       "Pf4+ Macrophages" ,   "Ccr2+ Macrophages"  , "Isg15+ Macrophages",  "Mcm6+ Macrophages" ,  "Birc5+ Macrophages" , "Gdf15+ Macrophages" ))

saveRDS(Macrophage2,'E:\\Lab\\GLM_CKW\\results\\Macrophages\\Macrophage2.RDS')
Macrophage2<-readRDS('E:\\Lab\\GLM_CKW\\results\\Macrophages\\Macrophage2.RDS')


library(reshape2)
#table(x)
pB2_df <- table(Macrophage2$subcelltype,Macrophage2$group) %>% melt()
#orig.ident为按照样本的分组进行展示
colnames(pB2_df) <- c("Cluster","Group","Number")

#pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster)
pB2_df=na.omit( pB2_df)   #去除NA值
sample_color <- my36colors[1:12] 
pB2 <- ggplot(data = pB2_df, aes(x = Group, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB2

prop.table(table(Macrophage2$group,Macrophage2$subcelltype),margin = 1)





pB3 <- ggplot(data = pB2_df, aes(x =Number, y = Cluster , fill = Group)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB3


sample_color <- c("Group1" = "#FF0000", "Group2" = "#0000FF")
pB3 <- ggplot(data = pB2_df, aes(x = Number, y = Cluster, fill = Group)) +
  geom_bar(stat = "identity", width=0.8, position="fill") +
  scale_fill_manual(values = sample_color) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Ratio") +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB3


unique(Macrophage2$group)
cluster_colors<-c('#5353ec' , '#E95C59')
sample_color <- c("Control" = "#1f77b4", "Motolimod" = "#ff7f0e")
sample_color <- c("HTH-01-015" = "#5353ec", "Vehicle" = "#E95C59")
pB3 <- ggplot(data = pB2_df, aes(x = Number, y = Cluster, fill = Group)) +
  geom_bar(stat = "identity", width=0.8, position="fill") +
  scale_fill_manual(values = sample_color) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Ratio") +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB3




TSNEPlot(scRNA_harmony2,label=T,group.by='celltype2',cols = my36colors, pt.size=0.5)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4',
               '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', 
               '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', 
               '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E',
               '#68A180', '#3A6963', '#968175')

color<-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#C1E6F3', '#E95C59','#58A4C3',
         '#AB3282', '#3A6963', '#968175')
UMAPPlot(scRNA_harmony2,label=T,group.by='celltype3',cols = color, pt.size=0.5)
UMAPPlot(scRNA_harmony2,label=T,group.by='group',cols = color, pt.size=0.5)
'#58A4C3','#F3B1A0'

'#E5D2DD', '#F3B1A0', '#57C3F3'

'#6778AE', '#E39A35', '#E95C59'
markers<-c('Ptprc','Il7r','Cd3e','Cd3d','Cd4','Foxp3','Cd8a',#T
           'Ncr1',#NK
           'Ms4a1','Cd19','Cd79a','Cd79b',#B
           'Siglech',#pDC
           'Ly6c2','Nr4a1','Itgal',#Monocyto
           'C1qa','C1qb','C1qc',#Macro
           'Clec10a',#cDC2
           'Clec9a','Ccl22',#cDC1
           'S100a9','S100a8','G0s2','Csf3r')
VlnPlot(scRNA_harmony2, features = markers,  
        stacked=T,pt.size=0,  
        cols = my36colors,#颜色  
        direction = "horizontal", #水平作图  
        x.lab = '', y.lab = '')+#横纵轴不标记任何东西  
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#




ccunique(scRNA_harmony2$celltype3)
Idents(scRNA_harmony2)<-scRNA_harmony2$celltype3
CD8<-subset(scRNA_harmony2,idents = c('CD8 T cells'))

CD82<-CD8

CD8_data<-GetAssayData(CD8,slot = 'counts')
CD8<- CreateSeuratObject(counts = CD8_data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = CD8@var.genes, npcs = 50, verbose = FALSE)

CD8@meta.data<-CD82@meta.data
colnames(CD8@meta.data)[1]<-'group'

unique(CD8$group)
CD8 <- CD8 %>% 
  RunHarmony("group", plot_convergence = TRUE)
ElbowPlot(CD8,ndims=50)
CD8 <- CD8 %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  RunTSNE(reduction = "harmony", dims = 1:40)%>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
Idents(CD8)<-CD8$group
UMAPPlot(CD8,label=T,split.by='group',cols = my36colors, pt.size=0.5)

library(viridis)
VlnPlot(CD8,features = c('Mki67','Ifng','Gzmb','Prf1'),group.by = 'group',ncol=2,pt.size = 0.5)
boxplot(CD8,features = c('Mki67','Ifng','Gzmb','Prf1','Pdcd1','Havcr2'))
FeaturePlot(CD8,features = c('Mki67'),reduction = 'tsne',cols = viridis(10), pt.size=1.5,split.by ='group')
saveRDS(CD8,'E:\\Lab\\GLM_CKW\\results\\CD8.RDS')
RidgePlot(CD8, features = c('Mki67','Ifng','Gzmb','Prf1','Pdcd1','Havcr2'), group.by = 'group')
# 提取表达矩阵
expression_data <- FetchData(CD8, vars = c('Mki67','Ifng','Gzmb','Prf1','Pdcd1','Havcr2', 'group'))

# 使用基础R的箱线图
boxplot(Mki67 ~ group, data = expression_data, main = 'Mki67 Expression by Group')



CD8<-readRDS('E:\\Lab\\GLM_CKW\\results\\CD8.RDS')
scRNA_harmony2$group2<-scRNA_harmony2$group
scRNA_harmony2$group[which(scRNA_harmony2$group=='HTH-01-015')]<-'HTH'
scRNA_harmony2$group <- factor(scRNA_harmony2$group, levels = c('Vehicle','HTH'))
VlnPlot(scRNA_harmony2, 
        features = c('Nuak1'), 
        group.by = 'celltype3', 
        ncol = 2, 
        pt.size = 0.5, 
        cols = c('#58A4C3','#F3B1A0'),split.by = 'group')
CD8$group[which(CD8$group=='HTH-01-015')]<-'HTH'
CD8$group <- factor(CD8$group, levels = c('Vehicle','HTH'))
VlnPlot(CD8, 
        features = c('Mki67', 'Ifng', 'Gzmb', 'Prf1'), 
        group.by = 'group', 
        ncol = 2, 
        pt.size = 0.5, 
        cols = c('#F3B1A0','#58A4C3'))+scale_x_discrete(limits = c('Vehicle','HTH'))

'#FF9900','#0066FF',


















Idents(scRNA_harmony2)<-scRNA_harmony2$celltype3
Mye<-subset(scRNA_harmony2,idents = c('Macrophages','cDC1','cDC2','Monocyto','Neutrophils'))
Mye_data<-GetAssayData(Mye,slot = 'counts')
Mye2<- CreateSeuratObject(counts = Mye_data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = Mye2@var.genes, npcs = 50, verbose = FALSE)

Mye2@meta.data<-Mye@meta.data


# unique(Mye2$group)
Mye2 <- Mye2 %>% 
  RunHarmony("group", plot_convergence = TRUE)
ElbowPlot(Mye2,ndims=50)
Mye2 <- Mye2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  RunTSNE(reduction = "harmony", dims = 1:40)%>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
Mye2_markers<-FindAllMarkers(Mye2, logfc.threshold = 0.5,
                             test.use = "roc", 
                             return.thresh = 0.25, 
                             min.pct = 0.3, only.pos = T)
write.csv(Mye2_markers,'E:\\Lab\\GLM_CKW\\results/Myeloid/Mye2_markers.csv')

unique(colnames(Mye2@meta.data))
UMAPPlot(Mye2,group.by='seurat_clusters',label=T,pt.size=0.5,split.by='group')


library(reshape2)
#table(x)
pB2_df <- table(Mye2$celltype3,Mye2$group) %>% melt()
#orig.ident为按照样本的分组进行展示
colnames(pB2_df) <- c("Cluster","Group","Number")

#pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster)
pB2_df=na.omit( pB2_df)   #去除NA值
sample_color <- my36colors 
pB2_df$Cluster <- as.factor(pB2_df$Cluster)

pB2 <- ggplot(data = pB2_df, aes(x = Group, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
pB2

Idents(Mye2)<-Mye2$group
unique(Mye2$group)
deg=FindMarkers(Mye2,ident.1 = c("HTH-01-015"),
                ident.2 = c("Vehicle"),      #FindMarkers是差异分析
                group.by = "group",logfc.threshold =0.25)
write.csv(deg,'E:\\Lab\\GLM_CKW\\results/Myeloid/deg.csv')
deg_gene<-rownames(deg)

library(clusterProfiler)
library(AnnotationHub)
library(dbplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
gsym.id <- bitr(deg_gene[1:50], #基因名
                fromType = "SYMBOL", #从gene symbol
                toType = "ENTREZID", #提取ENTREZ ID
                OrgDb = "org.Mm.eg.db") #相应物种的包，人
egohr <- enrichGO(gene = gsym.id$ENTREZID,
                  #小鼠用这行
                  OrgDb = org.Mm.eg.db,
                  #人类用这行
                  #OrgDb = org.Hs.eg.db,
                  #非模式生物用这行，例如玉米
                  #OrgDb = maize.db,
                  ont = "BP", #或MF或CC
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,readable = T) 
barplot(egohr)
egohrresult<-egohr@result





