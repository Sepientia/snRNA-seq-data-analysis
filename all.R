data_dir <- "./snRNA"
list.files(data_dir)
scRNAlist  <- list()
for (i in list.files(data_dir)) {
  counts <- Read10X(data.dir = paste0(data_dir,i))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=i,                                       
                                       min.cells=3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = i)
}
all <- merge(x=scRNAlist[[1]], y=scRNAlist[2:length(scRNAlist)])
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = '^Mt-')
all.final <- subset(all,subset = percent.mt<10)

all.final <- subset(all.final,subset = nFeature_RNA<4000)

all.final <- subset(all.final,subset = nFeature_RNA>500)

all.final <- subset(all.final,subset = nCount_RNA<8000)

all.final$seq <- "snRNA-seq"
VlnPlot(all.final, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "seq",ncol = 3,pt.size = 0)

VlnPlot(all.final, features = c("nFeature_RNA"),group.by = "seq",pt.size = 0)+
  scale_y_continuous( breaks = c(400,seq(1000, 4000, by = 1000)), limits = c(400, 4000)) 


all.final <- NormalizeData(all.final, normalization.method = "LogNormalize", 
                             scale.factor = 10000)
all.final <- FindVariableFeatures(all.final, selection.method = "vst", 
                                    nfeatures = 2000)
all.final <- ScaleData(all.final)
all.final <- RunPCA(all.final)
all.final <- FindNeighbors(all.final,  reduction = "pca",dims = 1:50)
all.final <- FindClusters(all.final, reduction = "pca", resolution = 0.4)


all.final <- RunUMAP(all.final,
                     umap.method = "uwot",
                     metric = "cosine",
                     n.components = 2,
                     n.neighbors = 50,
                     dims = 1:50,
                     spread = 0.3,
                     min.dist = 0.1,
                     reduction = "pca")
CellDimPlot(all.final,reduction = "umap",group.by = "RNA_snn_res.0.4",split.by = "group")

DimPlot(all.final,reduction = "umap",group.by = "RNA_snn_res.0.4",label = TRUE)


geneset <- list()
geneset[['fibroblasts']] <- c("Pdgfra", "Pcdh9", "Bmper")
geneset[['endothelial']] <- c("Pecam1", "Ccdc85a", "Btnl9","Vwf","Cd36")
geneset[['macrophages']] <- c("Fcgr1", "F13a1", "Adgre1")
geneset[['pericytes']] <- c("Pdgfrb", "Vtn", "Trpc3")
geneset[['adipocytes']] <- c("Adipoq", "Plin1", "Tshr")
geneset[['smooth muscle cells']] <- c("Acta2", "Myh11", "Cdh6")
geneset[['endocardial']] <- c("Npr3", "Tmem108", "Plvap")
geneset[['schwan']] <- c("Plp1", "Gfra3")
geneset[['tcell']] <- c("Itk","Cd3e","Nkg7")
geneset[['bcell']] <- c( "Pax5","Ly6d","Ms4a1")
geneset[["cardio"]] <- c("Mhrt", "Ryr2", "Slc8a1")


all.final <- AddModuleScore(all.final,features = geneset,name = names(geneset))
FeatureDimPlot(all.final,reduction = 'umap',features = "fibroblasts1")
FeatureDimPlot(all.final,reduction = 'umap',features = "endothelial2")
FeatureDimPlot(all.final,reduction = 'umap',features = "macrophages3")
FeatureDimPlot(all.final,reduction = 'umap',features = "pericytes4")
FeatureDimPlot(all.final,reduction = 'umap',features = "adipocytes5")
FeatureDimPlot(all.final,reduction = 'umap',features = "smooth muscle cells6")
FeatureDimPlot(all.final,reduction = 'umap',features = "epicardial7")
FeatureDimPlot(all.final,reduction = 'umap',features = "schwan8")
FeatureDimPlot(all.final,reduction = 'umap',features = "tcell9")
FeatureDimPlot(all.final,reduction = 'umap',features = "bcell10")
FeatureDimPlot(all.final,reduction = 'umap',features = "cardio11")

DimPlot(all.final,reduction = 'umap',group.by="RNA_snn_res.0.4")
cluster1 <- data.frame(all.final$seurat_clusters,all.final@reductions$umap@cell.embeddings,all.final$tcell12,
                       all.final$bcell14,all.final$cardiomyocytes7,all.final$endothelial13)
cluster1 <- cluster1 %>% mutate(cluster=0)
cluster1[which(cluster1$all.final.seurat_clusters%in%c(6,9,13)),8] <- "Cardiomyocytes"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(22)),8] <-"Adipocytes"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(8,20)),8] <- "Pericytes"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(18)),8] <- "Smooth Muscle Cells"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(1,3,12)),8] <- "Fibroblasts"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(0,2,7,4,14,21)),8] <- "Endothelial Cells"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(5,15,17,19)),8] <- "Macrophage"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(10)),8] <- "Epicardial Cells"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(16)),8] <-"T Cells"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(16)&cluster1$all.final.bcell14>cluster1$all.final.tcell12),8] <-"B Cells"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(1)&cluster1$umap_1< c(-2)),8] <-"Schwann Cells"
cluster1[which(cluster1$all.final.seurat_clusters%in%c(11)),8] <- "undefined"
all.final<- AddMetaData(object = all.final,
                        metadata = cluster1$cluster,
                        col.name = "cluster")
DimPlot(all.final,reduction = 'umap',group.by="cluster")
CellDimPlot(all.final,reduction = 'umap',group.by="cluster",raster = TRUE)

Idents(all.final) <-all.final$cluster
counts_data <- GetAssayData(all.final, layer = "counts")

for(group in unique(all.final$cluster)) {
  # 获取该组的细胞
  cells_in_group <- WhichCells(all.final, idents =  group)
  for(gene in key_gene) {
    # 获取基因表达
    gene_expr <- data.frame(counts_data[gene, cells_in_group])
    a <-data.frame(gene_expr)
    a$sample<-unlist(strsplit(rownames(a),"_"))[seq(1,(dim(a)[1]*2-1),2)]
    a[a$sample%in%c("control1","control2","control3","control4"),3] <- "Control"
    a[a$sample%in%c("treat1","treat2","treat3","treat4"),3] <- "DM"
    colnames(a)[3] <-"Group"
    p1 <- ggplot(a, aes(x=Group, y=gene_expr,color=Group)) +
      # 绘制箱线图
      geom_boxplot(aes(fill=Group),
                   alpha=0.1)+ # 设置透明度
      stat_boxplot(geom = "errorbar", linewidth = 0.8, width = 0.8) +
      # 绘制散点
      #geom_jitter( size = 1)+
      ylim(0,8)+
      # 设置颜色
      scale_color_manual(values = pal_npg('nrc')(9))+
      scale_fill_manual(values = pal_npg('nrc')(9))+
      # 设置主题
      theme_bw()+
      # 去除网格线
      theme(panel.grid = element_blank())+
      stat_compare_means(comparisons = list(c("DM","Control")),
                         ,method = "t.test",label.y = c(6),label="p.format")
    ggsave(filename = paste(gene,"_",group,".svg",sep = ""),p1)
  }
}
CellDimPlot(all.final,reduction = 'umap',group.by="group",raster = TRUE)
for (gene in key_gene) {
  FeaturePlot(all.final,reduction = "umap",features = gene,raster = TRUE,cols = my.colors)
  ggsave(filename = paste(gene,"_featureplot_umap_nosplit.svg"),height = 7,width = 7)
  FeaturePlot(all.final,reduction = "umap",features = gene,split.by = "group",raster = TRUE,cols = my.colors)
  ggsave(filename = paste(gene,"_featureplot_umap_split.svg"),height = 7,width = 14)
}

my.colors <- colorRampPalette(c( "grey", "red"))(100)
FeaturePlot(all.final,reduction = "umap",features = "Grk3",raster = TRUE,cols = my.colors)
FeaturePlot(all.final,reduction = "umap",features = "Grk3",split.by = "group",raster = TRUE,cols = my.colors)

names_geneset <- c("fibroblasts"  ,               "endothelial1"    ,    "macrophages"     ,    "pericytes"  ,         "adipocytes" ,        
                   "cardio" ,     "smooth muscle cells"    ,     "epicardial" ,        "schwan"        ,      "tcell","bcell" )
result_data <- data.frame()
for(func_name in names_geneset) {
  func_genes <- geneset[[func_name]]
  func_genes <- func_genes[func_genes %in% rownames(counts_data)]  # 确保基因存在
  
  for(group in unique(all.final$cluster)) {
    # 获取该组的细胞
    cells_in_group <- WhichCells(all.final, idents =  group)
    for(gene in func_genes) {
      # 获取基因表达
      gene_expr <- counts_data[gene, cells_in_group]
      #计算平均表达
      avg_expr <- mean(gene_expr)
      # 计算表达百分比
      pct_expr <- sum(gene_expr > 0) / length(gene_expr) * 100
      # 添加到结果
      result_data <- rbind(result_data, data.frame(
        Function = func_name,
        Gene = gene,
        Group = group,
        AvgExpr = avg_expr,
        PctExpr = pct_expr
      ))
    }
  }
}

# 标准化平均表达，每个基因中最大值为1
for(gene in unique(result_data$Gene)) {
  gene_idx <- result_data$Gene == gene
  max_expr <- max(result_data$AvgExpr[gene_idx])
  if(max_expr > 0) {
    result_data$AvgExpr[gene_idx] <- result_data$AvgExpr[gene_idx] / max_expr
  }
}
result_data <- result_data[which(result_data$Group!="undefined"),]
# 确保功能和基因有正确的顺序
all_genes <- geneset[names_geneset]
all_genes <- unlist(all_genes)
result_data$Function <- factor(result_data$Function, levels = rev(names_geneset))

result_data$Gene <- factor(result_data$Gene, levels = rev(all_genes))
result_data$Group <- factor(result_data$Group, levels = c("B Cells","T Cells","Schwann Cells","Epicardial Cells","Smooth Muscle Cells",
                                                          "Cardiomyocytes","Adipocytes","Pericytes","Macrophage","Endothelial Cells","Fibroblasts"))

ggplot(result_data, aes(y = Group, x = Gene, color = AvgExpr, size = PctExpr)) +
  geom_point() +
  scale_color_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),
    guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"),
    name = "Relative\nExpression"
  ) +
  scale_size_continuous(name = "Percentage\nExpressed", range = c(0, 6)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, color = "black",family = "serif",face = "bold"),  # Adjust x-axis text
    axis.text.y= element_text(angle = 0, hjust = 1, vjust = 0.5, color = "black",family = "serif",face = "bold"),
    panel.border = element_rect(fill = NA, color = "black", size = 1.5),  # Add thicker panel border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", size = 1.5),  # Add thicker facet border
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = "black", size = 1.5) # Add thicker plot area border
  ) +
  xlab("Genes") +  
  ylab("Cell Clusters")
ggsave("marker.svg",height = 5,width = 10)



FeatureDimPlot(all.final,reduction = 'umap',features = "fibroblasts1",raster = TRUE,raster.dpi = )
ggsave("marker_umap1.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "endothelial2",raster = TRUE)
ggsave("marker_umap2.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "macrophages3",raster = TRUE)
ggsave("marker_umap3.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "pericytes4",raster = TRUE)
ggsave("marker_umap4.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "adipocytes5")
ggsave("marker_umap5.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "cardiomyocytes7",raster = TRUE)
ggsave("marker_umap6.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "smooth muscle cells6",raster = TRUE)
ggsave("marker_umap7.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "epicardial7",raster = TRUE)
ggsave("marker_umap8.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "schwan8")

ggsave("marker_umap9.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "tcell9")
ggsave("marker_umap10.svg",height = 7,width = 7)

FeatureDimPlot(all.final,reduction = 'umap',features = "bcell10")
ggsave("marker_umap11.svg",height = 7,width = 7)
FeatureDimPlot(all.final,reduction = 'umap',features = "cardio11",raster = TRUE)
ggsave("marker_umap12.svg",height = 7,width = 7)

test_a <- as.data.frame(all.final[["umap"]]@cell.embeddings)
test_a$celltype <- all.final$cluster
FeatureDimPlot(all.final,reduction = 'umap',features = "bcell10")+
  stat_unchull(aes(x="umap_1",y="umap_2"),
               data=subset(test_a, celltype=="B Cells"),
               alpha = 0.2, 
               delta = 0.2,
               linewidth = 0.5,
               show.legend = F)



plotData <- as.data.frame(all.final[["umap"]]@cell.embeddings)
plotData$cluster <- all.final$cluster
plotData$seurat.cluster <- all.final$seurat_clusters

p1<-ggplot(plotData, aes(x = umap_1, y = umap_2, fill = seurat.cluster, color = seurat.cluster)) +
  stat_unchull(alpha = 0.1, size = 0.5,delta = 0.5) +
  geom_point(size = 0.5) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line()
  )

p2<-ggplot(plotData, aes(x = umap_1, y = umap_2, fill = cluster, color = cluster) )+
  stat_unchull(alpha = 0.1, size = 0.5,delta = 0.5) +
  geom_point(size = 0.5) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line()
  )


