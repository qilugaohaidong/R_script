umi = read.table("GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",)

library(GEOquery)

gse = getGEO("GSE69405")
exprSet <- exprs(gse[[1]])
gse1 = pData(gse[[1]])

library(data.table)
system.time(gse131907 <- fread(input = "GSE131907_Lung_Cancer_raw_UMI_matrix.txt", stringsAsFactors = T, encoding = "UTF-8"))

con <- file("GSE131907_Lung_Cancer_raw_UMI_matrix.txt", "r")
line=readLines(con,n=1)
while( length(line) != 0 ) {
  print(line)
  line=readLines(con,n=1)
}
close(con)

gse69405 = read.table("GSE69405_PROCESSED_GENE_TPM_ALL.txt",header = T,row.names = 1)
gse69405.1 = aggregate(gse69405[,3:203],
                       by = list(gse69405$gene_name),
                       FUN = median)
row.names(gse69405.1) = gse69405.1[,1]
gse69405.1 = gse69405.1[,-1]


GSE131907_Lung_Cancer_raw_UMI_matrix[1:3,1:3]

aa = names(GSE131907_Lung_Cancer_raw_UMI_matrix)
gse131907 = GSE131907_Lung_Cancer_raw_UMI_matrix[,c(grep("LUNG_T",aa),grep("EBUS_06",aa),grep("EBUS_28",aa),grep("EBUS_49",aa),
                                                    grep("BRONCHO_58",aa))]


ab = names(gse131907)

ab = gsub("_LUNG_T",".LUNG_T",ab)
ab = gsub("_EBUS_",".EBUS_",ab)
ab = gsub("_BRONCHO_",".BRONCHO_",ab)
ab = data.frame(ab)


library(tidyr)

library(dplyr)

abb <- separate(ab, 1,into = c("code","id"), "[.]")

abb = unite(abb, "vs_am", `id`, `code`,sep = "_", remove = FALSE)
names(gse131907) = abb[,1]

gse2 = 


grep()
cc = gsub("\\.","_",names(gse131928))
names(gse131928) = cc

index1=grep('Tumor',gse1$description)
gse1=gse1[index1,]
gse2 = gse1

library(tidyr)
gse2 = unite(gse2, "vs_am", `patient id:ch1`, `well:ch1`, sep = "__", remove = FALSE)
row.names(gse2) = gse2$description.1

gse57872.1 = gse57872[,row.names(gse2)]
names(gse57872.1) = gse2$vs_am


pbmc1 = pbmc

ll=cbind2(pbmc,pbmc1)

library(methods)


library(limma)
library(Seurat)
library(dplyr)
library(magrittr)


data=avereps(data)



pbmc <- CreateSeuratObject(counts = matrix(as.numeric(as.matrix(gse131907)),nrow=nrow(gse131907),
                                           dimnames=list(rownames(gse131907),colnames(gse131907))),
                           project = "seurat", min.cells = 3, min.features = 50, names.delim = "_")
c[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="figure1.pdf",width=18,height=6)           #???Plot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
tiff(file="figure1.tiff",width=9000,height=3000,res = 300)           #???Plot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #???
#???(file="figure2.pdf",width=20,height=10)              #???t1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
                        CombinePlots(plots = list(plot1, plot2))
                        dev.off()
                        tiff(file="figure2.tiff",width=4000,height=2000,res = 300)              #??                     plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
                        plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
                        CombinePlots(plots = list(plot1, plot2))
                        dev.off()
                        #???                     pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
                        #??�                     pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
                        #???                     top10 <- head(x = VariableFeatures(object = pbmc), 10)
                        pdf(file="04.featureVar.pdf",width=10,height=6)             #???                     plot1 <- VariableFeaturePlot(object = pbmc)
                        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
                        CombinePlots(plots = list(plot1, plot2))
                        dev.off()
                        tiff(file="04.featureVar.tiff",width=3000,height=1800,res = 300)            #??                     plot1 <- VariableFeaturePlot(object = pbmc)
                        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
                        CombinePlots(plots = list(plot1, plot2))
                        dev.off()

write.csv(pbmc@assays[["RNA"]]@meta.features[["vst.variance"]],"vst.variance.csv")
write.csv(pbmc@assays[["RNA"]]@meta.features[["vst.variance.standardized"]],"vst.variance.standard.csv")
write.csv(pbmc@assays$RNA@data@Dimnames[[1]],"dimnames.csv")

write.csv(pbmc@reductions[["pca"]]@feature.loadings,"pca.csv")

write.csv(pbmc@meta.data[["nCount_RNA"]],"ncountrna.csv")
write.csv(pbmc@meta.data$nFeature_RNA,"nfeaturerna.csv")
write.csv(pbmc@meta.data$percent.mt,"mt.csv")


                        #####CA???�                     pbmc=ScaleData(pbmc,do.scale = TRUE,
                                       do.center = FALSE,
                                       vars.to.regress = c("percent.mt"))                     #PCA               pbmc=RunPCA(object= pbmc,pc.genes=VariableFeatures(object = pbmc))     #PCA????
                
                        #????ÿ??PC??�               pdf(file="05.pcaGene.pdf",width=10,height=8)
                        VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
                        dev.off()
                        tiff(file="05.pcaGene.tiff",width=3000,height=2400,res = 300)
                        VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
                        dev.off()
                        #???ɷַ??               pdf(file="pcacluster.pdf",width=9,height=6)
                        DimPlot(object = pbmc, reduction = "pca")
                        dev.off()
                        tiff(file="pcacluster.tiff",width=4500,height=3000,res = 300)
                        DimPlot(object = pbmc, reduction = "pca")
                        dev.off()
                        #???ɷַ??               pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
                        DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
                        dev.off()
                        tiff(file="05.pcaHeatmap.tiff",width=2000,height=1600,res = 300)
                        DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
                        dev.off()
                        #ÿ??PC??               pbmc <- JackStraw(object = pbmc, num.replicate = 100)
                        pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
                        pdf(file="figure3.pdf",width=8,height=6)
                        JackStrawPlot(object = pbmc, dims = 1:20)
                        dev.off()
                        tiff(file="figure3.tiff",width=2400,height=1800,res = 300)
                        JackStrawPlot(object = pbmc, dims = 1:20)
                        dev.off()
                        pdf(file="figure3.pdf",width=8,height=6)
                        JackStrawPlot(object = pbmc, dims = 1:20, reduction = "pca", xmax = 0.1, ymax = 1)
                        dev.off()
                        tiff(file="figure3.tiff",width=2400,height=1800,res = 300)
                        JackStrawPlot(object = pbmc, dims = 1:20, reduction = "pca", xmax = 0.1, ymax = 1)
                        dev.off() 
                        
                        ElbowPlot(pbmc)
                        ###########???????               pcSelect=20
                        pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #?????ڽ�               pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #??ϸ?????               pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)                      #TSNE????               pdf(file="figure4-a.pdf",width=15,height=10)
                        TSNEPlot(object = pbmc, 
                                 #do.label = TRUE, 
                                 pt.size = 2, label = TRUE)    #TSNE???ӻ               dev.off()
                        tiff(file="figure4-a.tiff",width=3000,height=2000,res = 300)
                        TSNEPlot(object = pbmc, 
                                 #do.label = TRUE, 
                                 pt.size = 2, label = TRUE)    #TSNE???�               dev.off()
                        write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)
                        
                        
                        ##Ѱ?Ҳ???               logFCfilter=0.5
                        adjPvalFilter=0.05
                        pbmc.markers <- FindAllMarkers(object = pbmc,
                                                       min.pct = 0.25,        
                                                       only.pos = FALSE)
                        sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
                        write.table(sig.markers,file="06.markers.xls",sep="\t",row.names=F,quote=F)
                        top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                        top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
                        
                        
                        #????marke               pdf(file="figure5.pdf",width=20,height=15)
                        DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
                        dev.off()
                        tiff(file="figure5.tiff",width=6000,height=4500,res = 300)
                        DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
                        dev.off()
                        #????marke               pdf(file="06.markerViolin.pdf",width=15,height=6)
                        VlnPlot(object = pbmc, features = c("PLAUR", "MDK","IL27RA"))
                        dev.off()
                        pdf(file="06.markerViolin1.pdf",width=15,height=6)
                        VlnPlot(object = pbmc, features = c("NRCAM", "SERPINE2","TNFSF14"))
                        dev.off()
                        
                        pdf(file="06.markerViolin2.pdf",width=15,height=6)
                        VlnPlot(object = pbmc, features = c("PDGFA", "GZMB","RTN4"))
                        dev.off()
                        #????marker??clus               pdf(file="06.markerScatter.pdf",width=10,height=6)
                        FeaturePlot(object = pbmc, features = c("PLAUR", "MDK"),cols = c("green", "red"))
                        dev.off()
                        
                        pdf("cd8a.pdf",width = 10,height = 5)
                        FeaturePlot(object = pbmc, features = c("CD8A","CD8B"),cols = c("blue", "orange"))
                        dev.off()
                        tiff("cd8a.tiff",width = 2000,height = 1000,res = 300)
                        FeaturePlot(object = pbmc, features = c("CD8A","CD8B"),cols = c("blue", "orange"))
                        dev.off()
                        
                        
                        
                        FeaturePlot(object = pbmc, features = c("CD8A"),cols = c("white", "orange"))
                        FeaturePlot(object = pbmc, features = c("CD8B"),cols = c("white", "orange"))
                        
                        pbmc.clusters = data.frame(pbmc$seurat_clusters)
                        pbmc.clusters$id = row.names(pbmc.clusters)
                        
                        
                       cd8gse = gse131907[,pbmc.clusters[which(pbmc.clusters$pbmc.seurat_clusters == 0 |pbmc.clusters$pbmc.seurat_clusters == 1|pbmc.clusters$pbmc.seurat_clusters == 2),2]]
                     
                       cd8gse = data.frame(t(cd8gse))
                       cd8gse = cd8gse[which(cd8gse$CD3D > 0|cd8gse$CD3E > 0|cd8gse$CD3G > 0),]
                       
                       cd8gse = cd8gse[which(cd8gse$CD8A > 0|cd8gse$CD8B > 0),]
                       cd8gse = cd8gse[-which(cd8gse$CD4>0) ,]
                       
                       cd8gse = data.frame(t(cd8gse))
                       cd8gse = gse131907[,names(cd8gse)]
                       
                       
                    
                       
                       cd8gse = gse131907[c("CD8A","CD8B","CD3D","CD3E","CD3G","CD4","DCT","BNC2","CRABP1",
                                            "PTGDS","FILIP1L","DKK3","AHNAK2"),
                                         singler.other[which(singler.other$singler.other == "CD8 T cells"),2]]
                       singler.other$id = row.names(singler.other)
                  
                        #????marker??clus               pdf(file="06.markerBubble.pdf",width=12,height=6)
                        cluster10Marker=names(index)
                        DotPlot(object = pbmc, features = cluster10Marker)
                        dev.off()
                        
                        library(SingleR)
                        library(celldex)
                        
                        
                        counts<-pbmc@assays$RNA@counts
                        clusters<-pbmc@meta.data$seurat_clusters
                        ann=pbmc@meta.data$orig.ident
                        singler = CreateSinglerObject(counts, annot = ann, "pbmc", min.genes = 0,
                                                      species = "Human", citation = "",
                                                      ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                                      fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
                                                      reduce.file.size = T, numCores = 1)
                        singler$seurat = pbmc
                        singler$meta.data$xy = pbmc@reductions$tsne@cell.embeddings
                        clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
                        write.table(clusterAnn,file="07.clusterAnn.txt",quote=F,sep="\t",col.names=F)
                        write.table(singler$other,file="07.cellAnn.txt",quote=F,sep="\t",col.names=F)
                        write.csv(pbmc@reductions$tsne@cell.embeddings,"tsne.csv")
                        singler.other = data.frame(singler$other)
                        
                        tsne1 = read.csv("tsne.csv",header = T,row.names = 1)
                        
                        library(ggplot2)
                        
                        pdf("figure4-b.pdf",width = 15,height = 12)
                        ggplot(tsne1, aes(tSNE_1, tSNE_2,fill = ID)) +geom_point(size=4,colour="black",alpha=0.7,shape=21)
                        dev.off()
                        tiff("figure4-b.tiff",width = 4500,height = 3600,res = 300)
                        ggplot(tsne1, aes(tSNE_1, tSNE_2,fill = ID)) +geom_point(size=4,colour="black",alpha=0.7,shape=21)
                        dev.off()
                        
                        pdf("figure4-c.pdf",width = 15,height = 12)
                        ggplot(tsne1, aes(tSNE_1, tSNE_2,fill = cell)) +geom_point(size=4,colour="black",alpha=0.7,shape=21)
                        dev.off()
                        tiff("figure4-c.tiff",width = 4500,height = 3600,res = 300)
                        ggplot(tsne1, aes(tSNE_1, tSNE_2,fill = cell)) +geom_point(size=4,colour="black",alpha=0.7,shape=21)
                        dev.off()
                        pdf("figure4-d.pdf",width = 15,height = 12)
                        ggplot(tsne1, aes(tSNE_1, tSNE_2,fill = cluster)) +geom_point(size=4,colour="black",alpha=0.7,shape=21)
                        dev.off()
                        tiff("figure4-d.tiff",width = 4500,height = 3600,res = 300)
                        ggplot(tsne1, aes(tSNE_1, tSNE_2,fill = cluster)) +geom_point(size=4,colour="black",alpha=0.7,shape=21)
                        dev.off()
                        
                   

   
 pdf("figure412-c.pdf",width = 15,height = 12)
 ggplot(tsne1, aes(tSNE_1, tSNE_2,color = cell)) +geom_point()+ 
   scale_color_manual(values=c( "blue","blue","blue","orange","blue","blue","blue","blue"))+ggtitle("CD8 T cells") +
   
   theme(plot.title = element_text(vjust = -6)) + theme(legend.position='none')
 dev.off()
 tiff("figure412-c.tiff",width = 4500,height = 3600,res = 300)
 ggplot(tsne1, aes(tSNE_1, tSNE_2,color = cell)) +geom_point()+ 
   scale_color_manual(values=c( "blue","blue","blue","orange","blue","blue","blue","blue"))+ggtitle("CD8 T cells") +
   
   theme(plot.title = element_text(vjust = -6)) + theme(legend.position='none')
 dev.off()

 
 
 
 
 pbmc1 <- CreateSeuratObject(counts = matrix(as.numeric(as.matrix(cd8gse)),
                                             nrow=nrow(cd8gse),
                                             dimnames=list(rownames(cd8gse),
                                                           colnames(cd8gse))),
                             project = "seurat", min.cells = 3, min.features = 50, names.delim = "_")
 
 
 #ʹ??PercenFeatureSet??????????ent.mt"]] <- PercentageFeatureSet(object = pbmc1, pattern = "^MT-")
 pdf(file="figure111.pdf",width=12,height=6)           #????????????ct = pbmc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
 dev.off()
 tiff(file="figure111.tiff",width=3000,height=1500,res = 300)           #?????????????ct = pbmc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
 dev.off()
 pbmc1 <- subset(x = pbmc1, subset = nFeature_RNA > 50 & percent.mt < 5)    #?????ݽ??й?????ȵ???gure211.pdf",width=20,height=10)              #?????????????tureScatter(object = pbmc1, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
 plot2 <- FeatureScatter(object = pbmc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
 CombinePlots(plots = list(plot1, plot2))
 dev.off()
 tiff(file="figure211.tiff",width=4000,height=2000,res = 300)              #?????????????tureScatter(object = pbmc1, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
 plot2 <- FeatureScatter(object = pbmc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
 CombinePlots(plots = list(plot1, plot2))
 dev.off()
 #?????ݽ??б?�malizeData(object = pbmc1, normalization.method = "LogNormalize", scale.factor = 10000)
 #??ȡ??Щ??ϸ?dVariableFeatures(object = pbmc1, selection.method = "vst", nfeatures = 2000)
 #????????????�d(x = VariableFeatures(object = pbmc1), 10)
 pdf(file="04.featureVar11.pdf",width=10,height=6)             #?????????????iableFeaturePlot(object = pbmc1)
 plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
 CombinePlots(plots = list(plot1, plot2))
 dev.off()
 tiff(file="04.featureVar11.tiff",width=3000,height=1800,res = 300)            #?????????????iableFeaturePlot(object = pbmc1)
 plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
 CombinePlots(plots = list(plot1, plot2))
 dev.off()
 
 
 ###############??####pbmc1=Scata(pbmc1,do.scale = TRUE,
                do.center = FALSE,
                vars.to.regress = c("percent.mt")) 
                    #PCA??ά֮ǰ?�t= pbmc1,pc.genes=VariableFeatures(object = pbmc1))     #PCA????
 
 #????ÿ?�???ne11.pdf",width=10,height=8)
 VizDimLoadings(object = pbmc1, dims = 1:4, reduction = "pca",nfeatures = 20)
 dev.off()
 tiff(file="05.pcaGene11.tiff",width=3000,height=2400,res = 300)
 VizDimLoadings(object = pbmc1, dims = 1:4, reduction = "pca",nfeatures = 20)
 dev.off()
 #???ɷַ???ͼ??
 pdfer11.pdf",width=9,height=6)
 DimPlot(object = pbmc1, reduction = "pca")
 dev.off()
 tiff(file="pcacluster11.tiff",width=4500,height=3000,res = 300)
 DimPlot(object = pbmc1, reduction = "pca")
 dev.off()
 #???ɷַ?????ͼ
 pdfatmap11.pdf",width=10,height=8)
 DimHeatmap(object = pbmc1, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
 dev.off()
 tiff(file="05.pcaHeatmap11.tiff",width=2000,height=1600,res = 300)
 DimHeatmap(object = pbmc1, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
 dev.off()
 #ÿ??PC??pֵ?ֲ??;?(object = pbmc1, num.replicate = 100)
 pbmc1 <- ScoreJackStraw(object = pbmc1, dims = 1:20)
 pdf(file="figure311.pdf",width=8,height=6)
 JackStrawPlot(object = pbmc1, dims = 1:20)
 dev.off()
 tiff(fi
e="figure311.tiff",width=2400,height=1800,res = 300)
 JackStrawPlot(object = pbmc1, dims = 1:20)
 dev.off()
 pdf(file="figure311.pdf",width=8,height=6)
 JackStrawPlot(object = pbmc1, dims = 1:20, reduction = "pca", xmax = 0.1, ymax = 0.3)
 dev.off()
 tiff(file="figure311.tiff",width=2400,height=1800,res = 300)
 JackStrawPlot(object = pbmc1, dims = 1:20, reduction = "pca", xmax = 0.1, ymax = 0.3)
 dev.off() 
 
 #ElbowPlot(pbmc)
 
 ####################ker???? <- FindNeighbors(object = pbmc1, dims = 1:pcSelect)                #?????ers(object = pbmc1, resolution = 0.5)                  #??ϸ??????,?Ż???׼bject = pbmc1, dims = 1:pcSelect)                      #TSNE????
 pdf(file=1-a.pdf",width=12,height=8)
 TSNEPlot(object = pbmc1, 
          #do.label = TRUE, 
          pt.size = 2, label = TRUE)    #TSNE???ӻ?
 devle="figure411-a.tiff",width=2200,height=2000,res = 300)
 TSNEPlot(object = pbmc1, 
          #do.label = TRUE, 
          pt.size = 2, label = TRUE)    #TSNE???ӻ?
 dev.offable(pbmc1$seurat_clusters,file="06.tsneCluster11.txt",quote=F,sep="\t",col.names=F)
 
 write.csv(pbmc1@reductions$tsne@cell.embeddings,"tsne1.csv")
 
 
 ##Ѱ?Ҳ?????????????djPvalFilter=0.05
 pbmc.markers1 <- FindAllMarkers(object = pbmc1,
                                min.pct = 0.25,        
                                only.pos = FALSE)
 sig.markers1=pbmc.markers1[(abs(as.numeric(as.vector(pbmc.markers1$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers1$p_val_adj))<adjPvalFilter),]

 
 
 expr <- cd8gse
 expr$ID <- rownames(expr)
 s2e <- bitr(expr$ID, 
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)#??
 expr <- inner_j(expr,s2e,by=c("ID"="SYMBOL"))
 rownames(expr) <- expr$ENTREZID
 expr <- expr[,-5755]
 expr <- expr[,-5754]

 
 kegggmt <- read.gmt("c2.cp.kegg.v7.4.entrez.gmt")
 colnames(kegggmt)
 kegg_list = split(kegggmt$gene, kegggmt$term)
 
 expr=as.matrix(expr)
 kegg21 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)
 kegg21 = data.frame(kegg21)
 normalize = function(x){
   return((x-min(x))/(max(x)-min(x)))}
 #??ssGSEA score???н?
 kegg11111 
 kegg3=normalize(kegg21)
 


 kegggmt <- read.gmt("h.all.v7.4.entrez (1).gmt")
 colnames(kegggmt)
 kegg_list = split(kegggmt$gene, kegggmt$term)
 
 
 kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)
 kegg2 = data.frame(kegg2)
 normalize = function(x){
   return((x-min(x))/(max(x)-min(x)))}
 #??ssGSEA score???н?
 kegg111 =  
 
  write.table(sig.markers1,file="06.markers11.xls",sep="\t",row.names=F,quote=F)
 top10 <- pbmc.markers1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
 top5 <- sig.markers1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
 
 
 #????marker?ڸ???clus???ͼ
pdf",width=20,height=15)
 DoHeatmap(object = pbmc1, features = top10$gene) + NoLegend()
 dev.off()
 tiff(file="figure5.tiff",width=6000,height=4500,res = 300)
 DoHeatmap(object = pbmc1, features = top10$gene) + NoLegend()
 dev.off()
 #????marker??С????ͼf(filerViolin.pdf",width=15,height=6)
 VlnPlot(object = pbmc, features = c("PLAUR", "MDK","IL27RA"))
 dev.off()
 pdf(file="06.markerViolin1.pdf",width=15,height=6)
 VlnPlot(object = pbmc, features = c("NRCAM", "SERPINE2","TNFSF14"))
 dev.off()
 
 pdf(file="06.markerViolin2.pdf",width=15,height=6)
 VlnPlot(object = pbmc, features = c("CD8B"))
 dev.off()
 #????marker?ڸ???clus?ɢ??�rScatter.pdf",width=10,height=6)
 FeaturePlot(object = pbmc, features = c("PLAUR", "MDK"),cols = c("green", "red"))
 dev.off()
 FeaturePlot(object = pbmc1, features = c("GZMK","GZMA"),cols = c("blue", "orange"))
 

 
 VlnPlot(object = pbmc1, features = c("CCR7","TCF7","LEF1","SELL")) # cen Memory T 3 11
 VlnPlot(object = pbmc1, features = c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB","GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7")) # cen Memory T 3
 
 VlnPlot(object = pbmc1, features = c( 
                                  "CTLA4", "PDCD1", "LAG3", "HAVCR2", 
                                      "TIGIT"))

 
 write.csv(sig.markers1,"sig.markers1.csv")
 
 write.csv(sig.markers11,"sig.markers11.csv")
 
 #????marker?ڸ???clus?????�rBubble.pdf",width=12,height=6)
 cluster10Marker=names(index)
 DotPlot(object = pbmc, features = cluster10Marker)
 dev.off()
 
 twicegene = read.csv("gene.csv",header = F)
 
 new.cluster.ids <- c("Cytotoxic CD8 T cells1","Cytotoxic CD8 T cells2","Cytotoxic CD8 T cells3","Cytotoxic CD8 T cells4","Cytotoxic CD8 T cells5","Cytotoxic CD8 T cells6",
                      "naive/memory CD8 T cells","Exhausted CD8 T cells","Cytotoxic CD8 T cells7","Cytotoxic CD8 T cells8")
 
 new.cluster.ids <- c(0:9)
 
 
 names(new.cluster.ids) <- levels(pbmc1)
 pbmc1 <- RenameIdents(pbmc1, new.cluster.ids)
 pdf(file="figure6-b.pdf",width=12,height=8)
 TSNEPlot(object = pbmc1, 
          #do.label = TRUE, 
          pt.size = 2, label = TRUE,  label.size = 2.5)    #TSNE???ӻ?
 dev.off()
 tle="figure6-b.tiff",width=3000,height=2000,res = 300)
 TSNEPlot(object = pbmc1, 
          #do.label = TRUE, 
          pt.size = 2, label = TRUE,label.size = 2.5)    #TSNE???ӻ?
 dev.off()
 
te.csv(pbmc1@active.ident,"activa.csv")
 
 DimPlot(pbmc1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
 
 
 
 FeaturePlot(object = pbmc, features = 'PC_1')
 
                        hpca.se <- HumanPrimaryCellAtlasData() 
                        clusters1<-pbmc1@meta.data$seurat_clusters
                        write.csv(clusters1,"clusters1.csv")
                        
                        clusters2$pbmc1.seurat_clusters
                        mydata1 = data.frame(cell = c("Cytotoxic CD8 T cells1","Cytotoxic CD8 T cells2","Cytotoxic CD8 T cells3","Cytotoxic CD8 T cells4","Cytotoxic CD8 T cells5","Cytotoxic CD8 T cells6",
                                                     "naive/memory CD8 T cells","Exhausted CD8 T cells","Cytotoxic CD8 T cells7","Cytotoxic CD8 T cells8"),pbmc1.seurat_clusters= c(0,1,2,3,4,5,6,7,8,9))
                        

                      clusters3 =  merge(mydata1,clusters2)

                      
                      write.csv(clusters3[,2:3],"clusters3.csv")
 
                 table(pbmc1$seurat_clusters)/5509*100     
                      sum(table(pbmc1$seurat_clusters))
                      
                      
                      mydata = data.frame(cell = c("Cytotoxic CD8 T cells1","Cytotoxic CD8 T cells2","Cytotoxic CD8 T cells3","Cytotoxic CD8 T cells4","Cytotoxic CD8 T cells5","Cytotoxic CD8 T cells6",
                              "naive/memory CD8 T cells","Exhausted CD8 T cells","Cytotoxic CD8 T cells7","Cytotoxic CD8 T cells8"),`Fraction of cells`=                       c( table(pbmc1$seurat_clusters)/5509 )
)
         
 mydata1 = mydata
 mydata1 = data.frame()
 for (i in 1:10) {
   mydata1 = rbind(mydata1,sig.markers1[which(sig.markers1$cluster == c(i-1)),7])
 }
 
 mydata1 = data.frame(t(mydata1))
 
 for (i in 1:10) {
   mydata1[,i] [which(duplicated(mydata1[,i]) == "TRUE")] = ""
 }
 
 
 mydata1 = cbind(mydata[,1],rep("na",10),t(mydata1))
 
 write.csv(mydata1,"mydata11.csv")
 
 names(mydata)[2] = "Fraction of cells"               
               sum(mydata$`Frac of cells`)  
               library(ggplot2)
               library(RColorBrewer)
 pdf("figure7.pdf",width = 6,height = 4)                  
               ggplot(data=mydata,aes(cell,`Fraction of cells`,fill = cell))+
                           geom_bar(stat="identity", color="black",position=position_dodge(0.7), width=0.8,fill= c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                                                                                   "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF"
                          ),size=0.25) +#"#00AFBB"
                          # scale_fill_manual(values=brewer.pal(9,"YlOrRd")[c(6:2)])+
                           coord_flip()+xlab("")
                           theme(
                              axis.title=element_text(size=15,face="plain",color="black"),
                              axis.text = element_text(size=12,face="plain",color="black"),
                              legend.title=element_text(size=13,face="plain",color="white"),
                              legend.position = "right"# c(0.83,0.15)
                           )
                        
  dev.off()                      
                        
                        
  tiff("figure7.tiff",width = 1800,height = 1200,res = 300)                  
  ggplot(data=mydata,aes(cell,`Fraction of cells`,fill = cell))+
     geom_bar(stat="identity", color="black", position=position_dodge(0.5),width=0.8,fill= c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                                                             "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF"
     ),size=0.25) +#"#00AFBB"
     # scale_fill_manual(values=brewer.pal(9,"YlOrRd")[c(6:2)])+
     coord_flip()+xlab("")
  theme(
     axis.title=element_text(size=15,face="plain",color="black"),
     axis.text = element_text(size=12,face="plain",color="black"),
     legend.title=element_text(size=13,face="plain",color="white"),
     legend.position = "right"# c(0.83,0.15)
  )
  
  dev.off()                      
  
  pbmc1@assays$RNA@scale.data <- scale(pbmc1@assays$RNA@data, scale = TRUE)  
  #????marker?ڸ???clus???ͼ
.pdf",width=20,height=15)
  DoHeatmap(object = pbmc1,features = c("CCR7","TCF7","LEF1","SELL","PRF1", "IFNG", "GNLY", "NKG7", "GZMB","GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7",
                                        "CTLA4", "PDCD1", "LAG3", "HAVCR2", 
                                        "TIGIT")) 
  dev.off()
  tiff(file="figure8.tiff",width=4000,height=3000,res = 300)
  DoHeatmap(object = pbmc1,features = c("CCR7","TCF7","LEF1","SELL","PRF1", "IFNG", "GNLY", "NKG7", "GZMB","GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7",
                                        "CTLA4", "PDCD1", "LAG3", "HAVCR2", 
                                        "TIGIT"))  
  dev.off()                  
   
  
  pdf(file="figure8123.pdf",width=20,height=15)
  DoHeatmap(object = pbmc1,features = sig.markers11[,7]) 
  dev.off()
  tiff(file="figure8123.tiff",width=4000,height=3000,res = 300)
  DoHeatmap(object = pbmc1,features = sig.markers11[,7])  
  dev.off()                  
  

  
  write.csv(mydata,"mydata.csv")

  
  
  
  
  59 524
  
  
  tcga.KIRC.phenotype = read.csv("TCGA-LUAD.GDC_phenotype.csv",header = T,row.names = 1)
  tcga.KIRC.RNA = read.table(file = "TCGA-LUAD.htseq_fpkm.tsv.gz",sep = "\t",header = T,row.names = 1)
  tcga.KIRC.anno = read.table("gencode.v22.annotation.gene.probeMap",sep = "\t",header = T,row.names = 1)
  tcga.KIRC.sur = read.csv("TCGA-LUAD.survival.csv",header = T,row.names = 1)          
  
  
  tcga.KIRC.phenotype1 = tcga.KIRC.phenotype[which(tcga.KIRC.phenotype$sample_type_id.samples == "1"),]
  tcga.KIRC.phenotype1 = t(data.frame(t(tcga.KIRC.phenotype1)))
  tcga.KIRC.phenotype11 = tcga.KIRC.phenotype[which(tcga.KIRC.phenotype$sample_type_id.samples == "11"),]
  tcga.KIRC.phenotype11 = t(data.frame(t(tcga.KIRC.phenotype11)))
  
  tcga.KIRC.RNA1 = data.frame(tcga.KIRC.RNA[,intersect(row.names(tcga.KIRC.phenotype1),names(tcga.KIRC.RNA))],
                              tcga.KIRC.anno[row.names(tcga.KIRC.RNA),])
  
  tcga.KIRC.RNA11 = data.frame(tcga.KIRC.RNA[,c(intersect(row.names(tcga.KIRC.phenotype11),names(tcga.KIRC.RNA)),intersect(row.names(tcga.KIRC.phenotype1),names(tcga.KIRC.RNA)))],
                               tcga.KIRC.anno[row.names(tcga.KIRC.RNA),])
  
  
  tcga.KIRC.RNA21 <- aggregate(tcga.KIRC.RNA11[,1:583],
                               by = list(tcga.KIRC.RNA11$gene),
                               FUN = median)
  row.names(tcga.KIRC.RNA21) = tcga.KIRC.RNA21[,1]
  tcga.KIRC.RNA21 = tcga.KIRC.RNA21[,-1]
  tcga.KIRC.RNA31 = tcga.KIRC.RNA21[which(rowSums(tcga.KIRC.RNA21) > 0),]
  tcga.KIRC.RNA31 = normalizeBetweenArrays(tcga.KIRC.RNA31)
  
  
  
  
  
  inputFile= data.frame(tcga.KIRC.RNA31)                                        #?????ļ?
  inputFilFile[,intersect(row.names(tcga.KIRC.phenotype1),names(inputFile))]
  
  
  gmtFile =  "immune.gmt"                                       #GMT?ļ?
  #GMT?ļ??ð?
 Tary(GSary(limma)
  library(GSEABase)
  
  #??ȡ?????ļ????????.matrix(tcga.KIRC.RNA31[,60:583])
  

  dimnames=list(rownames(rt),colnames(rt))
  mat=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
  #mat=avereps(mat)
 # mat=mat[rowMeans(mat)>0,]
  geneSet=getGmt(gmtFile, 
                 geneIdType=SymbolIdentifier())
  
  #ssgsea????
  ssgseaScore =a(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  #????ssGSEA score????
  normalizeion(x){
     return((x-min(x))/(max(x)-min(x)))}
  #??ssGSEA score???н?
  ssgseaOutze(ssgseaScore)
  ssgseaOut.high=t(ssgseaOut)
  write.table(ssgseaOut.high,file="b-cells.txt",sep="\t",quote=F,col.names=F)
  
  rt2=data.frame(ssgseaOut) 
  annotation=data.frame(Type = factor(rep("B cells", ncol(rt2))))
  row.names(annotation) = names(rt2) 
  
  
  rt = data.frame(t(ssgseaOut))
  rt = data.frame(results[,1:10])  
  library(ggpubr)

  
  pdf("figure9-a1.pdf",width = 12,height = 9)
  par(mar = c(12,5,1,1))
  boxplot(rt,col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                     "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF"),par(las="2"),ylab = "Estimated Proportion")

  
  dev.off()
  
  tiff("figure9-a1.tiff",width = 3000,height = 2250,res = 300)
  par(mar = c(12,5,1,1))
  boxplot(rt,col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                     "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF"),par(las="2"),ylab = "Estimated Proportion")
  
  
  dev.off()
  
  
  
  rt2=data.frame(results)[,c(1:4,6,7,8)]
  library(pheatmap)  
  rt2 = data.frame(t(rt2))
  pdf("figure9-b11.pdf",width = 12,height = 9)
  pheatmap(rt2,
         
           fontsize_row=11,
          show_colnames = F,
           fontsize_col=3,
           color = colorRampPalette(c("green", "black", "red"))(50) )
  dev.off()
  tiff("figure9-b1.tiff",width = 3600,height =2000,res = 300)
  pheatmap(rt2, 
          # annotation=annotation, 
           #cluster_rows = FALSE,
           #cellwidth = 1.3,
           fontsize_row=11,
          show_colnames = F,
           fontsize_col=3,
           color = colorRampPalette(c("green", "black", "red"))(50) )
  dev.off()
  
  tcga.KIRC.RNA31 = data.frame(tcga.KIRC.RNA31)
 aa= data.frame(results[,1:10])
names(aa) = mydata[,1]
 aa$cc = aa$`Cytotoxic CD8 T cells1`+aa$`Cytotoxic CD8 T cells2`+aa$`Cytotoxic CD8 T cells3`+aa$`Cytotoxic CD8 T cells4`+aa$`Cytotoxic CD8 T cells5`+aa$`Cytotoxic CD8 T cells6`+aa$`Cytotoxic CD8 T cells7`+aa$`Cytotoxic CD8 T cells8`
 aa$cc = aa$cc/8
 aa$risk=as.vector(ifelse(aa$cc>median(aa$cc),"high","low"))
 tcga.KIRC.sur = data.frame(t(data.frame(t(tcga.KIRC.sur))))
 aa$OS = as.numeric(tcga.KIRC.sur[row.names(aa),1]) 
 aa$OS.time = as.numeric(tcga.KIRC.sur[row.names(aa),3]) 
 
 
 library(survival)
 
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
    pValue=signif(pValue,4)
    pValue=format(pValue, scientific = TRUE)
 }else{
    pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa)
 
 
 library(survminer)
 
 pdf(file="figure10.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,  xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 #aa = cbind(aa,results[row.names(aa),])
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells1`>median(aa$`Cytotoxic CD8 T cells1`),"high","low"))
 
 
 library(survival)
aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
    pValue=signif(pValue,4)
    pValue=format(pValue, scientific = TRUE)
 }else{
    pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-a.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,  xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-a.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells2`>median(aa$`Cytotoxic CD8 T cells2`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
   pValue=signif(pValue,4)
   pValue=format(pValue, scientific = TRUE)
 }else{
   pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-b.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,  xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-b.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells3`>median(aa$`Cytotoxic CD8 T cells3`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
   pValue=signif(pValue,4)
   pValue=format(pValue, scientific = TRUE)
 }else{
   pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-c.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,  xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-c.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells4`>median(aa$`Cytotoxic CD8 T cells4`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
   pValue=signif(pValue,4)
   pValue=format(pValue, scientific = TRUE)
 }else{
   pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-d.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,  xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-d.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells6`>median(aa$`Cytotoxic CD8 T cells6`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
    pValue=signif(pValue,4)
    pValue=format(pValue, scientific = TRUE)
 }else{
    pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-e.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-e.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 aa$risk=as.vector(ifelse(aa$`naive/memory CD8 T cells`>median(aa$`naive/memory CD8 T cells`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
    pValue=signif(pValue,4)
    pValue=format(pValue, scientific = TRUE)
 }else{
    pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-f.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-f.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE,xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 aa$risk=as.vector(ifelse(aa$`Exhausted CD8 T cells`>median(aa$`Exhausted CD8 T cells`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
    pValue=signif(pValue,4)
    pValue=format(pValue, scientific = TRUE)
 }else{
    pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-g.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-g.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells6`>median(aa$`Cytotoxic CD8 T cells6`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
   pValue=signif(pValue,4)
   pValue=format(pValue, scientific = TRUE)
 }else{
   pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-h.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-h.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells7`>median(aa$`Cytotoxic CD8 T cells7`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
   pValue=signif(pValue,4)
   pValue=format(pValue, scientific = TRUE)
 }else{
   pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-i.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-i.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells8`>median(aa$`Cytotoxic CD8 T cells8`),"high","low"))
 
 
 library(survival)
 aa1 = aa
 
 diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
 pValue = 1-pchisq(diff$chisq,df=1)
 if(pValue<0.001){
   pValue=signif(pValue,4)
   pValue=format(pValue, scientific = TRUE)
 }else{
   pValue=round(pValue,3)
 }
 
 fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
 
 
 library(survminer)
 
 pdf(file="figure10-j.pdf",
     width=7,
     height=6,onefile = FALSE)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 tiff(file="figure10-j.tiff",
      width=2333,
      height=2000,res = 300)
 ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
            legend.labs=c("High risk", "Low risk"), legend.title="Group",ggtheme = theme_bw(),
            palette=c("red", "blue"))
 dev.off()
 
 
 
 
 
 

 
 write.table( cbind(Name = row.names(cd8gse), DESCRIPTION = rep("na",29634), cd8gse[,cluster21[order(cluster21$pbmc1.seurat_clusters),2]]),"cd81.txt",sep = "\t",row.names = F, quote = F)
 cluster21 = data.frame(clusters2) 
write.csv(cluster21,"cluster21.csv")



ref = cd8gse[unique(sig.markers1$gene),]

clusters2 =  data.frame(pbmc1$seurat_clusters) 
clusters2$cc = row.names(clusters2)
for (i in 1:5509) {
   ref[,i] = as.numeric( ref[,i])
}

rrtt = data.frame(1:629)

for (i in 0:9) {
 ref1 = ref[,clusters2[which(clusters2$pbmc1.seurat_clusters == i),2]] 
# g = length(clusters2[which(clusters2$pbmc1.seurat_clusters == i),2])
 cc = apply(ref1, 1, mean)
   rrtt = cbind(rrtt,cc)
}
rrtt = rrtt[,-1]
names(rrtt) = mydata[,1]
 
 write.table(rrtt,"ref.txt",quote = F,sep = "\t")
 
 
  library(limma)
  
  data1 =  tcga.KIRC.RNA31[,60:583] 
  
  
  data1 = data1[rowMeans(data1)>0,]
  
  
  v <-voom(data1, plot = F, save.plot = F)
  out=v$E
  out=rbind(ID=colnames(out),out)
  write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)
  source("TMBimmune19.CIBERSORT.R")
  results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)
  write.csv(results,"results.csv")
  CIBERSORT.low = results[row.names(outTab11[which(outTab11$risk == "low"),]),]
  CIBERSORT.high = results[row.names(outTab11[which(outTab11$risk == "high"),]),]
  
  results = 
  
  tcga.KIRC.phenotype001 = tcga.KIRC.phenotype1[,c(41:43,97)]
  
  tcga.KIRC.phenotype001 = data.frame(tcga.KIRC.phenotype001)
  tcga.KIRC.phenotype001 = tcga.KIRC.phenotype001[intersect(row.names(results),row.names(tcga.KIRC.phenotype001)),]
  tcga.KIRC.phenotype001 =  tcga.KIRC.phenotype001[-which( tcga.KIRC.phenotype001$pathologic_M == ""),]
  tcga.KIRC.phenotype001 =  tcga.KIRC.phenotype001[-which( tcga.KIRC.phenotype001$tumor_stage.diagnoses == "not reported"),]
  
  tcga.KIRC.phenotype001$pathologic_M[which(tcga.KIRC.phenotype001$pathologic_M == "M1a")] = "M1"
  tcga.KIRC.phenotype001$pathologic_M[which(tcga.KIRC.phenotype001$pathologic_M == "M1b")] = "M1"
  tcga.KIRC.phenotype001$pathologic_T[which(tcga.KIRC.phenotype001$pathologic_T == "T1a")] = "T1"
  tcga.KIRC.phenotype001$pathologic_T[which(tcga.KIRC.phenotype001$pathologic_T == "T1b")] = "T1"
  tcga.KIRC.phenotype001$pathologic_T[which(tcga.KIRC.phenotype001$pathologic_T == "T2a")] = "T2"
  tcga.KIRC.phenotype001$pathologic_T[which(tcga.KIRC.phenotype001$pathologic_T == "T2b")] = "T2"
  tcga.KIRC.phenotype001$tumor_stage.diagnoses[which(tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage i"|
  tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage ia"|tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage ib")] = "stage I"
  tcga.KIRC.phenotype001$tumor_stage.diagnoses[which(tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage ii"|
                                                       tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage iia"|tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage iib")] = "stage II"
  tcga.KIRC.phenotype001$tumor_stage.diagnoses[which(tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage iiia"|tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage iiib")] = "stage III"
  tcga.KIRC.phenotype001$tumor_stage.diagnoses[which(tcga.KIRC.phenotype001$tumor_stage.diagnoses == "stage iv")] = "stage IV"
  
  tcga.KIRC.phenotype002 = cbind(tcga.KIRC.phenotype001,results[row.names(tcga.KIRC.phenotype001),])
  
 
  
  
  
  
   
 s1 = c(mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage I"),5]),mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage II"),5]),
 mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage III"),5]),mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage IV"),5]))
 
 s2 = c(mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage I"),8]),mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage II"),8]),
                                                                                                              mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage III"),8]),mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage IV"),8]))
 
 s3 = c(mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage I"),9]),mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage II"),9]),
                                                                                                              mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage III"),9]),mean(tcga.KIRC.phenotype002[which(tcga.KIRC.phenotype002$tumor_stage.diagnoses == "stage IV"),9])) 
 
s11 =   data.frame(s1,s2,s3)
 row.names(s11) =  
 names(s11) = names(tcga.KIRC.phenotype002)[c(5,8:9)]
 s11 = data.frame(ss = rep(c("stage I","stage II","stage III","stage IV"),3),`Cell proportions` = c(s1,s2,s3),dd = c(rep("Cytotoxic CD8 T cells1",4),
                                                                                                                   rep("naive/memory CD8 T cells",4),
                                                                                                                 rep("Exhausted CD8 T cells",4)) )
  
 
 names(tcga.KIRC.phenotype002)[4] = "Stage" 
 
 rt = tcga.KIRC.phenotype002
 
 library(ggpubr)
 
 clinical="pathologic_N" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$pathologic_N))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplotoxplot(rt, x=clinical, y="Exhausted CD8 T cells", color=clinical,
                   xlab=clinical,
                   ylab=paste("Exhausted CD8 T cells"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,".",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,".",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 library(ggpubr)
 rt = tcga.KIRC.phenotype002
 
 clinical="pathologic_M" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$pathologic_M))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplotoxplot(rt, x=clinical, y="Exhausted CD8 T cells", color=clinical,
                   xlab=clinical,
                   ylab=paste("Exhausted CD8 T cells"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,".",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,".",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 
 
 library(ggpubr)
 rt = tcga.KIRC.phenotype002
 
 clinical="pathologic_T" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$pathologic_T))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplotoxplot(rt, x=clinical, y="Exhausted CD8 T cells", color=clinical,
                   xlab=clinical,
                   ylab=paste("Exhausted CD8 T cells"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,".",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,".",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 
 
 library(ggpubr)
 rt = tcga.KIRC.phenotype002
 
 clinical="Stage" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$Stage))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplotoxplot(rt, x=clinical, y="Exhausted CD8 T cells", color=clinical,
                   xlab=clinical,
                   ylab=paste("Exhausted CD8 T cells"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,".",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,".",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 
 
 
 
 
 
 
 
 rt = tcga.KIRC.phenotype002
 
 library(ggpubr)
 
 clinical="pathologic_N" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$pathologic_N))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplotoxplot(rt, x=clinical, y="Cytotoxic CD8 T cells4", color=clinical,
                   xlab=clinical,
                   ylab=paste("Cytotoxic CD8 T cells4"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,"1.",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,"1.",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 library(ggpubr)
 rt = tcga.KIRC.phenotype002
 
 clinical="pathologic_M" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$pathologic_M))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplotoxplot(rt, x=clinical, y="Cytotoxic CD8 T cells4", color=clinical,
                   xlab=clinical,
                   ylab=paste("Cytotoxic CD8 T cells4"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,"1.",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,"1.",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 
 
 library(ggpubr)
 rt = tcga.KIRC.phenotype002
 
 clinical="pathologic_T" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$pathologic_T))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplot=ggboxplot(rt, x=clinical, y="Cytotoxic CD8 T cells4", color=clinical,
                   xlab=clinical,
                   ylab=paste("Cytotoxic CD8 T cells4"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,"1.",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,"1.",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 
 
 library(ggpubr)
 rt = tcga.KIRC.phenotype002
 
 clinical="Stage" 
 rt = rt[order(rt[,clinical]),]
 
 group=levels(factor(rt$Stage))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot
 boxplotoxplot(rt, x=clinical, y="Cytotoxic CD8 T cells4", color=clinical,
                   xlab=clinical,
                   ylab=paste("Cytotoxic CD8 T cells4"),
                   legend.title=clinical,
                   #title=paste0("Cancer: ",i),
                   add = "jitter")+ 
   stat_compare_means(comparisons = my_comparisons)
 pdf(file=paste0(clinical,"1.",".pdf"),width=7,height=5)
 print(boxplot)
 dev.off()
 tiff(file=paste0(clinical,"1.",".tiff"),width=1750,height=1250,res = 300)
 print(boxplot)
 dev.off()
 
 
 
 
 
 
 
 
 
 
 
  l proportions"
  
  pdf("zhe.pdf",width = 10,height = 5)
  ggplot(s11,aes(ss,`Cell proportions`,colour = dd ,group = dd,shape=dd))+geom_line()+geom_point(size=4)+xlab("")+
    theme(legend.title=element_text(face="italic", family="Times", colour="white",
                                    size=14))
  dev.off()
  tiff("zhe.tiff",width = 2600,height = 1300,res = 300)
  ggplot(s11,aes(ss,`Cell proportions`,colour = dd ,group = dd,shape=dd))+geom_line()+geom_point(size=4)+xlab("")+
    theme(legend.title=element_text(face="italic", family="Times", colour="white",
                                    size=14))
  dev.off()
  
  
  cd8gse
  write.csv(cd8gse,"cd8gse.csv",quote = F)
  
  library(clusterProfiler)
  library(biomaRt)
  library(org.Hs.eg.db)
  mart <- useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  

  eg1 <- bitr(sig.markers1$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  eg1 <- eg1[!duplicated(eg1[,2]),]
  go <- enrichGO(eg1[,2], OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, readable = T)
  
  kegg <- enrichKEGG(eg1[,2], organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 5,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
  
  pbmc.markers2 = sig.markers1[,6:7]
  names(pbmc.markers2)[2] =  "SYMBOL"
  de_gene_clusters = merge(pbmc.markers2,eg1)
  de_gene_clusters = de_gene_clusters[,2:3]
  de_gene_clusters = unique(de_gene_clusters)
  de_gene_clusters = data.frame(cluster = as.numeric(de_gene_clusters[,1]),ENTREZID = de_gene_clusters[,2])
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "1")] = "Cytotoxic CD8 T cells1"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "2")] = "Cytotoxic CD8 T cells2"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "3")] = "Cytotoxic CD8 T cells3"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "4")] = "Cytotoxic CD8 T cells4"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "5")] = "Cytotoxic CD8 T cells5"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "6")] = "Cytotoxic CD8 T cells6"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "7")] = "naive/memory CD8 T cells"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "8")] = "Exhausted CD8 T cells"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "9")] = "Cytotoxic CD8 T cells7"
  de_gene_clusters$cluster[which(de_gene_clusters$cluster == "10")] = "Cytotoxic CD8 T cells8"
  
  
  list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, 
                                 de_gene_clusters$cluster)
  
  #######
  # ???????ݽṹlist?? ?s <mpareCluster(list_de_gene_clusters,
                                fun="enrichGO",
                                OrgDb="org.Hs.eg.db", 
                                ont          = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff=0.05,
                                readable = T)
  
  # ???ӻ?
  pdf('fidf',width = 20,height = 8)

  dotplot(formula_res, showCategory=5) +theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
  dev.off()
  tiff('figure12-a.tiff',width = 5600,height = 2400,res = 300)
  dotplot(formula_res, showCategory=5) +theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
  dev.off()
  
  formula_res1 <- compareCluster(list_de_gene_clusters,
                                fun="enrichKEGG",
                                organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                                minGSSize = 5,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
  
  # ???ӻ?
  pdf('figudf',width = 18,height = 8)
  dotplot(formula_res1, showCategory=5) +theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
  dev.off()
  tiff('figure12-b.tiff',width = 5400,height = 2400,res = 300)
  dotplot(formula_res1, showCategory=5) +theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
  dev.off()
  
  
  write.csv(formula_res ,"formula_res.csv")
  write.csv(formula_res1 ,"formula_res1.csv")
  
  library(GSVA)
  
  library(GSEABase) 
  library(clusterProfiler)
  library(DOSE)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(stringr)
  options(stringsAsFactors = F)
  
  
  expr <- cd8gse[,clusters2[,2]]
  expr$ID <- rownames(expr)
  s2e <- bitr(expr$ID, 
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)#??
  expr <- inner_n(expr,s2e,by=c("ID"="SYMBOL"))
  rownames(expr) <- expr$ENTREZID
  expr <- expr[,-5755]
  expr <- expr[,-5754]
  meta <- pbmc1@meta.data[,c("celltype1")]
  
  kegggmt <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
  kegg_list = split(kegggmt$gene, kegggmt$term)
  
  expr=as.matrix(expr)
  kegg21 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)

  kegg21 = data.frame(kegg21)
  normalize = function(x){
    return((x-min(x))/(max(x)-min(x)))}
  #??ssGSEA score???н?
  kegg111111
  kegg3=normalize(kegg21)
  data = data.frame(1:186)
  kegg3=normalize(kegg2)
  data = data.frame(1:52)
  data$GC <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 0),2]], 1, mean)
  data$IGA1 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 1),2]], 1, mean)
  data$IGA2 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 2),2]], 1, mean)
  data$IGA3 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 3),2]], 1, mean)
  data$IGA4 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 4),2]], 1, mean)
  data$IGA5 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 5),2]], 1, mean)
  data$IGA6 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 6),2]], 1, mean)
  data$IGA7 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 7),2]], 1, mean)
  data$IGA8 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 8),2]], 1, mean)
  data$IGA9 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 9),2]], 1, mean)

 
  data = data[,-1]
  row.names(data) = row.names(kegg3)
  names(data) = mydata[,1]
  write.csv(data,"kegg2.csv")
  write.csv(data,"kegg21.csv")
  library(vioplot)
  pdf("vivi.pdf",width = 10,height = 6)
  par(mar = c(12,5,1,1),las="2")
  
  vioplot(data,
          names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                   "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))
  
  dev.off()
  tiff("vivi.tiff",width = 3000,height = 1800,res = 300)
  par(mar = c(12,5,1,1),las="2")
  vioplot(data,
          names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                   "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))
  
  dev.off()
  library(pheatmap)
  
  tiff("figure13-a.tiff",width = 7000,height = 2666,res = 300)
  pheatmap(t(data),cluster_rows = F,scale = "column",
           cluster_cols = F,legend = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  pdf("figure13-a.pdf",height = 10,width = 25)
  pheatmap(t(data),cluster_rows = F,scale = "column",
           cluster_cols = F,legend = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()

  
  pdf("vivi1.pdf",width = 10,height = 6)
  par(mar = c(12,5,1,1),las="2")
  
  vioplot(data,
          names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                   "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))
  
  dev.off()
  tiff("vivi1.tiff",width = 3000,height = 1800,res = 300)
  par(mar = c(12,5,1,1),las="2")
  vioplot(data,
          names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                   "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))
  
  dev.off()
  
  
  tiff("figure13-d.tiff",width = 7000,height = 2666,res = 300)
  pheatmap(t(data),cluster_rows = F,scale = "column",
           cluster_cols = F,legend = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  pdf("figure13-d.pdf",height = 10,width = 25)
  pheatmap(t(data),cluster_rows = F,scale = "column",
           cluster_cols = F,legend = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  
   library(dplyr)
 gty <- FindAllMarkers(object = pbmc1,
                                  min.pct = 0,    logfc.threshold = 0,
      
                                  only.pos = FALSE)
 
 gty1 <- FindMarkers(object = pbmc1, ident.1 = 0, min.pct = 0,logfc.threshold = 0)
 
 gty$SYMBOL = gty$gene
 expr <- inner_join(gty,s2e,by=c("gene"="SYMBOL"))
 expr = expr[,c(7,1:6,8)]
 de_gene_clusters = expr
  data = data[,-1]
  names(data) = mydata[,1]
  row.names(data) = row.names(kegg2)
  write.csv(data,"data.gsva.csv")
  write.csv(data,"data.gsva1.csv")
  
  gty1 <- FindMarkers(object = pbmc1, ident.1 = 0, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty1
  de_gene_clusters$gene = row.names(de_gene_clusters)

  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  geneset <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")  

  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps = 1e-10,
               pvalueCutoff = 1,verbose=F)

  y=data.frame(egmt)
  
  gty2 <- FindMarkers(object = pbmc1, ident.1 = 1, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty2
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y1=data.frame(egmt)
  
  gty3 <- FindMarkers(object = pbmc1, ident.1 = 2, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty3
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y2=data.frame(egmt)
 
   gty4 <- FindMarkers(object = pbmc1, ident.1 = 3, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty4
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y3=data.frame(egmt)
  
  
  gty5 <- FindMarkers(object = pbmc1, ident.1 = 4, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty5
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y4=data.frame(egmt)
  
  gty6 <- FindMarkers(object = pbmc1, ident.1 = 5, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty6
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y5=data.frame(egmt)
  
  gty7 <- FindMarkers(object = pbmc1, ident.1 = 6, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty7
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y6=data.frame(egmt)
  
  gty8 <- FindMarkers(object = pbmc1, ident.1 = 7, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty8
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y7=data.frame(egmt)
  
  gty9 <- FindMarkers(object = pbmc1, ident.1 = 8, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty9
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y8=data.frame(egmt)
  
  
  gty10 <- FindMarkers(object = pbmc1, ident.1 = 9, min.pct = 0,logfc.threshold = 0)
  
  de_gene_clusters = gty10
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  y9=data.frame(egmt)
  
 ggggggg = rep(c(mydata[,1]),c(12,12,12,12,12,12,12,12,12,12))
  
 write.csv(ggggggg,"gggggg.csv") 
 yyy = rbind(y,y1,y2,y3,y4,y5,y6,y7,y8,y9)
 yyy$cluster = rep(mydata[,1],rep(186,10))
  write.csv(yyy,"yyy.csv")
 
  
yi =  unique(  c(row.names(y[1:5,c(5,6)]),row.names(y1[1:5,c(5,6)]),row.names(y2[1:5,c(5,6)]),row.names(y3[1:5,c(5,6)]),row.names(y4[1:5,c(5,6)]),row.names(y5[1:5,c(5,6)]),
            row.names(y6[1:5,c(5,6)]),row.names(y7[1:5,c(5,6)]),row.names(y8[1:5,c(5,6)]),row.names(y9[1:5,c(5,6)]))) 

 nespv =  rbind(y[yi,c(1,5,6)],y1[yi,c(1,5,6)],y2[yi,c(1,5,6)],y3[yi,c(1,5,6)],y4[yi,c(1,5,6)],y5[yi,c(1,5,6)],y6[yi,c(1,5,6)],y7[yi,c(1,5,6)],y8[yi,c(1,5,6)],y9[yi,c(1,5,6)])
  
 nespv$cluster = rep(mydata[,1],c(26,26,26,26,26,26,26,26,26,26))
  
  c= c(9,1,3,10,3,4,1,0,3,4)
  
 c = merge(y, y1,by = c("ID"), all = T)
  c = merge(c,y2,by = c("ID"), all = T)
  c = merge(c,y3,by = c("ID"), all = T)
  c = merge(c,y4,by = c("ID"), all = T)
  c = merge(c,y5,by = c("ID"), all = T)
  c = merge(c,y6,by = c("ID"), all = T)
  c = merge(c,y7,by = c("ID"), all = T)
  c = merge(c,y8,by = c("ID"), all = T)
  c = merge(c,y9,by = c("ID"), all = T)
  
  c[is.na(c)] = 0
  
write.csv(c,"cp2.csv")
  nes = read.csv("nes.csv",header = T,row.names = 1)
  
  data = cbind(y$NES,y1$NES,y2$NES,y3$NES,y4$NES,y5$NES,y6$NES,y7$NES,y8$NES,y9$NES)
  
  data = data[,-1]
  row.names(data) = row.names(kegg3)
  
  row.names(data) = row.names(y)
  colnames(data) = mydata[,1]
 
  

  library(vioplot)
  library()
pdf("vivi.pdf",width = 10,height = 6)
par(mar = c(12,5,1,1),las="2")

vioplot(y$enrichmentScore,y1$enrichmentScore,y2$enrichmentScore,
        y3$enrichmentScore,y4$enrichmentScore,y5$enrichmentScore,y6$enrichmentScore,y7$enrichmentScore,y8$enrichmentScore,y9$enrichmentScore,
        names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                 "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))

 dev.off()
 tiff("vivi.tiff",width = 3000,height = 1800,res = 300)
 par(mar = c(12,5,1,1),las="2")
 vioplot(y$enrichmentScore,y1$enrichmentScore,y2$enrichmentScore,
         y3$enrichmentScore,y4$enrichmentScore,y5$enrichmentScore,y6$enrichmentScore,y7$enrichmentScore,y8$enrichmentScore,y9$enrichmentScore,
         names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                  "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))

 dev.off()

 pdf("vivi1.pdf",width = 10,height = 6)
 par(mar = c(12,5,1,1),las="2")
 
 vioplot(y10$enrichmentScore,y11$enrichmentScore,y12$enrichmentScore,
         y13$enrichmentScore,y14$enrichmentScore,y15$enrichmentScore,y16$enrichmentScore,y17$enrichmentScore,y18$enrichmentScore,y19$enrichmentScore,
         names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                  "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))
 
 dev.off()
 tiff("vivi1.tiff",width = 3000,height = 1800,res = 300)
 par(mar = c(12,5,1,1),las="2")
 vioplot(y10$enrichmentScore,y11$enrichmentScore,y12$enrichmentScore,
         y13$enrichmentScore,y14$enrichmentScore,y15$enrichmentScore,y16$enrichmentScore,y17$enrichmentScore,y18$enrichmentScore,y19$enrichmentScore,
         names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ,
                                                  "#0066FFFF", "#3300FFFF", "#CC00FFFF" ,"#FF0099FF" ))
 
 dev.off()
 
 
 
 
names(nespv)[1] = "pathway"
write.csv(nespv,"nespvkegg.csv") 

 tiff(file="go1dot11.tiff",width = 3000, height = 3000,res=300)
 ggplot(data=nespv,aes(cluster,pathway,group = cluster))+  geom_point()+geom_point(aes(size = -log10(pvalue),color=NES))+
   coord_flip() +scale_color_gradientn(colours = c("blue","white", "red"))+theme(legend.position="top")+xlab(NULL)+ylab(NULL)+
   theme(axis.text.x = element_text(angle=75, hjust = 1,vjust=1))
 
   
   
 dev.off()
 pdf(file="go1dot11.pdf",width = 12, height = 12)
 ggplot(data=nespv,aes(cluster,pathway,group = cluster))+  geom_point()+geom_point(aes(size = -log10(pvalue),color=NES))+
   coord_flip() +scale_color_gradientn(colours = c("blue","white", "red"))+theme(legend.position="top")+xlab(NULL)+ylab(NULL)+
   theme(axis.text.x = element_text(angle=75, hjust = 1,vjust=1))
 
 dev.off()
 
  dotplot(egmt2,split=".sign")+facet_grid(~.sign)
  library(GSVA)
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  m_df<- msigdbr(species = "Homo sapiens", category = "C2")
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  pbmc.genes <- wilcoxauc(pbmc1, 'seurat_clusters')
  
  pbmc.genes %>%
    dplyr::filter(group == "0") %>%
    arrange(desc(logFC), desc(auc)) %>%
    head(n = 10)
  
  fgseaRes <- fgsea(pathways = fgsea_sets,
                    stats = ranks ,
                    minSize=5,
                    maxSize=500,
                    nperm=10000)
  gsea_genes<-pbmc.markers1 %>%
    arrange(desc(myAUC), desc(avg_diff)) %>%
    dplyr::select(gene,avg_diff)
  
  expr <- cd8gse
  expr$ID <- rownames(expr)
  s2e <- bitr(expr$ID, 
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)#??
  expr <- inner_join(expr,s2e,by=c("ID"="SYMBOL"))
  rownames(expr) <- expr$ENTREZID
  expr <- expr[,-3070]
  expr <- expr[,-3069]
  meta <- pbmc1@meta.data[,c("celltype1")]
  
  kegggmt <- read.gmt("immune1.gmt")
  colnames(kegggmt)
  kegg_list = split(kegggmt$gene, kegggmt$term)
  
  expr=as.matrix(expr)
  kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)
  kegg2 = data.frame(kegg2)
  normalize = function(x){
    return((x-min(x))/(max(x)-min(x)))}
  #??ssGSEA score???н?
  kegg111 = kegg3=normalize(kegg2)
  data = data.frame(1:52)
  data$GC <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 0),2]], 1, mean)
  data$IGA1 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 1),2]], 1, mean)
  data$IGA2 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 2),2]], 1, mean)
  data$IGA3 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 3),2]], 1, mean)
  data$IGA4 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 4),2]], 1, mean)
  data$IGA5 <- apply(kegg3[,clusters2[which(clusters2$pbmc1.seurat_clusters == 5),2]], 1, mean)
  
  data = data[,-1]
  names(data) = mydata[,1]
  row.names(data) = row.names(kegg2)
  
  tiff("nes1.tiff",width = 2000,height = 3000,res = 300)
  pheatmap(data,cluster_rows = F,scale = "row",
           cluster_cols = F,legend = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  pdf("nes1.pdf",height = 9.5,width = 6)
  pheatmap(data,cluster_rows = F,scale = "row",
           cluster_cols = F,legend = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()

  
  geneset <- read.gmt("h.all.v7.4.entrez (1).gmt") 
  s2e = bitr(geneset$gene, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
  expr = geneset
  names(expr)[2] = "SYMBOL"
  expr <- inner_join(expr,s2e,by=c("ID"="SYMBOL")) 
  
  geneset <- read.gmt("immune1.gmt")
  
  de_gene_clusters = gty1
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  

  y10=data.frame(egmt)

  
  
  de_gene_clusters = gty2
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  #egmt2<- setReadable(egmt,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
  y11=data.frame(egmt)
  
  de_gene_clusters = gty3
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y12=data.frame(egmt)
  
  de_gene_clusters = gty4
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y13=data.frame(egmt)
  
  
  de_gene_clusters = gty5
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y14=data.frame(egmt)
  
  de_gene_clusters = gty6
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y15=data.frame(egmt)
  
  
  de_gene_clusters = gty7
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y16=data.frame(egmt)
  
  
  
  de_gene_clusters = gty8
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y17=data.frame(egmt)
  
  de_gene_clusters = gty9
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y18=data.frame(egmt)
  
  de_gene_clusters = gty10
  de_gene_clusters$gene = row.names(de_gene_clusters)
  
  geneList =de_gene_clusters$avg_log2FC
  names(geneList) = de_gene_clusters$gene
  geneList = sort(geneList,decreasing = T)
  
  egmt <- GSEA(geneList, TERM2GENE=geneset,minGSSize = 3,
               maxGSSize = 1000,
               eps =  1e-10,
               pvalueCutoff = 1,verbose=F)
  
  
  y19=data.frame(egmt)
  
  
  
  yyyyy = rbind(y10,y11,y12,y13,y14,y15,y16,y17,y18,y19)
  yyyyy$cluster = rep(mydata[,1],rep(52,10))
  write.csv(yyyyy,"yyyyy.csv")
  
  
  c = merge(y10, y11,by = c("ID"), all = T)
  c = merge(c,y12,by = c("ID"), all = T)
  c = merge(c,y13,by = c("ID"), all = T)
  c = merge(c,y14,by = c("ID"), all = T)
  c = merge(c,y15,by = c("ID"), all = T)
  c = merge(c,y16,by = c("ID"), all = T)
  c = merge(c,y17,by = c("ID"), all = T)
  c = merge(c,y18,by = c("ID"), all = T)
  c = merge(c,y19,by = c("ID"), all = T)
  
  c[is.na(c)] = 0
  
  write.csv(c,"cp21.csv")
  
 write.csv(yyyyy,"yyyyy.csv") 
  
 pdf("vivi1.pdf",width = 10,height = 6)
 par(mar = c(12,5,1,1),las="2")
 
 vioplot(y10$enrichmentScore,y11$enrichmentScore,y2$enrichmentScore,
         y9$enrichmentScore,y10$enrichmentScore,y11$enrichmentScore,
         names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ))
 
 dev.off()
 tiff("vivi1.tiff",width = 3000,height = 1800,res = 300)
 par(mar = c(12,5,1,1),las="2")
 vioplot(y6$enrichmentScore,0,y8$enrichmentScore,
         y9$enrichmentScore,y10$enrichmentScore,y11$enrichmentScore,
         names = mydata[,1], ylab="Score",col = c("#FF0000FF" ,"#FF9900FF" ,"#CCFF00FF","#33FF00FF" ,"#00FF66FF" ,"#00FFFFFF" ))
 
 dev.off()
 
 dotplot(nes)
 nespv = read.csv("nespv3.csv",header = T)
 
 yi =  unique(  c(row.names(y10[1:5,c(5,6)]),row.names(y11[1:5,c(5,6)]),row.names(y12[1:5,c(5,6)]),row.names(y13[1:5,c(5,6)]),row.names(y14[1:5,c(5,6)]),row.names(y15[1:5,c(5,6)]),
                  row.names(y16[1:5,c(5,6)]),row.names(y17[1:5,c(5,6)]),row.names(y18[1:5,c(5,6)]),row.names(y19[1:5,c(5,6)]))) 
 
 nespv =  rbind(y10[yi,c(1,5,6)],y11[yi,c(1,5,6)],y12[yi,c(1,5,6)],y13[yi,c(1,5,6)],y14[yi,c(1,5,6)],y15[yi,c(1,5,6)],y16[yi,c(1,5,6)],y17[yi,c(1,5,6)],y18[yi,c(1,5,6)],y19[yi,c(1,5,6)])
 
 nespv$class = rep(mydata[,1],c(21,21,21,21,21,21,21,21,21,21))
 names(nespv)[1] = "pathway" 
 write.csv(nespv,"nespvhall.csv") 
 tiff(file="go1dot1.tiff",width = 3000, height = 2000,res=300)
 ggplot(data=nespv,aes(class,pathway,group = class))+  geom_point()+geom_point(aes(size = -log10(pvalue),color=NES))+
   coord_flip() +scale_color_gradientn(colours = c("blue","white", "red"))+theme(legend.position="top")+xlab(NULL)+ylab(NULL)+
   theme(axis.text.x = element_text(angle=75, hjust = 1,vjust=1))
 
 
 
 dev.off()
 pdf(file="go1dot1.pdf",width = 12, height = 8)
 ggplot(data=nespv,aes(class,pathway,group = class))+  geom_point()+geom_point(aes(size = -log10(pvalue),color=NES))+
   coord_flip() +scale_color_gradientn(colours = c("blue","white", "red"))+theme(legend.position="top")+xlab(NULL)+ylab(NULL)+
   theme(axis.text.x = element_text(angle=75, hjust = 1,vjust=1))
 
 dev.off()
  
 monocle.matrix=cd8gse
 monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
 write.table(monocle.matrix,file="07.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
 monocle.matrix = monocle.matrix[,-2]
 library(clusterProfiler)
 library(biomaRt)
 library(org.Hs.eg.db)
 mart <- useMart("ensembl")
 mart <- useDataset("hsapiens_gene_ensembl", mart)
 
 
 
 eg1 <- bitr(row.names(monocle.matrix), fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
 eg1 = na.omit(eg1)
 eg1 <-eg1[!duplicated(eg1[,2]),]
 eg1 <-eg1[!duplicated(eg1[,3]),]
 
 names(eg1)[1]="id"
 row.names(eg1) = eg1[,1]
 monocle.matrix = monocle.matrix[row.names(eg1),]
 row.names(monocle.matrix) = eg1[,3]
 
 monocle.matrix=as.matrix(pbmc1@assays$RNA@data)
 monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
 write.table(monocle.matrix,file="07.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
 monocle.sample=as.matrix(pbmc1@meta.data)
 monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
 write.table(monocle.sample,file="07.monocleSample.txt",quote=F,sep="\t",row.names=F)
 monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
 monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
 write.table(monocle.geneAnn,file="07.monocleGene.txt",quote=F,sep="\t",row.names=F)
 write.table(singler$other,file="07.monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
 pbmc.markers1 <- FindAllMarkers(object = pbmc1,
                                 min.pct = 0.25,        
                                 only.pos = FALSE)
 write.table(pbmc.markers1,file="07.monocleMarkers.txt",sep="\t",row.names=F,quote=F)
 library(Rcpp)
 library(monocle)
  monocle.matrix=read.table("07.monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
 monocle.sample=read.table("07.monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
 monocle.geneAnn=read.table("07.monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
 
 marker=read.table("07.monocleMarkers.txt",sep="\t",header=T,check.names=F)
 
 #??Seurat????ת??Ϊmocle??�ix(monocle.matrix), 'sparseMatrix')
 pd<-new("AnnotatedDataFrame", data = monocle.sample)
 fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
 cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
 
 #??????һ???????????names(pData(cds))=="seurat_clusters"]="Cluster"
 pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
 
 #????ϸ??????????
 cle("07.monocleClusterAnn.txt",header=F,sep="\t",check.names=F)
 clusterAnn=as.character(clusterRt[,2])
 names(clusterAnn)=paste0("cluster",clusterRt[,1])
 pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)
 
 #αʱ??????????
 cdseFactors(cds)
 cds <- estimateDispersions(cds)
 
 cds <- setOrderingFilter(cds,c("CCR7","TCF7","LEF1","SELL","PRF1", "IFNG", "GNLY", "NKG7", "GZMB","GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7",
                                "CTLA4", "PDCD1", "LAG3", "HAVCR2", 
                                "TIGIT"))
 
 cds <- setOrderingFilter(cds,marker$gene)
 plot_ordering_genes(cds)
 cds <- reduceDimension(cds, max_components = 2,reduction_method = "DDRTree")
 cds <- orderCells(cds)
 pdf(file="cluster.trajectory.pdf",width=6.5,height=6)
 plot_cell_trajectory(cds, color_by = "Cluster")
 dev.off()
 tiff(file="cluster.trajectory.tiff",width=1950,height = 1800,res = 300)
 plot_cell_trajectory(cds, color_by = "Cluster")
 dev.off()
 pdf(file="cellType.trajectory.pdf",width=10,height=9.5)
 plot_cell_trajectory(cds,color_by = "cell_type2")
 dev.off()
 
 cds <- orderCells(cds, root_state = 7)
 plot_cell_trajectory(cds,color_by="State", size=1,show_backbone=TRUE) 
 
 
 tiff("cellType.trajectory.tiff",width=3000,height=2900,res = 300)
 plot_cell_trajectory(cds,color_by = "cell_type2")
 dev.off()
 plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
 
state = pData(cds)
  plot_cell_trajectory(cds,color_by = "State")
 
 plot_cell_clusters(cds)
 
 plot_genes_in_pseudotime(cds, color_by="Hours")
 cg = c("CCR7","TCF7","LEF1","SELL","PRF1", "IFNG", "GNLY", "NKG7", "GZMB","GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7",
        "CTLA4", "PDCD1", "LAG3", "HAVCR2", 
        "TIGIT")
 
 for (i in 1:25) {
   
   
  pp = plot_genes_in_pseudotime(cds[cg[i],],
                            color_by = "cell_type2")
   pdf(file=paste("cellType77",i,".pdf"),width=10,height=3.33)
   
   print(pp)
   
   dev.off()
   tiff(file=paste("cellType77",i,".tiff"),width=3000,height=1000,res = 300)
   print(pp)
   
   
   dev.off()
 }
 
 
 monocle.sample = data.frame(monocle.sample)
monocle.sample$seurat_clusters[which(monocle.sample$seurat_clusters == "0")] = "Cytotoxic CD8 T cells1"
monocle.sample$seurat_clusters[which(monocle.sample$seurat_clusters == "1")] = "Cytotoxic CD8 T cells2"
monocle.sample$seurat_clusters[which(monocle.sample$seurat_clusters == "2")] = "Cytotoxic CD8 T cells3"
monocle.sample$seurat_clusters[which(monocle.sample$seurat_clusters == "3")] = "naive/memory CD8 T cells"
monocle.sample$seurat_clusters[which(monocle.sample$seurat_clusters == "4")] = "Exhausted CD8 T cells"
monocle.sample$seurat_clusters[which(monocle.sample$seurat_clusters == "5")] = "Cytotoxic CD8 T cells4"
monocle.sample = monocle.sample[,c(1,7)]
 
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.sample,file="monocle.sample111.txt",quote=F,sep="\t",row.names=F)
 
  
inters = read.table("significant_means.txt",sep = "\t",header = T,row.names = 1)


write.csv(ao,"ao.csv")
write.csv(ao1,"ao1.csv")
write.csv(ao,"aopheatmap.csv")

ao1 = c()
for (i in 1:592) {
  ap1 = length(na.omit(t(inters[i,12:111])))
  ao1 = c(ao1,ap1)
}

ao1 = data.frame(ao1,pair = inters[,1])


ao = c()
for (i in 12:length(names(inters))) {
 ap = length(na.omit(inters[,i]))
  ao = c(ao,ap)
}

ao = data.frame(rep(mydata[order(mydata[,1]),1],c(10,10,10,10,10,10,10,10,10,10)),rep(mydata[order(mydata[,1]),1],10),ao)

dim(ao) <- c(10, 10)
ao = data.frame(ao)
names(ao) = mydata[order(mydata[,1]),1]
row.names(ao) = mydata[order(mydata[,1]),1]
  library(pheatmap)
  inters = read.csv("tretete.csv",header = T,row.names = 1)

  tiff("inter1.tiff",width = 2250,height = 1500,res = 300)
  pheatmap(ao,cluster_rows = F,
           cluster_cols = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  pdf("inter1.pdf",height = 6,width = 9)
  pheatmap(ao,cluster_rows = F,
           cluster_cols = F,color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  
write.csv(pbmc.markers1,"pbmc.markers1.csv")

chayi = sig.markers1[which(sig.markers1$cluster == 7),7]

chayi = unique(chayi)

chayi1 = pbmc.markers1[which(pbmc.markers1$cluster == 0),7]

chayi1 = unique(chayi1)

pa = read.csv("table_s1_ligand-receptor_pair_list_bbz040.csv",header = T,row.names = 1)

write.csv(state,"state1.csv")

sig1 = intersect(chayi,pa[,2])
sig2 = intersect(chayi,pa[,3])
sig11 = intersect(chayi1,pa[,2])
sig21 = intersect(chayi1,pa[,3])


sig1 = data.frame(sig1)
names(sig1) = "Ligand.Symbol"

sig2 = data.frame(sig2)
names(sig2) = "Receptor.Symbol"


trrust_rawdata = read.csv("trrust_rawdata.human.csv",header = F)


sig.markers11 = sig.markers1[which(abs(sig.markers1$avg_log2FC) > 1),]
sig.markers11 = sig.markers11[which(sig.markers11$cluster == 7),]
write.csv(sig.markers11,"sigmarker11.csv")
write.csv(union(union(sig1,sig2),intersect(trrust_rawdata[,1],chayi)),"exactr.csv")

write.csv(union(union(sig11,sig21),intersect(trrust_rawdata[,1],chayi1)),"exactr1.csv")




aa = data.frame(t(tcga.KIRC.RNA31[intersect(unique(c(union(union(sig1,sig2),intersect(trrust_rawdata[,1],chayi)),row.names(dif2))),row.names(tcga.KIRC.RNA31)),]) )  
aa = data.frame(t(tcga.KIRC.RNA31[union(sig.markers11[,7],row.names(dif2)),]) )  

aa1 = aa[intersect(row.names(tcga.KIRC.phenotype1),row.names(aa)),]
aa = aa[row.names(aa1),]


tcga.KIRC.sur = data.frame(t(data.frame(t(tcga.KIRC.sur))))
aa$OS = as.numeric(tcga.KIRC.sur[row.names(aa),1]) 
aa$OS.time = as.numeric(tcga.KIRC.sur[row.names(aa),3]) 

library(survival)
rt = aa[,c(66,67,1:65)]

outTab2=data.frame()


for(gene in colnames(rt[,3:ncol(rt)])){
  a = rt[,gene] <= median(rt[,gene])
  diff=survdiff(Surv(OS.time, OS) ~ c(rt[,gene] <= median(rt[,gene])),data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab2=rbind(outTab2,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,4)
  #pVa?ue=format(pValue, scientific = TRUE)
  
  fit <- survfit(Surv(OS.time, OS) ~ a, data = rt)
  summary(fit)}


aa = data.frame(t(tcga.KIRC.RNA31[intersect(unique(c(union(union(sig11,sig21),intersect(trrust_rawdata[,1],chayi1)),row.names(dif2))),row.names(tcga.KIRC.RNA31)),]) )  
aa1 = aa[intersect(row.names(tcga.KIRC.phenotype1),row.names(aa)),]
aa = aa[row.names(aa1),]


tcga.KIRC.sur = data.frame(t(data.frame(t(tcga.KIRC.sur))))
aa$OS = as.numeric(tcga.KIRC.sur[row.names(aa),1]) 
aa$OS.time = as.numeric(tcga.KIRC.sur[row.names(aa),3]) 

library(survival)
rt = aa[,c(183,184,1:182)]

outTab21=data.frame()


for(gene in colnames(rt[,3:ncol(rt)])){
  a = rt[,gene] <= median(rt[,gene])
  diff=survdiff(Surv(OS.time, OS) ~ c(rt[,gene] <= median(rt[,gene])),data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab21=rbind(outTab21,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,4)
  #pVa?ue=format(pValue, lsientific = TRUE)
  
  fit <- survfit(Surv(OS.time, OS) ~ a, data = rt)
  summary(fit)}




library(pheatmap)

pdf("tf.pdf",width = 8,height = 10)
pheatmap(rrtt1[intersect(trrust_rawdata[,1],chayi),],show_rownames = T, show_colnames = T,scale = "row",
         angle_col = c("45"))
dev.off()
tiff("tf.tiff",width = 2400,height = 3000,res = 300)
pheatmap(rrtt1[intersect(trrust_rawdata[,1],chayi),],show_rownames = T, show_colnames = T,scale = "row",
         angle_col = c("45") )
dev.off()


for (i in 1:3068) {
  cd8gse[,i] = as.numeric(cd8gse[,i])
}

pair = ao1[-which(ao1$ao1 == 0),]
c = as.numeric( pair[,1])
pdf("pair111.pdf",width = 15,height = 7.5)
par(mar = c(13,5,1,1),las="2")
barplot(c, names.arg = pair[,2],ylab = "Ligand-receptor Pairs",xlab = "",col = "turquoise")
dev.off()
tiff("pair111.tiff",width = 4000,height = 2000,res = 300)
par(mar = c(13,5,1,1),las="2")

barplot(c, names.arg = pair[,2],ylab = "Ligand-receptor Pairs",xlab = "",col = "turquoise")
dev.off()


rrtt1 = data.frame(1:29634)

for (i in 0:9) {
 cd8gse11 = cd8gse[,clusters2[which(clusters2$pbmc1.seurat_clusters == i),2]] 
  # g = length(clusters2[which(clusters2$pbmc1.seurat_clusters == i),2])
  cc = apply(cd8gse11, 1, mean)
  rrtt1 = cbind(rrtt1,cc)
}
rrtt1 = rrtt1[,-1]
names(rrtt1) = mydata[,1]



string4 = read.table("string_interactions (23).tsv",header = F)

intersect(string4[,1],outTab2[which(outTab2$pvalue < 0.05),1])


  aa$risk=as.vector(ifelse(aa$`Cytotoxic CD8 T cells1`>median(aa$`Cytotoxic CD8 T cells1`),"high","low"))
  
 qqq = aa[which(aa$risk == "low"),]
 qqq = rbind(aa[which(aa$risk == "low"),],aa[which(aa$risk == "high"),])
 
 
 
 
  library(limma)
  
  group_list = c(rep("low",262), rep("high",262))
  
  design <- model.matrix(~0+factor(group_list))
  design
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(tcga.KIRC.RNA31[,row.names(qqq)])
  
  design
  
  contrast.matrix<-makeContrasts(high-low,levels = design)
  contrast.matrix ##?????????????????? <- lmFit(tcga.KIRC.RNA31[,row.names(qqq)],design)
  ##step2 ???ݶԱ?ģ?ͽ??в?�.fit(fit, contrast.matrix) 
  ##step3 ??Ҷ˹????
  fit2 <-t2) 
  ##step4 ???????л????ļ?????Table(fit2, coef=1, n=Inf)
  ##step5 ??P.Value????ɸѡ????ȫ????t1[tempOutput1[, "P.Value"]<0.05,]
  dif1<- dif1[abs(dif1[,"logFC"])>1,]
  
  
  
  library(pheatmap)
  map3 <- tcga.KIRC.RNA31[row.names(dif1),row.names(qqq)]
  
  annotation_col = data.frame(Type = factor(rep(c("low", "high"), c(262,262))))
  
  rownames(annotation_col) = row.names(design)
  tiff("figure2-A-11.tiff",width = 4000,height = 2666,res = 300)
  pheatmap(map3,cellwidth = 1.5,scale = "row",show_rownames = F,show_colnames = F)
  dev.off()
  pdf("figure2-A-11.pdf",height = 10,width = 15)
  pheatmap(map3,cellwidth = 1.5,scale = "row",show_rownames = F,show_colnames = F)
  dev.off()
  
  
  write.csv(dif1,"tcga.dif3.csv")
  
  write.csv(tempOutput1,"tcga.temp3.csv")
  
  library("ggplot2")
  library("ggthemes")
  library(ggpubr)
  temp1<-tempOutput1
  temp1$log10pvalue <- -log10(temp1$P.Value)#??pvalue????ת??
  
temp1$ "not-signficant"
  temp1$group[which((temp1$P.Value < 0.05)&(temp1$logFC < -1))] = "down"
  temp1$group[which((temp1$P.Value < 0.05)&(temp1$logFC > 1))] ="up"
  
  table(temp1$group)
  
  temp1$labels = ""
  #temp1 <- temp1[order(temp1$P.Value),]
  upgenes <- head(row.names(temp1[temp1$group == 'up',]),30)
  downgenes <- head(row.names(temp1[temp1$group == 'down',]),30)
  temp1.30genes <- c(as.character(upgenes),as.character(downgenes))
  temp1$labels[match(temp1.30genes,row.names(temp1))] <- temp1.30genes
  
  tiff("figure112-a.tiff",width = 4000,height = 2450,res = 300)
  ggscatter(temp1,x="logFC",y="log10pvalue",
            color="group",
            palette=c("grey","red"),
            size = 1,
            label = temp1$labels,
            font.label = 8,
            repel = T,
            xlab = "log2FoldChange",
            ylab = "-log10pvalue",)+theme_base()+
    xlim(-3,3)+
    geom_hline(yintercept = 1.3,linetype="dashed")+geom_vline(xintercept = c(-1,1),linetype="dashed")
  dev.off()
  

  
  pdf("figure112-a.pdf",width = 15,height = 10)
  ggscatter(temp1,x="logFC",y="log10pvalue",
            color="group",
            palette=c("grey","red"),
            size = 1,
            label = temp1$labels,
            font.label = 8,
            repel = T,
            xlab = "log2FoldChange",
            ylab = "-log10pvalue",)+theme_base()+
    xlim(-3,3)+
    geom_hline(yintercept = 1.3,linetype="dashed")+geom_vline(xintercept = c(-1,1),linetype="dashed")
  dev.off()
  
  aa = data.frame(results)[intersect(row.names(tcga.KIRC.phenotype1),names(tcga.KIRC.RNA31)),]
  aa$risk=as.vector(ifelse(aa$Exhausted.CD8.T.cells>median(aa$Exhausted.CD8.T.cells),"high","low"))
  
  qqq = rbind(aa[which(aa$risk == "low"),],aa[which(aa$risk == "high"),])
  
  
  
  
  library(limma)
  
  group_list = c(rep("low",262), rep("high",262))
  
  design <- model.matrix(~0+factor(group_list))
  design
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(tcga.KIRC.RNA31[,row.names(qqq)])
  
  design
  
  contrast.matrix<-makeContrasts(high-low,levels = design)
  contrast.matrix ##??????????????????t <- lmFit(tcga.KIRC.RNA31[,row.names(qqq)],design)
  ##step2 ???ݶԱ?ģ?ͽ??в?ֵ.fit(fit, contrast.matrix) 
  ##step3 ??Ҷ˹????
  fit2 <t2) 
  ##step4 ???????л????ļ?????Table(fit2, coef=1, n=Inf)
  ##step5 ??P.Value????ɸѡ????ȫ????t2[tempOutput2[, "P.Value"]<0.05,]
  dif2 <- dif2[abs(dif2[,"logFC"])>1,]
  
  
  
  library(pheatmap)
  map3 <- tcga.KIRC.RNA31[row.names(dif2),row.names(qqq)]
  
  annotation_col = data.frame(Type = factor(rep(c("low", "high"), c(262,262))))
  
  rownames(annotation_col) = row.names(design)
  tiff("figure3-A-11.tiff",width = 4000,height = 2666,res = 300)
  pheatmap(map3,cellwidth = 1.5,scale = "row",show_colnames = F)
  dev.off()
  pdf("figure3-A-11.pdf",height = 10,width = 15)
  pheatmap(map3,cellwidth = 1.5,scale = "row",show_colnames = F)
  dev.off()
  
  
  write.csv(dif2,"tcga1.dif3.csv")
  
  write.csv(tempOutput2,"tcga1.temp3.csv")
  
  library("ggplot2")
  library("ggthemes")
  library(ggpubr)
  temp1<-tempOutput2
  temp1$log10pvalue <- -log10(temp1$P.Value)#??pvalue????ת??
  
temp1$ "not-signficant"
  temp1$group[which((temp1$P.Value < 0.05)&(temp1$logFC < -1))] = "down"
  temp1$group[which((temp1$P.Value < 0.05)&(temp1$logFC > 1))] ="up"
  
  table(temp1$group)
  
  temp1$labels = ""
  #temp1 <- temp1[order(temp1$P.Value),]
  upgenes <- head(row.names(temp1[temp1$group == 'up',]),30)
  downgenes <- head(row.names(temp1[temp1$group == 'down',]),30)
  temp1.30genes <- c(as.character(upgenes),as.character(downgenes))
  temp1$labels[match(temp1.30genes,row.names(temp1))] <- temp1.30genes
  
  tiff("figure1-a11.tiff",width = 4000,height = 2450,res = 300)
  ggscatter(temp1,x="logFC",y="log10pvalue",
            color="group",
            palette=c("blue","grey","red"),
            size = 1,
            label = temp1$labels,
            font.label = 8,
            repel = T,
            xlab = "log2FoldChange",
            ylab = "-log10pvalue",)+theme_base()+
    xlim(-2.5,2.5)+
    geom_hline(yintercept = 1.3,linetype="dashed")+geom_vline(xintercept = c(-1,1),linetype="dashed")
  dev.off()
  
  
  
  pdf("figure1-a11.pdf",width = 15,height = 10)
  ggscatter(temp1,x="logFC",y="log10pvalue",
            color="group",
            palette=c("blue","grey","red"),
            size = 1,
            label = temp1$labels,
            font.label = 8,
            repel = T,
            xlab = "log2FoldChange",
            ylab = "-log10pvalue",)+theme_base()+
    xlim(-2.5,2.5)+
    geom_hline(yintercept = 1.3,linetype="dashed")+geom_vline(xintercept = c(-1,1),linetype="dashed")
  dev.off()
  
  
  write.csv(sig.markers1,"sigmarkers1.csv")
  
  pdf(file="figure17.pdf",width=20,height=15)
  DoHeatmap(object = pbmc1,features = sig.markers11[,7]) 
  dev.off()
  tiff(file="figure17.tiff",width=4000,height=3000,res = 300)
  DoHeatmap(object = pbmc1,features =sig.markers11[,7] ) 
  dev.off()  
  
  
  pdf(file="hot1.pdf",width=20,height=15)
  DoHeatmap(object = pbmc1,features = c("IL7R","KLRG1","TBX21","CX3CR1")) 
  dev.off()
  tiff(file="hot1.tiff",width=4000,height=3000,res = 300)
  DoHeatmap(object = pbmc1,features =c("IL7R","KLRG1","TBX21","CX3CR1")) 
  dev.off() 
  
  
  pdf(file="hot2.pdf",width=20,height=15)
  DoHeatmap(object = pbmc1,features = c("TCF7","TOX","PDCD1")) 
  dev.off()
  tiff(file="hot2.tiff",width=4000,height=3000,res = 300)
  DoHeatmap(object = pbmc1,features =c("TCF7","TOX","PDCD1")) 
  dev.off() 
  
  pdf(file="hot3.pdf",width=20,height=15)
  DoHeatmap(object = pbmc1,features = c("IL7R","KLRG1","TBX21","CX3CR1","TCF7","TOX","PDCD1")) 
  dev.off()
  tiff(file="hot3.tiff",width=4000,height=3000,res = 300)
  DoHeatmap(object = pbmc1,features =c("IL7R","KLRG1","TBX21","CX3CR1","TCF7","TOX","PDCD1")) 
  dev.off() 
  
 string4 = read.table("string_interactions (54).tsv",header = F)
  unique() intersect(string4[,1],)
      
  aa = data.frame(t(tcga.KIRC.RNA31[as.character(a[which(a$color == "red"),1]),]) )  
  
 a = data.frame(a) 
  
  aa1 = aa[intersect(row.names(tcga.KIRC.phenotype1),row.names(aa)),]
  
  aa = aa[row.names(aa1),]

  tcga.KIRC.sur = data.frame(t(data.frame(t(tcga.KIRC.sur))))
  aa$OS = as.numeric(tcga.KIRC.sur[row.names(aa),1]) 
  aa$OS.time = as.numeric(tcga.KIRC.sur[row.names(aa),3]) 
  
  
  library(survival)
  library(survminer)
 r = c(1:79)
  
  for (i in 1:32) {
    aa$risk=as.vector(ifelse(aa[,i]>median(aa[,i]),"high","low"))
    
    aa1 = aa
    diff = survdiff(Surv(OS.time, OS) ~risk,data = aa1)
    pValue = 1-pchisq(diff$chisq,df=1)
    if(pValue<0.001){
      pValue=signif(pValue,4)
      pValue=format(pValue, scientific = TRUE)
    }else{
      pValue=round(pValue,3)
    }
    
    fit <- survfit(Surv(OS.time, OS) ~risk,data = aa1)
    ggsu =  ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
                       legend.labs=c("High", "Low"), legend.title=names(aa)[i],ggtheme = theme_bw(),
                       palette=c("red", "blue"))
    
    pdf(file=paste(names(aa)[i],".pdf",sep = ""),
        width=7,
        height=6,onefile = FALSE)
   print( ggsu)
    dev.off()
    tiff(file=paste(names(aa)[i],".tiff",sep = ""),
         width=2333,
         height=2000,res = 300)
   print( ggsu)
    dev.off()
    
  }
  
  
  
  
  library(survminer)
  
  pdf(file="figure19-a.pdf",
      width=7,
      height=6)
  ggsurvplot(fit, conf.int=TRUE, pval=TRUE, xlab="Days",ylab="Overall Survival",
             legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
             palette=c("red", "blue"))
  dev.off()
  tiff(file="figure19-a.tiff",
       width=2333,
       height=2000,res = 300)
  ggsurvplot(fit, conf.int=TRUE, pval=TRUE,xlab="Days",ylab="Overall Survival",
             legend.labs=c("High", "Low"), legend.title="Group",ggtheme = theme_bw(),
             palette=c("red", "blue"))
  dev.off()
  
  
  
  
  
 a =  data.frame(table(unlist(c(string4$V1,string4$V2)))) 

 a$color = a$Freq
 
 a$color[which(a$color >= 5)] = "red"
 a$color[which(a$color <= 4)] = "blue"
  links = string4[,1:2]
  nodes = a
  
  
  links = links[,c(1,3,2)]
  write.csv(links,"links0823.csv",quote = F)
  write.csv(nodes,"nodes0823.csv",quote = F)
  library(igraph) #????igraph??
  
  no= uniq(ao[,2])
  
  
  links = ao
  (net <- graph_from_data_frame(list(links),list(nodes),directed = F)) #????????????????
  s <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ('#0096ff', "lightblue", "azure3","firebrick1")
  V(net)$color<- V(net)$color  #???ýڵ?????ɫ
  Vt)$Freq + 2#???ñߵĿ??ȣ??????_from_data_frame(links)
  E(net)$width  <- E(net)$ao/5
  pdf("figureigra.pdf",height = 15,width = 15)
  
  set.seed(100)#?趨???ӽڵ㣬ͬһ?)
  
  plot(net,layout=layout.fruchterman.reingold,          #?ڵ??Ĳ???
       edge0.5,               #???ü?ͷ??С
       #v,
       # vertex.color=map(degree(net),c(1,20)),
       #vertex.size=1,
      # vertex.color="blue",
       vertex.frame.color="transparent",  #?ڵ??߿?͸??
       veror="black",        #?ڵ???ǩ??ɫ
        vex=1.2,             #?ڵ???ǩ??С
       # vamily="B",           #?ڵ???ǩ????
       edg                 #???Ƿ???????ȡֵ0-1??0k"  #?ߵ???ɫ
  )
  # legen=1.2, c("1","2","3","4","5","6","7","8","9","10","10-20",">=20"), 
  #        pch=21, col=unique(nodes[order(nodes$Freq),3]), pt.bg=unique(nodes[order(nodes$Freq),3]),pt.cex=3, cex=2)
  dev.off()
  
  tiff("figureigra.tiff",height = 3500,width = 3500,res = 300)
  set.seed(100)#?趨???ӽڵ㣬ͬһ?ֲ??  
  plot(net,layout=layout.fruchterman.reingold,          #?ڵ??Ĳ???
       edge.arro               #???ü?ͷ??С
       #vertex     # vertex.color=map(degree(net),c(1,20)),
       #vertex.size=1,
     #  vertex.color="blue",
       vertex.frame.color="transparent",  #?ڵ??߿?͸??
       vertex.black",        #?ڵ???ǩ??ɫ
        verte2,             #?ڵ???ǩ??С
       # vertey="B",           #?ڵ???ǩ????
       edge.c             #???Ƿ???????ȡֵ0-1??0Ϊ??#?ߵ???ɫ
  )
  # legend(x=, c("1","2","3","4","5","6","7","8","9","10","10-20",">=20"), 
  #        pch=21, col=unique(nodes[order(nodes$Freq),3]), pt.bg=unique(nodes[order(nodes$Freq),3]),pt.cex=3, cex=2)
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  aa = data.frame(t(tcga.KIRC.RNA31[c("AGER","BIRC5","CD69","CDC20","FOXM1","GAPDH","IL7R","SFTPD"),] ),class = rep(c("Normal","Tumor"),c(59,524))) 
  
  aa = data.frame(t(tcga.KIRC.RNA31[c("AGER","CD69","FOSB","GAPDH","IL7R","PIGR","SCGB3A1","SFTPD","SLC34A2"),] ),class = rep(c("Normal","Tumor"),c(59,524))) 
  library(ggpubr)
  rt = aa
  clinical="class" 
  rt = rt[order(rt[,clinical]),]
  
  c = c("a","b","c","d","e","f","g","h","i")
  
  group=levels(factor(rt$class))
  comp=combn(group,2)
  my_comparisons=list()
  for (i in 1:9) {
    for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
    #????boxplot
    boxplot=ggbot(rt, x="class", y=names(rt)[i], color="class",
                      xlab="",
                      ylab=paste(names(rt)[i]),
                      legend.title="",
                      #title=paste0("Cancer: ",i),
                      add = "jitter")+ 
      stat_compare_means(comparisons = my_comparisons)
    pdf(file=paste0("figure28","-",c[i],".pdf",sep = ""),width=7,height=5)
    print(boxplot)
    dev.off()
    tiff(file=paste0("figure28","-",c[i],".tiff",sep = ""),width=1750,height=1250,res = 300)
    print(boxplot)
    dev.off() 
  }
  
  aa2 = cbind(aa[row.names(tcga.KIRC.phenotype001),],tcga.KIRC.phenotype001)
    
  library(ggpubr)
  rt = aa2
  c = c(1:99)
  for (i in 11:14) {
    rt = rt[order(rt[,i]),]
    
    group=levels(factor(rt[,i]))
    comp=combn(group,2)
    my_comparisons=list()
    for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]} 
    
    for (n in 1:9) {
      #????boxplot
      boxplot=gglot(rt, x=names(rt)[i], y=names(rt)[n], color=names(rt)[i],
                        xlab="",
                        ylab=paste(names(rt)[n]),
                        legend.title="",
                        #title=paste0("Cancer: ",i),
                        add = "jitter")+ 
        stat_compare_means(comparisons = my_comparisons)
      pdf(file=paste0("figure","-",names(rt)[i],"-",names(rt)[n],".pdf",sep = ""),width=7,height=5)
      print(boxplot)
      dev.off()
      tiff(file=paste0("figure","-",names(rt)[i],"-",names(rt)[n],".tiff",sep = ""),width=1750,height=1250,res = 300)
      print(boxplot)
      dev.off() 
    }
   
    
  }
  
  
 
  
  write.table( cbind(Name = row.names(cd8gse), DESCRIPTION = rep("na",29634), cd8gse[,cluster21[order(cluster21$pbmc1.seurat_clusters),2]]),"cd81.txt",sep = "\t",row.names = F, quote = F)
  cluster21 = data.frame(clusters2) 
  write.csv(cluster21,"cluster21.csv")
  
  
     
  
  
 cbind(id=row.names(as.matrix(tcga.KIRC.RNA31[,row.names(aa1)])),as.matrix(tcga.KIRC.RNA31[,row.names(aa1)]))
  write.table(cbind(id=row.names(as.matrix(tcga.KIRC.RNA31[,row.names(aa1)])), as.matrix(tcga.KIRC.RNA31[,row.names(aa1)])),file="symbol.txt",quote=F,sep="\t",row.names=F)
  
  names(aa1)[8:9] = c("fustat","futime")
  aa1$risk=as.vector(ifelse(aa1$AGER>median(aa1$AGER),"high","low"))
 aa1$riskScore = aa1[,1]
  
  
 
  
  write.table(cbind(id=row.names(aa1),aa1[,8:11]),file="tcgaRisk.txt",quote=F,sep="\t",row.names=F)
 
  links = read.csv("lin.csv",header = T)
  nodes = read.csv("nodes.csv",header = T )
  
 links = links[,c(1,3,2)]
  library(igraph) #????igraph??
  
  nodes$color[which(nodes$color == "blue")] = "turquoise"
  (net <- graph_from_data_frame(list(links),list(nodes),directed = F)) #????????????????
  rainbow(#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ('#0096ff', "lightblue", "azure3","firebrick1")
  V(net)$color<- colrs[factor(V(net)$Freq)]  #???ýڵ?????ɫ
  V(net)$si+ 2#???ñߵĿ??ȣ???????????weita_frame(links)
  E(net)$width  <- E(net)$value/5
  pdf("figureigra.pdf",height = 10,width = 15)
 
   set.seed(1000)#?趨???ӽڵ㣬ͬһ?ֲ??ֻ??  plot(net,layout=layout.fruchterman.reingold,          #?ڵ??Ĳ???
       edge.arrow.si           #???ü?ͷ??С
       #vertex.la # vertex.color=map(degree(net),c(1,20)),
       #vertex.size=1,
       #vertex.color="blue",
       vertex.frame.color="transparent",  #?ڵ??߿?͸??
       vertex.labk",        #?ڵ???ǩ??ɫ
       # vertex.la
       # vertex.label.cex=1.2,             
       # vertex.label.family="B",                   #???Ƿ???????ȡֵ0-1??0Ϊ??????????ɫ
  )
 # legend(x=1.2, y1","2","3","4","5","6","7","8","9","10","10-20",">=20"), 
 #        pch=21, col=unique(nodes[order(nodes$Freq),3]), pt.bg=unique(nodes[order(nodes$Freq),3]),pt.cex=3, cex=2)
  dev.off()
  
  tiff("figureigra.tiff",height = 2000,width = 3000,res = 300)
  set.seed(1000)#?趨???ӽڵ㣬ͬһ?ֲ??ֻ???��?(net,layout=layout.fruchterman.reingold,          #?ڵ??Ĳ???
       edge.arrow.size=
       edge.arrow.size=0.5,              ertex.color=map(degree(net),c(1,20)),
       #vertex.size=1,
       #vertex.color="blue",
       vertex.frame.color="transparent",  #?ڵ??߿?͸??
       vertex.label.c       #?ڵ???ǩ??ɫ
       # vertex.label        #?ڵ???ǩ??С
       # vertex.label.         #?ڵ???ǩ????
       edge.curved=0,      #???Ƿ???????ȡֵ0-1??0Ϊ???????
   �
  )
  # legend(x=1.2, y="2","3","4","5","6","7","8","9","10","10-20",">=20"), 
  #        pch=21, col=unique(nodes[order(nodes$Freq),3]), pt.bg=unique(nodes[order(nodes$Freq),3]),pt.cex=3, cex=2)
  
  dev.off()
  
  
  
  library(plyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  setwd("C:\\Users\\dugu\\gsea_home\\output\\oct22\\my_analysis.Gsea.1634961400583\\cc")
  files=grep(".tsv",dir(),value=T)                                         #??ȡĿ¼?µ?????xls?ļ?
  data = la                                        #??ȡÿ???ļ?
  names(data) = files= ldply(data, data.frame)
  dataSet$pathway = gsub(".tsv","",dataSet$.id)                            #???ļ???׺ɾ??
  
  gseaCol=c("#58C"#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
    geom_point(shape=21) + scale_fill_manual(values = gseaCol[1:nrow(dataSet)]) + 
    labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0),limits =c(min(dataSet$RUNNING.ES-0.02), max(dataSet$RUNNING.ES+0.02))) +   
    theme_bw() + theme(panel.grid =element_blank()) + theme(panel.border = element_blank()) + 
    theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
    geom_hline(yintercept = 0) + guides(fill=guide_legend(title = NULL)) + 
    theme(legend.background = element_blank()) + theme(legend.key = element_blank())
  pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
    scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
    labs(x = "High <----------------->Low ", y = "", title = "") + 
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)
  
  gGsea = ggplot_gtable(ggplot_build(pGsea))
  gGene = ggplot_gtable(ggplot_build(pGene))
  maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
  gGsea$widths = as.list(maxWidth)
  gGene$widths = as.list(maxWidth)
  dev.off()
  
  #??ͼ?ο??ӻ?????????"multipleGSEA.pdA.pdf',      #???ͼƬ???ļ?
      width=9,     ????ͼƬ?߶?
      height=5)  
      height=5)               e(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
  dev.off()
  
  tiff('multipleGSEA.tiff',      #????ͼƬ???ļ?
       width=2700, ????????ͼƬ?߶?
       height=1500         #????????ͼƬ?߶?
  par(mar=c(5,5,2e(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
  dev.off()
  dataSet$pathway = gsub("KEGG_","",dataSet$pathway)                            #???ļ???׺ɾ??
  
  write.csv(data")
  
  
  
  
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(dplyr)
  

  data <- read.csv('loc_cran_packages.csv')
  
  # most popular programming languages from TIOBE Index (Nov. 2019) found in data
  # (only languages with position <= 16 are considered)
  popular_languages <- unique(data$language)
  
  # number of packages to display
  number_of_pkgs <- length(unique(data$pkg_name))
  
  # find largest packages written in popular languages
  top_packages <- data %>%
    filter(language %in% popular_languages) %>%
    group_by(pkg_name) %>% summarize(total_code = sum(code)) %>%
    arrange(desc(total_code)) %>%
    head(number_of_pkgs) %>%
    select(pkg_name, total_code)
  
  # all popular languages per package
  top_languages_per_pkg <- data %>%
    filter(
      pkg_name %in% top_packages$pkg_name,
      language %in% popular_languages
    ) %>%
    arrange(pkg_name, desc(code)) %>%
    group_by(pkg_name) %>%
    mutate(
      main = row_number() == 1, # main language of package should be opaque
      total_code = sum(code)
    ) %>%
    ungroup() %>%
    select(language, pkg_name, code, total_code, main)
  
  # only following languages found in given packages
  (top_languages <- top_languages_per_pkg %>%
      pull(language) %>%
      sort)
  
  top_language_colors <- c(
    '#efb306',
    '#cd023d',
    '#0f8096', '#000000'
    
  )
  
  names(top_language_colors) <- unique(data$language)
  
  edges1 <- top_languages_per_pkg %>%
    transmute(from = language, to = pkg_name, total_code = code, main)
  
  edges2 <- top_languages_per_pkg %>%
    count(language, wt = code, name = 'total_code') %>%
    transmute(
      from = '',
      to = language,
      total_code,
      main = TRUE
    )
  
  edges <- bind_rows(edges1, edges2)
  
  vertices1 <- top_languages_per_pkg %>%
    filter(main) %>%
    transmute(
      node = pkg_name, language, total_code, level = 1
    )
  
  vertices2 <- edges2 %>%
    transmute(
      node = to, language = to, total_code, level = 2
    )
  
  vertices3 <- tibble(
    node = '', language = NA, total_code = 0, level = 3
  )
  
  vertices <- bind_rows(vertices1, vertices2, vertices3) %>%
    mutate(
      radius = total_code**(1.8), # scaling circles
      language = factor(language, names(top_language_colors))
    ) %>%
    arrange(level, language, node)
  
  graph <- graph_from_data_frame(edges, vertices = vertices)
  
  # create custom layout by updating existing circle layout
  layout <- create_layout(graph, layout = 'circle')
  
  outer_circle <- layout %>%
    filter(level == 1) %>%
    mutate(language = factor(language, names(top_language_colors))) %>%
    arrange(language, desc(name)) %>%
    mutate(
      x = cos((row_number() - 1.5) / number_of_pkgs * 2 * pi),
      y = sin((row_number() - 1.5) / number_of_pkgs * 2 * pi)
    )
  
  # positioning circle centers manually by specifying polar coords
  angles <- c(10,  180,   290,190)
  radii <- c(0.5, 0.3, 0.35,0.7)
  centers <- tibble(
    x = radii * cos(angles / 180 * pi),
    y = radii * sin(angles / 180 * pi)
  )
  inner_circle <- bind_cols(centers, select(filter(layout, level != 1), -x, -y))
  
  
  layout[] <- bind_rows(outer_circle, inner_circle) %>% arrange(.ggraph.index)
  
  
  
  #layout[87:89,7] = c(800,1000,600)
  
  
  tiff("cfdd.tiff",width = 9500,height = 9500,res = 300)
  
  ggraph(layout) +
    geom_edge_diagonal(
      aes(edge_color = node1.language, edge_alpha = as.factor(main)),
      edge_width = 2, show.legend = FALSE
    ) +
    geom_node_point(
      aes(size = (100*radius+6), color = language),
      alpha = 0.5, show.legend = FALSE
    ) + geom_node_text(
      aes(
        x = 1.4 * x,
        y = 1.4 * y,
        label = "",
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        filter = !(name %in% top_languages)
      ),
      size = 5, hjust = 'outward', family = 'Oswald'
    )+
    geom_node_text(
      aes(
        x = 1.0175 * x,
        y = 1.0175 * y,
        label = name,
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        filter = !(name %in% top_languages)
      ),
      size = 6, hjust = 'outward', family = 'Oswald'
    ) +
    geom_node_text(
      aes(
        x =x,
        y = y,
        label = name,
        filter = name %in% top_languages
      ),
      size = 8, hjust = 0.5, family = 'Oswald'
    ) +
    # geom_node_text(
    #  aes( x = x, y = y - 0.045,label = ifelse(total_code > 1000,format(total_code, big.mark = ','), total_code),filter = name %in% top_languages
    # ), size = 7, hjust = 0.5, family = 'Oswald') +
    scale_edge_color_manual(values = top_language_colors) +
    scale_color_manual(values = top_language_colors) +
    scale_size_area(max_size = 150) +
    scale_edge_alpha_manual(values = c(0.15, 1)) +
    coord_fixed()+
    # labs( title = 'ccc', subtitle = '',caption = '') +
    theme_void() 
  
  #+
  #theme(text = element_text(family = 'Oswald'), legend.position = c(0.645, 0.51), plot.title = element_text(
  #  face = 'bold', hjust = 0.5, size = 20, margin = margin(t = 45, b = 3) ), plot.subtitle = element_text(
  #  face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),plot.caption = element_text(
  #  face = 'plain', color = '#dedede', size = 8, hjust = 1, margin = margin(b = 20)))
  dev.off()
  
  
  
  
  
  #layout[87:89,7] = c(800,1000,600)
  pdf("cfdd.pdf",width = 35,height = 35)
  
  ggraph(layout) +
    geom_edge_diagonal(
      aes(edge_color = node1.language, edge_alpha = as.factor(main)),
      edge_width = 2, show.legend = FALSE
    ) +
    geom_node_point(
      aes(size = (100*radius+6), color = language),
      alpha = 0.5, show.legend = FALSE
    ) + geom_node_text(
      aes(
        x = 1.4 * x,
        y = 1.4 * y,
        label = "",
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        filter = !(name %in% top_languages)
      ),
      size = 5, hjust = 'outward'#, family = 'Oswald'
    )+
    geom_node_text(
      aes(
        x = 1.0175 * x,
        y = 1.0175 * y,
        label = name,
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        filter = !(name %in% top_languages)
      ),
      size = 6, hjust = 'outward'#, family = 'Oswald'
    ) +
    geom_node_text(
      aes(
        x =x,
        y = y,
        label = name,
        filter = name %in% top_languages
      ),
      size = 8, hjust = 0.5#, family = 'Oswald'
    ) +
    # geom_node_text(
    #  aes( x = x, y = y - 0.045,label = ifelse(total_code > 1000,format(total_code, big.mark = ','), total_code),filter = name %in% top_languages
    # ), size = 7, hjust = 0.5, family = 'Oswald') +
    scale_edge_color_manual(values = top_language_colors) +
    scale_color_manual(values = top_language_colors) +
    scale_size_area(max_size = 150) +
    scale_edge_alpha_manual(values = c(0.15, 1)) +
    coord_fixed()+
    # labs( title = 'ccc', subtitle = '',caption = '') +
    theme_void() 
  
  #+
  #theme(text = element_text(family = 'Oswald'), legend.position = c(0.645, 0.51), plot.title = element_text(
  #  face = 'bold', hjust = 0.5, size = 20, margin = margin(t = 45, b = 3) ), plot.subtitle = element_text(
  #  face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),plot.caption = element_text(
  #  face = 'plain', color = '#dedede', size = 8, hjust = 1, margin = margin(b = 20)))
  dev.off()
  
  
  pancan = read.csv("pancan.csv",header = T,row.names = 1)
  #pancan = log2(pancan)
  pdf("figure271111slc34a2.pdf",width = 12,height = 8)
  par(las = 2)
  barplot(as.matrix(t(pancan)),beside = T,ylim = c(0,10),col = c("red","blue"),legend=names(pancan),ylab = "log2(TPM+1)")
  dev.off()
  tiff("figure271111slc34a2.tiff",width = 2100,height = 1400,res = 300)
  par(las = 2)
  barplot(as.matrix(t(pancan)),beside = T,ylim = c(0,10),col = c("red","blue"),legend=names(pancan),ylab = "log2(TPM+1)")
  dev.off()
  
  
  
  gse43458 = read.csv("GSE43458_series_matrix.csv",header = T,row.names = 1)
  
  gse43458.phe = read.csv("gse4345.csv",header = T,row.names = 2)
  
  
  GPL6244 = read.csv("GPL6244-17930.csv",header = T,row.names = 1)
  
  grep("IL7R",GPL6244$gene_assignment,value = T)
  
  grep("AGER",GPL6244$gene_assignment,value = T)
  
  grep("GAPDH",GPL6244$gene_assignment,value = T)
  grep("CD69",GPL6244$gene_assignment,value = T)
  
  
  table(gse43458.phe$X.Sample_characteristics_ch1)
  row.names(GPL6244)[grep("IL7R",GPL6244$gene_assignment)]
  
  row.names(GPL6244)[grep("AGER",GPL6244$gene_assignment)]
  
  row.names(GPL6244)[grep("GAPDH",GPL6244$gene_assignment)]
  
  row.names(GPL6244)[grep("CD69",GPL6244$gene_assignment)]
  
  row.names(gse43458) = paste0("tt",row.names(gse43458))
  
  
  
  mean(as.numeric(gse43458["tt7953385",]),as.numeric(gse43458["tt8104901",]),as.numeric(gse43458["tt8104901",]))
  
  
 kj = data.frame(CD69 = as.numeric(gse43458["tt7961075",]),class = rep(c("LUAD(never-smoker)","Normal","LUAD(smoker)"),c(40,30,40)))
 
 
 
 
 
 
 library(ggpubr)
 
 group=levels(factor(kj$class))
 comp=combn(group,2)
 my_comparisons=list()
 for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
 #????boxplot

   boxplot=ggboxplot(kj"class", y = names(kj)[1],color="class",fill = "class",
                     xlab="",
                     ylab=paste(names(kj)[1]),
                     legend.title="",
                     title=paste0(""),cor.coef.size = 10#ͼƬ??pֵ?ȱ?ǩ?Ĵ?С
             r"
   )+ stat_compare_means(comparisons = my_comparisons)+theme(text = element_text(size = 10))#xlab??ylab?Լ???С
   #+stat_comparf(file=paste0(names(kj)[1],".pdf"),width=5,height=5)
   print(boxplot)
   dev.off()
   tiff(file=paste0(names(kj)[1],".tiff"),width=1500,height=1500,res = 300)
   print(boxplot)
   dev.off()
 
   
   kj = data.frame(IL7R = as.numeric(gse43458["tt8104901",]),class = rep(c("LUAD(never-smoker)","Normal","LUAD(smoker)"),c(40,30,40)))
   
   
   
   
   
   
   library(ggpubr)
   
   group=levels(factor(kj$class))
   comp=combn(group,2)
   my_comparisons=list()
   for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
   #????boxplot
   
   boxplot=ggboxplot x="class", y = names(kj)[1],color="class",fill = "class",
                     xlab="",
                     ylab=paste(names(kj)[1]),
                     legend.title="",
                     title=paste0(""),cor.coef.size = 10#ͼƬ??pֵ?ȱ?ǩ?Ĵ?С
             r"
   )+ stat_compare_means(comparisons = my_comparisons)+theme(text = element_text(size = 10))#xlab??ylab?Լ???С
   #+stat_comparf(file=paste0(names(kj)[1],".pdf"),width=5,height=5)
   print(boxplot)
   dev.off()
   tiff(file=paste0(names(kj)[1],".tiff"),width=1500,height=1500,res = 300)
   print(boxplot)
   dev.off()
   
   
   
   
   
   
   
   kj = data.frame(GAPDH = as.numeric(gse43458["tt7953385",]),class = rep(c("LUAD(never-smoker)","Normal","LUAD(smoker)"),c(40,30,40)))
   
   
   
   
   
   
   library(ggpubr)
   
   group=levels(factor(kj$class))
   comp=combn(group,2)
   my_comparisons=list()
   for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
   #????boxplot
   
   boxplot=ggboxplot x="class", y = names(kj)[1],color="class",fill = "class",
                     xlab="",
                     ylab=paste(names(kj)[1]),
                     legend.title="",
                     title=paste0(""),cor.coef.size = 10#ͼƬ??pֵ?ȱ?ǩ?Ĵ?С
             r"
   )+ stat_compare_means(comparisons = my_comparisons)+theme(text = element_text(size = 10))#xlab??ylab?Լ???С
   #+stat_comparf(file=paste0(names(kj)[1],".pdf"),width=5,height=5)
   print(boxplot)
   dev.off()
   tiff(file=paste0(names(kj)[1],".tiff"),width=1500,height=1500,res = 300)
   print(boxplot)
   dev.off()
   
   
   
   row.names(GPL6244)[grep("AGER",GPL6244$gene_assignment)]
   
   
   
   
   
   
   
   kj = data.frame(AGER = (as.numeric(gse43458["tt8125341",])+as.numeric(gse43458["tt8178771",])+as.numeric(gse43458["tt8179967",]))/3,
                   class = rep(c("LUAD(never-smoker)","Normal","LUAD(smoker)"),c(40,30,40)))
   
   
   
   
   
   
   library(ggpubr)
   
   group=levels(factor(kj$class))
   comp=combn(group,2)
   my_comparisons=list()
   for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
   #????boxplot
   
   boxplot=ggboxplot x="class", y = names(kj)[1],color="class",fill = "class",
                     xlab="",
                     ylab=paste(names(kj)[1]),
                     legend.title="",
                     title=paste0(""),cor.coef.size = 10#ͼƬ??pֵ?ȱ?ǩ?Ĵ?С
             r"
   )+ stat_compare_means(comparisons = my_comparisons)+theme(text = element_text(size = 10))#xlab??ylab?Լ???С
   #+stat_comparf(file=paste0(names(kj)[1],".pdf"),width=5,height=5)
   print(boxplot)
   dev.off()
   tiff(file=paste0(names(kj)[1],".tiff"),width=1500,height=1500,res = 300)
   print(boxplot)
   dev.off()
   
   