library(BoolNet)


#A)


redjoo <- loadNetwork("01_RawData/redjoo.txt")
redjoo
atrajoo <- getAttractors(redjoo)

plotAttractors(atrajoo)

plotStateGraph(atrajoo, drawLabels = T)


#######Discusion

########Grafica 2


immu <- loadNetwork("01_RawData/immuno.txt") #Agregar un renglon mas para el txt
immu
atraimm <- getAttractors(immu)

plotAttractors(atraimm)




###########STARWARS#############

library(igraph)

starwars <- read.csv("01_RawData/star-wars-network-edges.csv")
View(starwars)




starwarsmat <- as.matrix(starwars)

redstarwars <- graph.data.frame(starwarsmat, directed = T)

plot(redstarwars)

##########Distribucion de conectivisades
distlegostar <- degree.distribution(redstarwars)
plot(distlegostar)

####diametros de la red

get.diameter(redstarwars)
diameter(redstarwars)


###clusterizar con los pesos
gruposstr <- cluster_edge_betweenness(redstarwars, weights = NULL, modularity = T, edge.betweenness = F)

plot(gruposstr, redstarwars)



###visualizar la red donde se observa una propiedad topologica

pesos <- E(redstarwars)$weight
pesos <- as.numeric(pesos)


E(redstarwars)$width <- log(pesos) + 1

plot(redstarwars)
###############modificar la red para que se vea mas bonito

######modificar estilo de la red tanto de nodos como de conectividades






















############WGCNA##############
BiocManager::install("GO.db")
BiocManager::install("impute")
library(GO.db)
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("preprocessCore")
library(WGCNA)

datosmodif <- read.table("WGCNA_tutorial-main/gene_counts_table_WGCNA_LC.txt",header=T,row.names=1,sep="\t")
#ron.names=1 pine en automatico los nombres y ya no hay que modificar de colnames <- rownames
#creo

View(datosmodif)

datasample <- read.csv("WGCNA_tutorial-main/sample_info.csv")
View(datasample)

dataExpr <- DESeqDataSetFromMatrix(countData = datosmodif[,-181],colData = datasample,design = ~ Zone) #modifica los datos para
#que puedna ajustarse a los otros datos?
mcols(dataExpr)$basepairs <- datosmodif$geneLengt1
fpkm_matrix <- fpkm(dataExpr) #hace el ajuste/normalizacion con log(FPKM+1)
datExpr_1 <- t(log2(fpkm_matrix+1))


head(datExpr_1[1:5,1:5]) # samples in row, genes in column
match(datasample$sample_ID, colnames(datosmodif))
datExpr_1 <- datExpr_1[,1:5000]


# Calculate sample distance and cluster the samples
sampleTree <- hclust(dist(datExpr_1), method = "average");
# plot sample tree
pdf(file = "05_Images/1-n-sampleClustering.pdf", width = 40, height = 9);
par(cex = 1.3);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

# Choose a set of soft threshold parameters
powers <- c(c(1:20), seq(from = 22, to=30, by=2))

powers


######

sft <- pickSoftThreshold(datExpr_1, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#####NO IMPRIME LA GRAFICA, es la grafica para ver el free scale?
#la grafica de la linea 78



###########Turn data expression into topological overlap matrix############
#opcion A

# Turn data expression into topological overlap matrix
power=sft$powerEstimate 

# Option 1: automatic
cor <- WGCNA::cor
net = blockwiseModules(datExpr_1, power = power,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)

cor<- stats::cor
# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file = "05_Images/4-module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


######OPCION 2
# Option 2a: step-by-step
power = power
adjacency = adjacency(datExpr_1, power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM


# Option 2b: ESTE ES EL IMPORTANTEEEEEE O CON EL SIGUE EL TUTOTIAL
TOM = TOMsimilarityFromExpr(datExpr_1, power = power) ######difeentes metodos para sacar la similaridad
dissTOM = 1-TOM 
dim(dissTOM)




#####################  Construct modules (proceed with the genetree from option 2b)#########


# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");

pdf(file = "05_Images//3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()  ######NO PLOTEA, SOLO SI LO PASA A pdf

# Module identification using dynamic tree cut
# We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "05_Images//4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()


##########MERGE MODULES JUNTAR LOS QUE QUEREMOS

# Calculate eigengenes
MEList = moduleEigengenes(datExpr_1, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.40
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr_1, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "05_Images//5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
write.table(merge$oldMEs,file="oldMEs.txt");
write.table(merge$newMEs,file="newMEs.txt");
#juntatr los que s ejuntaron en una tabla?

############Export of networks to external software##################


# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
  modules = c(substring(names(merge$oldMEs)[i], 3));
  genes = colnames(datExpr_1)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/orign_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}
# Export the gene list of new modules 
for (i in 1:length(merge$newMEs)){
  modules = c(substring(names(merge$newMEs)[i], 3));
  genes = colnames(datExpr_1)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}


#############PART 1: Correlate module eigen-genes and samples (or other discrete data)###############
# Heatmap of old module eigen-genes and samples
#pdf(file="oldMEs.pdf",heigh=80,width=20)
library("pheatmap")
rownames(merge$oldMEs)=names(data0[,-181])
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=F,show_colnames=T,fontsize=6)
dev.off()


# Heatmap of new module eigen-genes and sample trait (e.g. Zone)
col_ann <- datasample[,c(1,3)]
rownames(col_ann) <- col_ann[,1]
col_ann <- data.frame(col_ann)
col_ann$Zone <- as.factor(col_ann$Zone)
col_ann <- col_ann[order(col_ann$Zone),]
col_ann$sample_ID <- NULL
head(col_ann)
ann_color <- list("col_ann" = c("Z1" = "yellow",
                                "Z2" = "red",
                                "Z3" = "green"))

data <- data.frame(merge$newMEs)
data <- data[order(match(rownames(data), rownames(col_ann))),]
dim(merge$newMEs)

install.packages("pheatmap")
library(pheatmap)

pdf(file="05_Images//newMEs.pdf",heigh=60,width=20)
rownames(merge$newMEs)=names(datosmodif[,-181])
pheatmap(data,cluster_col=T,cluster_row=F,show_rownames=F,
         show_colnames=T,fontsize=6,
         annotation_row = col_ann, annotation_colors = ann_color)
dev.off()



#=====================================================================================
#
#############PART 2: Correlation between gene modules and microbial traits (continuous data)############

# Define numbers of genes and samples
nGenes = ncol(datExpr_1);
nSamples = nrow(datExpr_1);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_1, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# Read microbial data as traits
bac_traits = read.table("Yu2021NaturePlants-master/traits_file/b_order_234.txt", header = T, sep = "\t")
rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
####
write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");


#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("05_Images//5-module-traits-bacteria-order.pdf", width = 80, height = 15)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(bac_traits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor[1:25,1:25], 2), "\n(",
                    signif(moduleTraitPvalue[1:25,1:25], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[1:25,1:25])
pdf("05_Images//5-module-traits-bacteria-order1.pdf", width = 20, height = 10)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor[1:25,1:25],
               xLabels = colnames(bac_traits[1:25,1:25]),
               yLabels = colnames(MEs[1:25,1:25]),
               ySymbols = colnames(MEs[1:25,1:25]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
###########Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignficance###########

# Define variable Verru containing the Verrucomicrobiales column of bac_traits
Verru = as.data.frame(bac_traits$Verrucomicrobiales);
names(Verru) = "Verrucomicrobiales"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = orderMEs(cbind(MEs, Verru))

geneModuleMembership = as.data.frame(cor(datExpr_1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr_1, Verru, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Verru), sep="");
names(GSPvalue) = paste("p.GS.", names(Verru), sep="");

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


module = "lightgreen"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Verrucomicrobiales",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")


## Draw bubble plot for particular module
colsum_bac_traits <- colSums(bac_traits)
colsum_bac_traits <- data.frame(colsum_bac_traits)
colsum_bac_traits$b_order <- rownames(colsum_bac_traits)
library(tidyr)
moduleTraitCor_long <- data.frame(moduleTraitCor)
moduleTraitCor_long$module <- rownames(moduleTraitCor)
moduleTraitCor_long <- moduleTraitCor_long[,c(235,1:234)]
moduleTraitCor_long <- gather(moduleTraitCor_long, b_order, PCC, Pseudomonadales:Others, factor_key = TRUE)

moduleTraitPvalue_long <- data.frame(moduleTraitPvalue)
moduleTraitPvalue_long$module <- rownames(moduleTraitPvalue)
moduleTraitPvalue_long <- moduleTraitPvalue_long[,c(235,1:234)]
moduleTraitPvalue_long <- gather(moduleTraitPvalue_long, b_order, pval, Pseudomonadales:Others, factor_key = TRUE)

moduleTrait_long <- merge(moduleTraitCor_long, moduleTraitPvalue_long, by = c("module","b_order"))

bubble_Data <- merge(moduleTrait_long, colsum_bac_traits, by = "b_order")
#just want module = "lightgreen"
bubble_Data_lightgreen <- bubble_Data[which(bubble_Data$module == "MElightgreen"),]

library(ggplot2)
ggplot(bubble_Data_lightgreen, aes(x= colsum_bac_traits, y= PCC, size = colsum_bac_traits,
                                   color = PCC, label = b_order)) +
  geom_text(hjust = 1, size=3) +
  geom_point(alpha=1) + ylab("Module-taxon correlation") + xlab("Relative abundance (sum)") +
  theme_bw()


############# Summary ###################################

head(datExpr_1)[1:5,1:5] # transcriptome data

head(datasample)[1:5,] # metadata (sample info)
head(bac_traits)[1:5,1:5] # external trait







