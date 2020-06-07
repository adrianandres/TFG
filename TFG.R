# SCRIPT GSE6731

# ADRIAN ANDRES LECICA // TFG


##############################
# INSTALACIONES
##############################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("oligo")
BiocManager::install("limma")
BiocManager::install("annotate")
BiocManager::install("hgu95av2")
BiocManager::install("pd.hg.u95av2")
BiocManager::install("KEGG.db")
BiocManager::install("GSEABase")
BiocManager::install("clusterProfiler")

install.packages("openxlsx")
install.packages("scatterplot3d")
install.packages("gplots")

library(GEOquery)
library(oligo) # Analiza las muestras cargadas desde el ordenador

directorio<-"C:/Users/Universidad/Desktop/TFG"
setwd(directorio)

##############################
# DESCARGA DATOS DESDE GEO
##############################

gseid <-"GSE6731" # Descargamos los datos gracias al geoquery, paquete para bajar datos de geo (datos comprimidos)
gse <- getGEO(gseid,GSEMatrix=TRUE) # ExpressionSet
GEO <- getGEOSuppFiles(gseid,makeDirectory=TRUE) # Descargamos datos crudos
setwd(paste("C:/Users/Universidad/Desktop/TFG",gseid,sep="/")) # Guardamos los datos en el escritorio
untar("GSE6731_RAW.tar", exdir = getwd()) # Descomprimimos los datos (rar)
list.files() # Lista de archivos
celFiles = list.files( pattern = "CEL.gz") # CEL o EXP
GEOFS <- read.celfiles(celFiles) # Cargamos los datos crudos en R

##############################
# ANOTACIÓN FENOTÍPICA
##############################

GEOFS # ExpressionFeatureSet
sample(GEOFS)
head(exprs(GEOFS)) 
phenoData(GEOFS)
pData(GEOFS) # PhenoData
pD <- pData(gse[[1]])
pD

##############################
# CONTROL DE CALIDAD
##############################

image(GEOFS) # Imágenes de los arrays

# Evaluacion de los controles de calidad
colos<-rainbow(6) # Escala de colores
hist(GEOFS,target="core",col=colos) 
boxplot(GEOFS,target="core",col=colos,las=3,cex.axis=0.5)
MAplot(GEOFS) 
MAplot(GEOFS,pairs=TRUE) 

##############################
# NORMALIZACIÓN                    
##############################

#GEOFS.rma<-rma(GEOFS,target="core")
load("GEOFS.rma.RData")
dim(exprs(GEOFS.rma))
GEOFS.rma
class(GEOFS.rma)
pData(GEOFS.rma)=pData(gse[[1]])

# Graficos después de la normalización
hist(GEOFS.rma,target="core",col=colos) 
boxplot(GEOFS.rma,target="core", col=colos,las=3,cex.axis=0.5)
MAplot(GEOFS.rma) 
MAplot(GEOFS.rma,pairs=TRUE) 

##############################
# DISTRIBUCIÓN MUESTRAL
##############################

x<-exprs(GEOFS.rma)
names=as.character(GEOFS.rma$title)  # Establecemos nombres muestras
sampleNames(GEOFS.rma)


# Agrupacion jerarquica
clust.cor.ward<- hclust(as.dist(1-cor(x)),method="ward.D2") # Coeficiente de correlación de Pearson y método de linkage (Ward.D2)
plot(clust.cor.ward, labels=names,main="hierarchical clustering", hang=-1,cex=0.6)

library(scatterplot3d)

# Creamos vector de condiciones
cond=gsub("[0-9]","",names)
cond=gsub("--","_",cond)
cond=gsub("s|-","",cond)
cond=as.factor(cond)

cols=c("blue","lightblue","yellow","yellow","orange","green","black","grey")[cond] # Colores ordenados por condición
summary(pca.filt <- prcomp(t(x), cor = TRUE ))

# Principal component analysis
pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3], color=cols, xlab='PC1', ylab='PC2', zlab='PC3', main='PCA', pch=16,col.grid="lightblue")
text(pca3d$xyz.convert(pca.filt$x+0.5), labels=names, cex=0.5,pos = 2,
     offset = 1)
legend("topright",      
       legend = levels(cond),
       col = c("blue","lightblue","yellow","yellow","orange","green","black","grey"), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)

# Obtención % variabilidad explicado por cada componente
vars <- apply(pca.filt$x, 2, var)  
props <- round(vars / sum(vars),2)

# PCA bidimensional
m <- as.data.frame(pca.filt$x[, 1:3])
m$cols=cols
m$names=names
m$cond=cond
library(ggplot2)
ggplot(data = m, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(cond))) + 
  geom_text(aes(label=names),hjust=0, vjust=0, size=3)+
  scale_color_manual(values=c("blue","lightblue","yellow","yellow","orange","green","black","grey"))+
  theme_bw()

ggplot(data = m, aes(x = PC2, y = PC3)) +
  geom_point(aes(color = factor(cond))) + 
  geom_text(aes(label=names),hjust=0, vjust=0, size=3)+
  scale_color_manual(values=c("blue","lightblue","yellow","yellow","orange","green","black","grey"))+
  theme_bw()

##############################
# ANALISIS DE EXPRESIÓN DIFERENCIAL
##############################

library(limma)

# Diseño de la matriz para las condiciones
design<-model.matrix(~0+cond)
design
colnames(design)<-gsub("cond","",colnames(design))
rownames(design)<-sampleNames(GEOFS)
design

# Comparación de las condiciones seleccionadas
fit<-lmFit(GEOFS.rma,design)
contrast.matrix<-makeContrasts(CD_Aff-N,levels=design) # Condiciones a comparar (CD_Aff-N, UC_Aff-N y CD_Aff-UC_Aff)
fit2<-contrasts.fit(fit,contrast.matrix)
fite<-eBayes(fit2) # Metodo empírico de Bayes para moderar los errores estándar de los cambios estimados de log-fold.
top.table<-topTable(fite,coef=1,number=Inf,adjust="BH") # Método de corrección del p-valor aplicado (Benjamini & Hochberg)
results<-decideTests(fite)
table(results) 

#Distribución del p-valor
hist(top.table$P.Value,breaks=100,main="results P")

# Resultados con p-valor ajustado < 0,05
results.p0.05<-top.table[top.table$adj.P.Val<0.05,]
dim(results.p0.05)

# Resultados logFC>1
results.p0.05.logFC1<-top.table[top.table$adj.P.Val<0.05 & abs(top.table$logFC)>1,]
dim(results.p0.05.logFC1)

results

##############################
# ANOTACIÓN GENETICA
##############################

library(annotate) 
library(hgu95av2.db)

# Resultados
logFC <- results.p0.05.logFC1$logFC
pval <- results.p0.05.logFC1$P.Value
adj.pval<-results.p0.05.logFC1$adj.P.Val

# Plataforma usada para la anotación
entrezid<-mget(rownames(results.p0.05.logFC1), env=hgu95av2ENTREZID) # SIMBOLO
entrezid

dat <- exprs(GEOFS.rma)[rownames(results.p0.05.logFC1),]
sym<-mget(rownames(results.p0.05.logFC1), env=hgu95av2SYMBOL) # SIMBOLO
name<-mget(rownames(results.p0.05.logFC1), env=hgu95av2GENENAME) # NOMBRE
chr<-mget(rownames(results.p0.05.logFC1), env=hgu95av2CHR) # POSICI?N

affyids<-rownames(results.p0.05.logFC1)
genelist <- list(affyids) # ENTREZ ID
genelist
filename <- "Resultados.html"
title <- "GENES DIFERENCIALMENTE EXPRESADOS"
othernames <- list(sym,name,chr,round(logFC, 1), round(pval, 4), round(adj.pval, 4), round(dat, 2)) 
head <- c("Sonda", "Simbolo", "Nombre gen", "Posici?n cromos?mica","logFC", "p-valor","p-value ajustado",names)
repository <- list("affy")

htm<-htmlpage(genelist, filename, title, othernames, head, repository = repository) # Fichero en el directorio
df<-data.frame(htm)
write.table(df, file="GENESDIFERENCIALMENTEEXPRESADOS.txt")

# Heatmap solo de biopsias afectadas
library(gplots)

data.clus<-exprs(GEOFS.rma[rownames(results.p0.05.logFC1),])
rownames(data.clus)<-sym
colnames(data.clus)<-names

# Seleccionamos las condiciones para el heatmap
#i=grep("^N|CD-[0-9]+-Aff",colnames(data.clus))
#i=grep("^N|UC-[0-9]+-Aff",colnames(data.clus))
#i=grep("^UC.+-Aff|^CD.+-Aff",colnames(data.clus))
#i=grep("^UC.+Aff|^CD.+Aff|^CD.+Aff",colnames(data.clus))

data.clus=data.clus[,i]

clust.rows <- hclust(as.dist(1-cor(t(data.clus))),method="ward.D2")
clust.cols <- hclust(as.dist(1-cor(data.clus)),method="ward.D2")  
heatcol<-colorRampPalette(c("green", "Black","red"), space = "rgb") # Establecer escala de color
cols1=cols[i]

pdf("heatmap.result.pdf")
heatm<-heatmap.2(as.matrix(data.clus[ ,]), col = heatcol(256), # Utiliza la matriz de dadas normalizadas. 
                 dendrogram="both", Colv=as.dendrogram(clust.cols),
                 Rowv=as.dendrogram(clust.rows), # Escalamos por filas
                 scale="row",labRow = NULL, cexCol=0.8, 
                 main="",key=TRUE,keysize=1,
                 density.info="none",trace="none", 
                 ColSideColors = cols1 , sepcolor = TRUE)
legend("topright",      
       legend =unique(as.character(cond[i])),
       col = unique(cols1),
       lty= 1,             
       lwd = 5,           
       cex=.7
)
dev.off()

#### Heatmap biopsias afectadas + biopsias no afectadas
data.clus<-exprs(GEOFS.rma[rownames(results.p0.05.logFC1),])
rownames(data.clus)<-sym
colnames(data.clus)<-names

#i=grep("^N|^CD",colnames(data.clus))
#i=grep("^N|^UC",colnames(data.clus))
#i=grep("^UC|^CD",colnames(data.clus))
data.clus=data.clus[,i]

clust.rows <- hclust(as.dist(1-cor(t(data.clus))),method="ward.D2")
clust.cols <- hclust(as.dist(1-cor(data.clus)),method="ward.D2")  
heatcol<-colorRampPalette(c("green", "Black","red"), space = "rgb") # Establecer escala de color
cols1=cols[i]

pdf("heatmap.UN.result.pdf")
heatm<-heatmap.2(as.matrix(data.clus[ ,]), col = heatcol(256), # Utiliza la matriz de dadas normalizadas. Pasamos del negro al blaue en 256 tramos
                 dendrogram="both", Colv=as.dendrogram(clust.cols),
                 Rowv=as.dendrogram(clust.rows), # Escalamos por filas
                 scale="row",labRow = NULL, cexCol=0.8, 
                 main="",key=TRUE,keysize=1,
                 density.info="none",trace="none", 
                 ColSideColors = cols1 , sepcolor = TRUE)
legend("topright",      
       legend =unique(as.character(cond[i])),
       col = unique(cols1),
       lty= 1,             
       lwd = 5,           
       cex=.7
)
dev.off()

library(clusterProfiler)
library(openxlsx)

# Base de datos KEEG
kegg <- read.gmt( "c2.cp.kegg.v7.1.symbols.gmt")

##############################
# KEGG ENRIHMENT ANALYSIS
##############################

# Creamos objeto para guardar el archivo excel
wb <- createWorkbook()

# Asignamos los genes up-regulados y los genes down-regulados 
symUP<-mget(rownames(results.p0.05.logFC1[results.p0.05.logFC1$logFC>0,]), env=hgu95av2SYMBOL) # SIMBOLO
symDOWN<-mget(rownames(results.p0.05.logFC1[results.p0.05.logFC1$logFC<0,]), env=hgu95av2SYMBOL) # SIMBOLO

pdf("KEGG.enrichment.pdf")

# Enriquecimiento
enrchmt <- enricher(symUP, TERM2GENE=kegg)
dotplot(enrchmt, showCategory=30,font.size=8,title=paste("KEGG enrichment of",length(symUP),"genes"))
addWorksheet(wb, "KEGG UP") 
writeData(wb, "KEGG UP", enrchmt@result,  rowNames = F) # Añadimos los datos a la hoja excell
enrchmt <- enricher(symDOWN, TERM2GENE=kegg)
dotplot(enrchmt, showCategory=30,font.size=8,title=paste("KEEG enrichment of",length(symDOWN),"genes"))
addWorksheet(wb, "KEGG DOWN")
writeData(wb, "KEGG DOWN", enrchmt@result,  rowNames = F)
dev.off()

# Guardamos enrichment
saveWorkbook(wb, file = "ENRICHMENT_KEEG.xlsx", overwrite = TRUE)

##############################
# GENE SET ENRICHMENT ANALYSIS
##############################

geneList=sign(top.table$logFC)*(-log10(top.table$P.Value)) # Estadistico utilizado
symGSEA<-mget(rownames(top.table), env=hgu95av2SYMBOL) # Símbolo GSEA
names(geneList)=symGSEA 
geneList = sort(geneList ,decreasing = TRUE)
gseaobject<-GSEA(
  geneList,
  exponent = 1,
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE=kegg,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)

View(gseaobject@result)

# Guardamos GSEA en excel
gseawb <- createWorkbook()
addWorksheet(gseawb, "RUTAS METABOLICAS")
writeData(gseawb, "RUTAS METABOLICAS", gseaobject@result,  rowNames = F) 
saveWorkbook(gseawb, file = "ENRICHMENT_GSEA.xlsx", overwrite = TRUE)

