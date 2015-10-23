#### Heatmap code with legend bars from expression data sets

#setwd()

## Load libraries
library("gplots")
library("devtools")
library("Biobase")
library("RColorBrewer")
library("limma")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


### Functions for color mapping

color.map_oa <- function(oa) { if (oa==0) "grey" else "black" }
color.map_er <- function(er) { if (er==0) "grey" else "black" }
color.map_erb <- function(erb) { if (erb==0) "grey" else "black" }
color.map_edmfs <- function(edmfs) { if (edmfs ==0) "grey" else "black" }


## Load and clean data; make into an expression dataset

load(file="dat_compact.RData")

dat_compact2 <- na.omit(dat_compact)

OA_pred <- dat_compact2$OncogeneAddiction_Pred
er_stat <- dat_compact2$er
erb_stat <- dat_compact2$erb
edmfs_stat <- dat_compact2$e.dmfs

# make eset

pData <- cbind(OA_pred,er_stat,erb_stat, edmfs_stat)
rownames(pData) <- dat_compact2$tumorID
pData <- as.data.frame(pData)

metadata <- data.frame(labelDescription=
                         c("Oncogene Addiction",
                           "ER status",
                           "ERB status",
                           "Events distant metastasis free survival"),
                       row.names=c("OA_pred", "er_stat", "erb_stat", "edmfs_stat"))

phenoData <- new("AnnotatedDataFrame",
                 data=pData, varMetadata=metadata)

experimentData <- new("MIAME",
                      name = "Ola",
                      lab = "Van 't Veer Lab",
                      title = "Oncogene Addiction Project")

eset <- ExpressionSet(assayData=t(temp3),
                      phenoData=phenoData,
                      experimentData=experimentData)

### make legend bar
addict_colors <- unlist(lapply(eset$OA_pred, color.map_oa))
er_colors <- unlist(lapply(eset$er_stat, color.map_er))
erb_colors <- unlist(lapply(eset$erb, color.map_erb))
edmfs_colors <- unlist(lapply(eset$edmfs, color.map_edmfs))

clab <- cbind(edmfs_colors, er_colors, addict_colors, erb_colors)
colnames(clab)=c("e.dmfs", "er","addiction","erb")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidean")}
myclust=function(c) {hclust(c,method="ward.D2")}


# extract expression data from eset
y <- exprs(eset)

# plot heatmap
main_title="EMC344"
par(cex.main=1)


heatmap.3(y,
          hclustfun=myclust, 
          distfun=mydist, 
          scale="row",
          na.rm = TRUE, 
          dendrogram="column", 
          margins=c(8, 15.5),
          Rowv=TRUE, Colv=TRUE, 
          ColSideColors=clab, 
          symbreaks=FALSE, key=TRUE, 
          keysize=0.8,
          symkey=FALSE,
          density.info="none", trace="none", main=main_title, 
          #labCol=samples, 
          #labRow=lbl, 
          cexRow=0.9, col=rev(redgreen(75)),
          ColSideColorsSize=4
          #   KeyValueName="Prob. Response"
)
dev.off()

# add legend
legend("topright",
       legend=c(unique(samples),"",
                "DMSO", "LAP","",
                "Non-addicted", "Addicted", "",
                "Parent-Cell-derived","Cell-derived", "Mouse-derived"),
       fill=c(brewer.pal(8, "Set1"),
              "white",
              "lightblue","orange", "white",
              "royalblue","orangered", "white",
              "lightgrey","darkgrey","black"),
       border=FALSE, bty="n", 
       #y.intersp = 0.7, 
       cex=0.6)

#dev.off()

