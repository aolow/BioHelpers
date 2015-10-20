## GOstats code to do Hypergeometric Testing
## adapted from vignettes towards my oncogene addiction dataset

# load list of genes of interest
dat <- read.csv("~/Analysis/Interaction_model_Addiction_p=0.01.csv")

require("GOstats")
library("lumiHumanAll.db")
library(annotate)

# load list of all genes in the microarray
load(file="./eset.RData")
all_genes <- unique(as.character(getEG(row.names(exprs(eset)), "lumiHumanAll.db")))

# p-value threshold for significance
p <-0.0001

dat$geneID <- as.character(getEG(as.character(unlist(dat$X)), "lumiHumanAll.db"))


# set parameters
params <- new("GOHyperGParams",
              geneIds = as.character(dat$GeneID[1:200]), # my gene list 
              universeGeneIds = all_genes, # total gene list
              pvalueCutoff = p,
              annotation = "org.Hs.eg.db", # expression array platform
              ontology = "BP", # GO category
              conditional = TRUE, # take into account GO categ. dependencies
              testDirection = "over") # test overrepresentation

# calculate enrichment
res <- hyperGTest(params)
pval <- -log(summary(res)[,"Pvalue"])
terms <- summary(res)[,"Term"]
names(pval) <- terms

# quick plot
library(lattice)
par(oma=c(15,1,1,1))
barplot(pval[1:30], las=2, col="blue", ylim=c(0,16), ylab="-10log(P-Value))")

#htmlReport(res, file="")

####### KEGG
p <- 0.05
# set parameters
params <- new("KEGGHyperGParams",
              geneIds = as.character(dat$GeneID[1:200]), # my gene list 
              universeGeneIds = all_genes, # total gene list
              pvalueCutoff = p,
              annotation = "org.Hs.eg.db", # expression array platform
              testDirection = "over") # test overrepresentation

# calculate enrichment
res <- hyperGTest(params)

# investigate results
names(summary(res))
summary(res)[,c("Pvalue","Count","Size","Term")]

#htmlReport(res, file="KEGG_Culture_p05.html")

pval <- -log(summary(res)[,"Pvalue"])
terms <- summary(res)[,"Term"]
names(pval) <- terms

# quick plot
par(oma=c(15,1,1,1))
barplot(pval, las=2, col="red", ylim=c(0,15), ylab="-10log(P-Value))", main="KEGG Pathway Enrichment")
