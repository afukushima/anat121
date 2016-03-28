#############################################################################
## R-script for reconstruction of osteoblast-specific coexpression networks 
##
## [for ANAT 2016, http://www.procomu.jp/anat121/]
##
## Please note that you should have all dependent packages.
##
## written by Atsushi Fukushima [atsushi.fukushima AT riken.jp]
#############################################################################
library(GEOquery)
library(affy)
library(mouse4302.db)
library(igraph)

##############################################
## pre-processing function (downloads raw data and normalizes/summarizes them by RMA)
##
## Please note that the following process is time-consuming.
## Alternatively, you can use pre-processed data matrix (filename: mydata.txt).
## in R prompt,
##  res <- read.table("mydata.txt", header=TRUE, row.names=1, sep="\t")
##  dim(res) ## 11277  34 
## Then, please go to "constructing correlation networks".
pre.processing <- function(GEO.id) {
    data <- getGEOSuppFiles(GEO.id)
    tarfile <- paste(GEO.id, "/", GEO.id, "_RAW.tar", sep="")
    untar(tarfile, exdir=GEO.id)
    ## Data pre-processing
    ## target files
    tgt <- list.files(GEO.id, pattern="*.CEL.gz|*.cel.gz", full.names=TRUE, )
    ## RMA normalization (Bolstad et al., Bioinformatics 2003)
    eset <- justRMA(filenames=tgt)
    ## centering
    sc <- scale(exprs(eset), scale=FALSE)
    return(sc)
}

## target files
files <- read.table("target_GSE.txt", header=FALSE, sep="\t")
files <- as.character(t(files))

res <- matrix(nrow=45101)
i <- 1
for(GEO.id in files) {
    tmp <- pre.processing(GEO.id)
    print(paste(i, GEO.id, sep=":"))
    res <- cbind(res, tmp)
    i <- i+1
}
res <- res[,-1]

dim(res)  ## 45101 314   (36 datasets, see "target_GSE.txt")


###############################################
## Filtering cross-hybridization probesets and control probes
###############################################
probesets <- rownames(res)
rmv1 <- grep("s_at$", probesets)
rmv2 <- grep("x_at$", probesets)
rmv3 <- grep("^AFF", probesets)
rmv <- c(rmv1, rmv2, rmv3)
res2 <- res[-rmv,]
dim(res2) ## 40508 314


###############################################
## filtering with variance
###############################################
var.thr <- 0.3
vars <- apply(res2, 1, var)
res3 <- res2[vars >= var.thr,]
dim(res3) ## 27467  314


###############################################
## filtering with mean values
###############################################
library(genefilter)
k.samp <- round(dim(res3)[2] * 0.2)    ## sample number
sig.A <- 2**-3                         ## signal intensity A
f1 <- kOverA(k.samp, A=sig.A)          ## expression measure above A in at least k samples
ffun <- filterfun(f1)                  
res4 <- res3[genefilter(res3, ffun),]  ## filter out
dim(res4) ## 18191  314


###############################################
## Summarization of probesets and annotations
###############################################
x <- mouse4302SYMBOL
# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- unlist(as.list(x[mapped_probes]))

symbols <- unlist(as.list(x[rownames(res4)]))

## unification of annotation
unique.sym <- unique(symbols)        
unique.sym <- unique.sym[unique.sym != ""]
unique.sym <- unique.sym[!is.na(unique.sym)]
unique.sym <- unique.sym[!is.nan(unique.sym)]

mydata <- t(apply(as.matrix(unique.sym), 1,
                  function(i, d = data.frame(res4), s, p) {
                      apply(d[which(s == i), ], 2, p, na.rm = TRUE)
                  }, data.frame(res4), symbols, mean)
            )

rownames(mydata) <- unique.sym 
tmp <- cbind(rownames(mydata), mydata)
write.table(tmp, file="mydata.txt", quote=FALSE, row.names=FALSE, sep="\t")

dim(mydata)  ## 11277  314

###############################################
## constructing correlation networks
###############################################

## mydata
mydata.cor <- cor(t(mydata), method="spearman")
mydata.cor <- round(mydata.cor, 1)
dim(mydata.cor)  ## 11277 11277

cor.thr <- 0.7
g <- simplify( graph.adjacency( (mydata.cor >= cor.thr), mode="undirected") )
summary(g)
## output
outfile <- paste("CorrNet_sp_thr", cor.thr, ".txt", sep="")
write.table(get.edgelist(g),
            file=outfile,
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE,
            sep="\t")

## You should import this output file by Cytoscape (http://www.cytoscape.org/).
