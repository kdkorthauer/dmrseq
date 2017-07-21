## The following script downloads and constructs the BS.chr21 dataset,
## included in the dmrseq package
## This script was modified from bsseq - includes two additional samples
## compared to the example data in bsseq

library(dmrseq)

## First we download.  Each file is slightly less than 200 MB
# you may need to change the 'method' in 'download.file' to suit
# the utilities available on your OS

# cell type imr90, replicates 1 and 2
download.file(url = "ftp://ftpuser3:s3qu3nc3@neomorph.salk.edu/mc/mc_imr90_r1.tar.gz",
              destfile = "mc_imr90_r1.tar.gz", method='curl')
untar("mc_imr90_r1.tar.gz",  "mc_imr90_r1/mc_imr90_r1_21", compressed = TRUE)
download.file(url = "ftp://ftpuser3:s3qu3nc3@neomorph.salk.edu/mc/mc_imr90_r2.tar.gz",
              destfile = "mc_imr90_r2.tar.gz", method='curl')
untar("mc_imr90_r2.tar.gz",  "mc_imr90_r2/mc_imr90_r2_21", compressed = TRUE)

# cell type h1, replicates 1 and 2
download.file(url = "ftp://ftpuser3:s3qu3nc3@neomorph.salk.edu/mc/mc_h1_r1.tar.gz",
              destfile = "mc_h1_r1.tar.gz", method='curl')
untar("mc_h1_r1.tar.gz",  "mc_h1_r1/mc_h1_r1_21", compressed = TRUE)
download.file(url = "ftp://ftpuser3:s3qu3nc3@neomorph.salk.edu/mc/mc_h1_r2.tar.gz",
              destfile = "mc_h1_r2.tar.gz", method='curl')
untar("mc_h1_r2.tar.gz",  "mc_h1_r2/mc_h1_r2_21", compressed = TRUE)


## Now the workhorse function

read.lister <- function(file) {
  dat <- read.table(file, skip = 1, row.names = NULL,
                    col.names = c("chr", "pos", "strand", "context", "M", "Cov"),
                    colClasses = c("character", "integer", "character",
                                   "character", "integer", "integer"))
  ## we remove all non-CpG calls.  This includes SNPs
  dat <- dat[dat$context == "CG",]
  dat$context <- NULL
  dat$chr <- paste("chr", dat$chr, sep = "")
  ## Now we need to handle that the data has separate lines for each strand
  ## We join these
  tmp <- dat[dat$strand == "+",]
  BS.forward <- BSseq(pos = tmp$pos, chr = tmp$chr, 
                      M = as.matrix(tmp$M, ncol = 1),
                      Cov = as.matrix(tmp$Cov, ncol = 1), 
                      sampleNames = "forward")
  tmp <- dat[dat$strand == "-",]
  BS.reverse <- BSseq(pos = tmp$pos - 1L, chr = tmp$chr, 
                      M = as.matrix(tmp$M, ncol = 1),
                      Cov = as.matrix(tmp$Cov, ncol = 1), 
                      sampleNames = "reverse")
  BS <- combine(BS.forward, BS.reverse)
  BS <- collapseBSseq(BS, columns = c("a", "a"))
  BS
}

BS.imr90.r1 <- read.lister("mc_imr90_r1/mc_imr90_r1_21")
sampleNames(BS.imr90.r1) <- "imr90.r1"
BS.imr90.r2 <- read.lister("mc_imr90_r2/mc_imr90_r2_21")
sampleNames(BS.imr90.r2) <- "imr90.r2"

BS.h1.r1 <- read.lister("mc_h1_r1/mc_h1_r1_21")
sampleNames(BS.h1.r1) <- "h1.r1"
BS.h1.r2 <- read.lister("mc_h1_r2/mc_h1_r2_21")
sampleNames(BS.h1.r2) <- "h1.r2"

BS.chr21 <- combine(BS.imr90.r1, BS.imr90.r2, BS.h1.r1, BS.h1.r2)
pData(BS.chr21)$CellType <- c(rep("imr90",2), rep("h1",2))
pData(BS.chr21)$Rep <- rep(c("replicate1", "replicate2"),2)
validObject(BS.chr21)
pData(BS.chr21)

#remove any loci which have no coverage in one or more samples
meth.levels.raw = getMeth(BS.chr21, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
BS.chr21 = BS.chr21[-no.hits]

save(BS.chr21, file = "./data/BS.chr21.rda")
library(tools)
resaveRdaFiles("./data/BS.chr21.rda")


