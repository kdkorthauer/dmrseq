## The following script uses the dmrseq package to
## constructs the annot.chr21 dataset,
## an example of annotatr annotation to use in plotting
## included in the dmrseq package for example purposes

library(dmrseq)

# load example data 
data(BS.chr21)

# get annotation information for hg19
annot.chr21 <- getAnnot("hg19")

# only keep this information for chromosome 21 (for example dataset)
annot.chr21[[1]] <- annot.chr21[[1]][seqnames(annot.chr21[[1]])=="chr21",]
annot.chr21[[2]] <- annot.chr21[[2]][seqnames(annot.chr21[[2]])=="chr21",]

save(annot.chr21 , file = "./data/annot.chr21.rda")
library(tools)
resaveRdaFiles("./data/annot.chr21.rda")
