## The following script uses the dmrseq package to
## constructs the dmrs.ex dataset,
## an example of dmrseq output
## included in the dmrseq package for example purposes

library(dmrseq)

# load example data 
data(BS.chr21)

# the covariate of interest is the "CellType" column of pData(BS.chr21)
testCovariate <- "CellType"

# run dmrseq on a subset of the chromosome (20K CpGs)
dmrs.ex <- dmrseq(bs=BS.chr21[240001:260000,],
                  cutoff = 0.05,
                  testCovariate=testCovariate)

save(dmrs.ex, file = "./data/dmrs.ex.rda")
library(tools)
resaveRdaFiles("./data/dmrs.ex.rda")
