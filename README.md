# dmrseq: Inference for differentially methylated regions (DMRs) from bisulfite sequencing

A central question in the analysis of bisulfite sequencing data 
is to detect regions (collections of 
neighboring CpGs) with systematic differences between conditions, 
as compared to within-condition variability. These so-called *Differentially
Methylated Regions* (DMRs) are thought to be more informative than single CpGs 
in terms of of biological function. 

<p align="center">
  <img src="/inst/sticker/dmrseq.png" height="300"/>
</p>

The package **dmrseq** 
provides a rigorous permutation-based approach to
detect and perform inference for differential methylation by use of 
generalized least squares models that account for inter-individual and 
inter-CpG variability to generate region-level statistics that can be
comparable across the genome. The framework performs well even
on samples as small as two per group. 

## Installation

**dmrseq** is available on 
[Bioconductor](https://bioconductor.org/packages/dmrseq). You can install
it with R version 3.5.0 or higher with the following commands:

```
source("https://bioconductor.org/biocLite.R")
biocLite("dmrseq")
```

## Getting started

See the vignette for information on how to use the package to perform
typical methylation analysis workflows.

## Learn more

More details of the **dmrseq** framework can be found in the manuscript

> Korthauer, K., Chakraborty, S., Benjamini, Y., and Irizarry, R.A.
> Detection and accurate False Discovery Rate control of differentially 
methylated regions from Whole Genome Bisulfite Sequencing
> *Biostatistics*, 2018 (in press).
> [BioRxiv:10.1101/183210](http://www.biorxiv.org/content/early/2017/08/31/183210)


## License/Copyright
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 
This package is made available under an MIT license.  
