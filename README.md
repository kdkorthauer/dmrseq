# dmrseq: Inference for differentially methylated regions (DMRs) from bisulfite sequencing

An central question in the analysis of bisulfite sequencing data 
is to detect regions (collections of 
neighboring CpGs) with systematic differences between conditions, 
as compared to within-condition variability. These so-called *Differentially
Methylated Regions* (DMRs) are thought to be more informative than single CpGs 
in terms of of biological function. The package dmrseq 
provides a rigorous permutation-based approach to
detect and perform inference for differential methylation by use of 
generalized least squares models that account for inter-individual and 
inter-CpG variability to generate region-level statistics that can be
comparable across the genome. The framework performs well even
on samples as small as two per group. The vignette explains the
use of the package and demonstrates typical workflows.
