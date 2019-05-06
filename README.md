## SynMut: Tools for Designing Synonymously Mutated Sequences with Different Genomic Signatures
*Haogao Gu, Leo L.M. Poon*

 <img src="https://raw.githubusercontent.com/Koohoko/Koohoko.github.io/master/SynMut/images/sph_logo.png" alt="drawing" width="200" ALIGN="LEFT" /> 
 
 ##### This work was conducted in School of Public Health, The University of Hong Kong under the supervison of Prof. Leo Poon.

***
[DOI: 10.18129/B9.bioc.SynMut](https://doi.org/doi:10.18129/B9.bioc.SynMut)  
<img border="0" src="https://bioconductor.org/shields/availability/3.9/SynMut.svg">
<img border="0" src="https://bioconductor.org/shields/build/devel/bioc/SynMut.svg">

### Introduction

*SynMut* designs synonymous mutants for DNA sequences. 

There are increasing demands on designing virus mutants with specific dinucleotide or codon composition. This tool can take both dinucleotide preference and/or codon usage bias into account while designing mutants. It also works well for desinging mutants with extremely over-/under- represented dinucleotides. 

This tool was originally designed for generating recombinant virus sequences in influenza A virus to study the effects of different dinucleotide usage and codon usage, yet the functions provided in this package can be generic to a variety of other biological researches.

### Components of the package

![image](https://raw.githubusercontent.com/Koohoko/Koohoko.github.io/master/SynMut/images/component.png)

### Installation 
Use the below code to install the package.

```r
# Stable version
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

if (!requireNamespace("SynMut"))
    BiocManager::install("SynMut")

# Development version
if (!requireNamespace("devtools"))
    install.packages("devtools")

if (!requireNamespace("SynMut"))
    devtools::install_github("Koohoko/SynMut")
```

### Example and methods

Details tutorial please refer to the [vignette](https://koohoko.github.io/SynMut/index.html).

The strategies and functionalities of the `codom_mimic` and `dinu_to` functions can be found at [here](https://koohoko.github.io/SynMut/algorithm.html).

### Find it at Bioconductor
https://bioconductor.org/packages/devel/bioc/html/SynMut.html
***
