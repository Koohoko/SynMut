## SynMut: Tools for Designing Synonymously Mutated Sequences with Different Genomic Signatures

### Intro

*SynMut* designs synonymous mutants for DNA sequences. 
There are increasing demands on designing virus mutants with specific dinucleotide or codon composition. This tool can take both dinucleotide preference and/or codon usage bias into account while designing mutants. It also works well for desinging mutants with extremely over-/under- represented dinucleotides. 
This tool was originally designed for generating recombinant virus sequences in influenza A virus to study the effects of different dinucleotide usage and codon usage, yet the functions provided in this package can be generic to a variety of other biological researches.

### Components of the package

![image](https://raw.githubusercontent.com/Koohoko/Koohoko.github.io/master/SynMut/images/component.png)

### Installation 
Use the below code to install the development version of this package.

```r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

if (!requireNamespace("SynMut", quietly = TRUE))
    devtools::install_github("Koohoko/SynMut")

library(SynMut)
```

### Example and methods

Details tutorial please refer to the [vignette](https://koohoko.github.io/SynMut/index.html).

The strategies and functionalities of the `codom_mimic` and `dinu_to` functions can be found at [here](https://koohoko.github.io/SynMut/algorithm.html).

