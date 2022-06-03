# Class definition --------------------------------------------------------

#' An S4 class to record DNA sequences and variable regions for mutations
#'
#' Recording codon DNA sequences and region.
#'
#' @slot dnaseq a DNAStingSet object recording the sequence(s)
#' @slot region a list specifying paticular regions in the sequences
#'   allowed to be mutated
#'
#' @importFrom methods setClass
#' @import Biostrings
#' @seealso \code{\link{input_seq}}, \code{\link{get_cu}},
#'   \code{\link{get_region}}
#' @author Haogao Gu
#' @name regioned_dna-class
#' @rdname regioned_dna-class
#' @exportClass regioned_dna

setClass(Class = "regioned_dna",
    slots = c(dnaseq = "DNAStringSet",
        region = "list"))

setValidity("regioned_dna",
    function(object) {
        object.region <- object@region
        object.seq <- object@dnaseq
        check.na <- all(is.na(object.region))
        dnaseq.n <- length(object.seq)
        region.n <- length(object.region)
        check.3 <-
            all(vapply(object.seq, length, numeric(1)) %% 3 == 0)
        if (!check.na) {
            check.length <- dnaseq.n == region.n
            check.num <-
                vapply(object.seq, length, numeric(1)) / 3 ==
                vapply(object.region, length, numeric(1))
        } else {
            check.num <- TRUE
            check.length <- TRUE
        }
        if (!check.3) {
            stop("the length of the dna sequences must be divisible by 3")
        }
        if (all(check.num, check.length)) {
            TRUE
        } else{
            "dnaseq and region must have the same length."
        }
    })

##viewer
setMethod(
    f = "show",
    signature = "regioned_dna",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat("Number of sequences: ", length(object@dnaseq) - 1, "\n", sep = "")
    }
)


# Accessor ----------------------------------------------------------------

#' Get the variable region
#'
#' Access the variable regions
#'
#' @param object regioned_dna
#' @param ... ...
#'
#' @return list
#' @seealso \code{\link{input_seq}}, \code{\link{get_cu}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_region(rgd.seq)
#' @name get_region
#' @rdname get_region-methods
#' @exportMethod get_region
setGeneric(
    name = "get_region",
    def = function(object, ...) standardGeneric("get_region")
)

#' @rdname get_region-methods
setMethod(
    f = "get_region",
    signature = "regioned_dna",
    definition = function(object) {
        object@region[seq_len(length(object@region) - 1)]
    }
)

#' Get the DNAStringSet data
#'
#' Access the DNA sequence data in DNAStringSet.
#'
#' @param object A regioned_dna object.
#' @param ... ...
#'
#' @return DNAStringSet
#' @exportMethod get_dna
#' @name get_dna
#' @rdname get_dna-methods
#'
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_dna(rgd.seq)
setGeneric(
    name = "get_dna",
    def = function(object, ...) standardGeneric("get_dna")
)

#' @rdname get_dna-methods
setMethod(
    f = "get_dna",
    signature = "regioned_dna",
    definition = function(object) {
        return(object@dnaseq[seq_len(length(object@dnaseq) - 1)])
    }
)


# simple calculation ------------------------------------------------------

#' Get codon usage matrix
#'
#' Access the codon usage matrix
#'
#' @param object regioned_dna / DNAStringSet
#' @param ... ...
#'
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_nu}}, \code{\link{get_du}}, \code{\link{get_freq}},
#'   \code{\link{get_rscu}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_cu(rgd.seq)
#' @exportMethod get_cu
#' @name get_cu
#' @rdname get_cu-methods
setGeneric(
    name = "get_cu",
    def = function(object, ...) standardGeneric("get_cu")
)

#' @rdname get_cu-methods
setMethod(
    f = "get_cu",
    signature = "regioned_dna",
    definition = function(object) {
        dnaseq <- object@dnaseq[seq_len(length(object@dnaseq) - 1)]
        oligonucleotideFrequency(dnaseq, width = 3, step = 3)
    }
)

#' @rdname get_cu-methods
setMethod(
    f = "get_cu",
    signature = "DNAStringSet",
    definition = function(object) {
        oligonucleotideFrequency(object, width = 3, step = 3)
    }
)

#' Get dinucleotide usage matrix
#'
#' Access the dinucleotide usage matrix
#'
#' @param object regioned_dna / DNAStringSet
#' @param ... ...
#'
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_nu}}, \code{\link{get_cu}}, \code{\link{get_freq}},
#'   \code{\link{get_rscu}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_du(rgd.seq)
#' @exportMethod get_du
#' @name get_du
#' @rdname get_du-methods
setGeneric(
    name = "get_du",
    def = function(object, ...) standardGeneric("get_du")
)

#' @rdname get_du-methods
setMethod(
    f = "get_du",
    signature = "regioned_dna",
    definition = function(object) {
        dnaseq <- object@dnaseq[seq_len(length(object@dnaseq) - 1)]
        oligonucleotideFrequency(dnaseq, width = 2, step = 1)
    }
)

#' @rdname get_du-methods
setMethod(
    f = "get_du",
    signature = "DNAStringSet",
    definition = function(object) {
        oligonucleotideFrequency(object, width = 2, step = 1)
    }
)

#' Get nucleotide usage matrix
#'
#' Access the nucleotide usage matrix
#'
#' @param object regioned_dna / DNAStringSet
#' @param ... ...
#'
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_cu}}, \code{\link{get_du}}, \code{\link{get_rscu}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_nu(rgd.seq)
#' @exportMethod get_nu
#' @name get_nu
#' @rdname get_nu-methods
setGeneric(
    name = "get_nu",
    def = function(object, ...) standardGeneric("get_nu")
)

#' @rdname get_nu-methods
setMethod(
    f = "get_nu",
    signature = "regioned_dna",
    definition = function(object) {
        dnaseq <- object@dnaseq[seq_len(length(object@dnaseq) - 1)]
        oligonucleotideFrequency(dnaseq, width = 1, step = 1)
    }
)

#' @rdname get_nu-methods
setMethod(
    f = "get_nu",
    signature = "DNAStringSet",
    definition = function(object) {
        oligonucleotideFrequency(object, width = 1, step = 1)
    }
)

#' Get codon usage frequency of synonymous codons
#'
#' Access the synonymous codon usage frequency
#'
#' @param object regioned_dna / DNAStringSet / codon usage matrix (vector)
#' @param numcode The ncbi genetic code number for translation. Default value: \code{1}. Details please refer to \code{?seqinr::translate} ("https://rdrr.io/cran/seqinr/man/translate.html").
#' @param ... ...
#'
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_cu}}, \code{\link{get_du}}, \code{\link{get_rscu}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_freq(rgd.seq)
#' @exportMethod get_freq
#' @name get_freq
#' @rdname get_freq-methods
setGeneric(
    name = "get_freq",
    def = function(object, numcode = 1, ...) standardGeneric("get_freq")
)

#' @rdname get_freq-methods
setMethod(
    f = "get_freq",
    signature = "regioned_dna",
    definition = function(object, numcode) {
        dnaseq <- object@dnaseq[seq_len(length(object@dnaseq) - 1)]
        tmp <- oligonucleotideFrequency(dnaseq,
            width = 3, step = 3)
        tmp <- apply(tmp, 1, function(x) {
            r <- freq(x, numcode)
            names(r) <- NULL
            r <- unlist(r)[names(x)]
            return(r)
        })
        return(t(tmp))
    }
)

#' @rdname get_freq-methods
setMethod(
    f = "get_freq",
    signature = "DNAStringSet",
    definition = function(object, numcode) {
        tmp <-
            oligonucleotideFrequency(object, width = 3, step = 3)
        tmp <- apply(tmp, 1, function(x) {
            r <- freq(x, numcode)
            names(r) <- NULL
            r <- unlist(r)[names(x)]
            return(r)
        })
        return(t(tmp))
    }
)

#' @rdname get_freq-methods
setMethod(
    f = "get_freq",
    signature = "matrix",
    definition = function(object, numcode) {
        tmp <- apply(object, 1, function(x) {
            r <- freq(x, numcode)
            names(r) <- NULL
            r <- unlist(r)[names(x)]
            return(r)
        })
        return(t(tmp))
    }
)

#' @rdname get_freq-methods
setMethod(
    f = "get_freq",
    signature = "vector",
    definition = function(object, numcode) {
        r <- freq(object, numcode)
        names(r) <- NULL
        r <- unlist(r)[names(object)]
        return(r)
    }
)

#' Get Relative Synonymous Codon Usage (rscu) of synonymous codons
#'
#' Access the Relative Synonymous Codon Usage rscu
#'
#' @param object regioned_dna / DNAStringSet / codon usage matrix (vector)
#' @param numcode The ncbi genetic code number for translation. Default value: \code{1}. Details please refer to \code{?seqinr::translate} ("https://rdrr.io/cran/seqinr/man/translate.html").
#' @param ... ...
#'
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_cu}}, \code{\link{get_du}}, \code{\link{get_freq}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_rscu(rgd.seq)
#' @exportMethod get_rscu
#' @name get_rscu
#' @rdname get_rscu-methods
setGeneric(
    name = "get_rscu",
    def = function(object, numcode=1, ...) standardGeneric("get_rscu")
)

#' @rdname get_rscu-methods
setMethod(
    f = "get_rscu",
    signature = "regioned_dna",
    definition = function(object, numcode) {
        dnaseq <- object@dnaseq[seq_len(length(object@dnaseq) - 1)]
        tmp <- oligonucleotideFrequency(dnaseq,
            width = 3, step = 3)
        tmp <- apply(tmp, 1, function(x) {
            r <- freq(x, numcode)
            r <- lapply(r, function(x) {
                return(x / ((1 / length(x)) * sum(x)))
            })
            names(r) <- NULL
            r <- unlist(r)[names(x)]
            return(r)
        })
        return(t(tmp))
    }
)

#' @rdname get_rscu-methods
setMethod(
    f = "get_rscu",
    signature = "DNAStringSet",
    definition = function(object, numcode) {
        tmp <-
            oligonucleotideFrequency(object, width = 3, step = 3)
        tmp <- apply(tmp, 1, function(x) {
            r <- freq(x, numcode)
            r <- lapply(r, function(x) {
                return(x / ((1 / length(x)) * sum(x)))
            })
            names(r) <- NULL
            r <- unlist(r)[names(x)]
            return(r)
        })
        return(t(tmp))
    }
)
