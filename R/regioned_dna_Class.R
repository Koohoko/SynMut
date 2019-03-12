.onAttach <- function(libname, pkgname){
  packageStartupMessage("SynMut: Tools for designing Synonymous codon/dinucleotide Mutants in Virus Research
-----------------------------------------------------------------------------------
For disccussion, please access: https://github.com/Koohoko/SynMut/issues ")
}

# Class definition --------------------------------------------------------

#' An S4 class to record DNA sequences and variable regions for mutations
#'
#' Recording codon DNA sequences and region.
#'
#' @slot \code{dnaseq} a DNAStingSet object recording the sequence(s)
#' @slot \code{region} a list specifying paticular regions in the sequences
#'   allowed to be mutated
#'
#' @importFrom methods setClass
#' @seealso \code{\link{input_seq}}, \code{\link{get_cu}},
#'   \code{\link{get_region}}
#' @export
#' @author Haogao Gu
setClass(Class = "regioned_dna",
  slots = c(
    dnaseq = "DNAStringSet",
    region = "list"
  )
)

setValidity("regioned_dna",
  function(object){
    check.na <- all(is.na(object@region))
    dnaseq.n <- length(object@dnaseq)
    region.n <- length(object@region)
    check.3 <- all(sapply(object@dnaseq, length) %% 3 == 0)
    if(!check.na){
      check.length <- dnaseq.n == region.n
      check.num <- sapply(object@dnaseq, length)/3 == sapply(object@region, length)
    } else {
      check.num <- TRUE
      check.length <- TRUE
    }
    if(!check.3){
      stop("the length of the dna sequences must be divisible by 3")
    }
    if(all(check.num, check.length)){TRUE} else{
      "dnaseq and region must have the same length."
    }
  })

##viewer
setMethod(f = "show", signature = "regioned_dna",
  definition = function(object){
    cat("An object of class ", class(object), "\n", sep = "")
    cat("Number of sequences: ", length(object@dnaseq), "\n", sep = "")
  })


# Accessor ----------------------------------------------------------------

#' Get the variable region
#'
#' Access the variable regions
#'
#' @param object regioned_dna
#' @return list
#' @seealso \code{\link{input_seq}}, \code{\link{get_cu}}
#' @examples
#' get_region(rgd.seq)
#' @export
#' @name get_region
#' @rdname get_region-method
setGeneric(name = "get_region",
  def = function(object, ...){standardGeneric("get_region")})

#' @name get_region
#' @rdname get_region-method
setMethod(f = "get_region", signature = "regioned_dna",
  definition = function(object){
    object@region
  })

#' Get the DNAStringSet data
#'
#' Access the DNA sequence data in DNAStringSet.
#'
#' @param object A regioned_dna object.
#' @param ...
#'
#' @return DNAStringSet
#' @export
#' @name get_dna
#' @rdname get_dna-method
#'
#' @examples
#' get_dna(rgd.seq)
setGeneric(name = "get_dna",
  def = function(object, ...){standardGeneric("get_dna")})

#' @name get_dna
#' @rdname get_dna-method
setMethod(f = "get_dna", signature = "regioned_dna",
  definition = function(object){
    return(object@dnaseq)
  })


# simple calculation ------------------------------------------------------

#' Get codon usage matrix
#'
#' Access the codon usage matrix
#'
#' @param object regioned_dna / DNAStringSet
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_du}}, \code{\link{get_freq}}, \code{\link{get_rscu}}
#' @examples
#' get_cu(rgd.seq)
#' @export
#' @name get_cu
#' @rdname get_cu-method
setGeneric(name = "get_cu",
  def = function(object, ...){standardGeneric("get_cu")})

#' @name get_cu
#' @rdname get_cu-method
setMethod(f = "get_cu", signature = "regioned_dna",
  definition = function(object){
    Biostrings::oligonucleotideFrequency(object@dnaseq, width = 3, step = 3)
  })

#' @name get_cu
#' @rdname get_cu-method
setMethod(f = "get_cu", signature = "DNAStringSet",
  definition = function(object){
    Biostrings::oligonucleotideFrequency(object, width = 3, step = 3)
  })

#' Get dinucleotide usage matrix
#'
#' Access the dinucleotide usage matrix
#'
#' @param object regioned_dna / DNAStringSet
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_cu}}, \code{\link{get_freq}}, \code{\link{get_rscu}}
#' @examples
#' get_du(rgd.seq)
#' @export
#' @name get_du
#' @rdname get_du-method
setGeneric(name = "get_du",
  def = function(object, ...){standardGeneric("get_du")})

#' @name get_du
#' @rdname get_du-method
setMethod(f = "get_du", signature = "regioned_dna",
  definition = function(object){
    Biostrings::oligonucleotideFrequency(object@dnaseq, width = 2, step = 2)
  })

#' @name get_du
#' @rdname get_du-method
setMethod(f = "get_du", signature = "DNAStringSet",
  definition = function(object){
    Biostrings::oligonucleotideFrequency(object, width = 2, step = 2)
  })


#' Get Relative Codon Usage Frequency of synonymous codons
#'
#' Access the Relative Synonymous Codon Usage
#'
#' @param object regioned_dna / DNAStringSet / codon usage matrix (vector)
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_cu}}, \code{\link{get_du}}, \code{\link{get_rscu}}
#' @examples
#' get_freq(rgd.seq)
#' @export
#' @name get_freq
#' @rdname get_freq-method
setGeneric(name = "get_freq",
  def = function(object, ...){standardGeneric("get_freq")})

#' @name get_freq
#' @rdname get_freq-method
setMethod(f = "get_freq", signature = "regioned_dna",
  definition = function(object){
    tmp <- Biostrings::oligonucleotideFrequency(object@dnaseq,
      width = 3, step = 3)
    tmp <- apply(tmp, 1, function(x){
      r <- freq(x)
      names(r) <- NULL
      r <- unlist(r)[names(x)]
      return(r)
    })
    return(t(tmp))
  })

#' @name get_freq
#' @rdname get_freq-method
setMethod(f = "get_freq", signature = "DNAStringSet",
  definition = function(object){
    tmp <- Biostrings::oligonucleotideFrequency(object, width = 3, step = 3)
    tmp <- apply(tmp, 1, function(x){
      r <- freq(x)
      names(r) <- NULL
      r <- unlist(r)[names(x)]
      return(r)
    })
    return(t(tmp))
  })

#' @name get_freq
#' @rdname get_freq-method
setMethod(f = "get_freq", signature = "matrix",
  definition = function(object){
    tmp <- apply(object, 1, function(x){
      r <- freq(x)
      names(r) <- NULL
      r <- unlist(r)[names(x)]
      return(r)
    })
    return(t(tmp))
  })

#' @name get_freq
#' @rdname get_freq-method
setMethod(f = "get_freq", signature = "vector",
  definition = function(object){
    r <- freq(object)
    names(r) <- NULL
    r <- unlist(r)[names(object)]
    return(r)
  })

#' Get Relative Synonymous Codon Usage (rscu) of synonymous codons
#'
#' Access the Relative Synonymous Codon Usage rscu
#'
#' @param object regioned_dna / DNAStringSet / codon usage matrix (vector)
#' @return matrix
#' @seealso \code{\link{input_seq}}, \code{\link{get_region}},
#'   \code{\link{get_cu}}, \code{\link{get_du}}, \code{\link{get_freq}}
#' @examples
#' get_rscu(rgd.seq)
#' @export
#' @name get_rscu
#' @rdname get_rscu-method
setGeneric(name = "get_rscu",
  def = function(object, ...){standardGeneric("get_rscu")})

#' @name get_rscu
#' @rdname get_rscu-method
setMethod(f = "get_rscu", signature = "regioned_dna",
  definition = function(object){
    tmp <- Biostrings::oligonucleotideFrequency(object@dnaseq,
      width = 3, step = 3)
    tmp <- apply(tmp, 1, function(x){
      r <- freq(x)
      r <- sapply(r, function(x){
        return(x/((1/length(x)) * sum(x)))
      })
      names(r) <- NULL
      r <- unlist(r)[names(x)]
      return(r)
    })
    return(t(tmp))
  })

#' @name get_rscu
#' @rdname get_rscu-method
setMethod(f = "get_rscu", signature = "DNAStringSet",
  definition = function(object){
    tmp <- Biostrings::oligonucleotideFrequency(object, width = 3, step = 3)
    tmp <- apply(tmp, 1, function(x){
      r <- freq(x)
      r <- sapply(r, function(x){
        return(x/((1/length(x)) * sum(x)))
      })
      names(r) <- NULL
      r <- unlist(r)[names(x)]
      return(r)
    })
    return(t(tmp))
  })


# internal function -------------------------------------------------------

codon.count <- function(x){
  base::split(x, seqinr::translate(s2c(c2s(names(x)))))
}

freq <- function(x){
  lapply(codon.count(x), function(x){
    x/sum(x)
  })
}