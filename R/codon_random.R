#' Generate random synonymous mutations
#'
#' Generating radom synonymous mutations (in user-defined region), with
#' optionally keeping/not keeping the original codon usage bias.
#'
#' @param \code{object} A regioned_dna object.
#' @param \code{n} Optional n parameter specifying what proportion of the codons
#'   to be mutate. Default value: \code{NA}.
#' @param \code{keep} Logical parameter controling whether keeping the codon
#'   usage bias of the original sequence. Default value: \code{FALSE}.
#'
#' @return A regioned_dna object containing the mutants
#' @seealso \code{\link{input_seq}}, \code{\link{dinu_to}},
#'   \code{\link{codon_to}}, \code{\link{codon_mimic}},
#'   \code{\link{codon_mimic_dinu_to}}
#' @examples
#' setseed(2019)
#' get_cu(codon_random(rgd.seq, n = 0.5))
#' get_cu(codon_random(rgd.seq))
#' @name codon_random
#' @rdname codon_random-method
#'
#' @export
#' @importFrom seqinr c2s s2c synsequence splitseq syncodons
#' @include regioned_dna_Class.R input_seq.R
setGeneric(name = "codon_random",
  def = function(object, n = NA, keep = FALSE, ...){
    standardGeneric(f = "codon_random")
  })

#' @name codon_random
#' @rdname codon_random-method
setMethod(f = "codon_random", signature = "regioned_dna",
  definition = function(object, n, keep){
    if(!is.logical(keep)){
      stop("'keep' should be either TRUE or FALSE")
    }
    check.n <- all(n > 0, n < 1)
    if(is.na(n)){
      n = 1
    } else if(!check.n){
      stop("n must be at range (0, 1)")
    } else {
      n = n
    }

    check.region <- all(is.na(object@region))
    if(!check.region){
      seq <- sapply(as.character(object@dnaseq), function(x){splitseq(s2c(x))})
      seq.region <- mapply(function(x, y){
        return(x[y])
      }, seq, object@region)
    } else {
      seq.region <- sapply(as.character(object@dnaseq),
        function(x){splitseq(s2c(x))})
    }

    if(n != 1){
      id <- sapply(seq.region, function(x){
        id <- sample(seq_len(length(x)), round(length(x)*n))
      })
      seq.mut <- mapply(function(x, y){
        x[y]
      }, seq.region, id)
    } else {
      seq.mut <- seq.region
    }

    if(keep == FALSE){
      seq.mut <- sapply(seq.mut, function(x){
        toupper(sapply(syncodons(x), function(x){sample(x, 1)}))
      })
    } else {
      seq.mut <- sapply(seq.mut, function(x){
        splitseq(synsequence(s2c(c2s(x))))
      })
    }

    if(n != 1){
      seq.mut <- mapply(function(x, y, z){
        x[y] <- z
        return(x)
      }, seq.region, id, seq.mut)
    }

    if(!check.region){
      seq.mut <- mapply(function(x, y, z){
        x[y] <- z
        return(x)
      }, seq, object@region, seq.mut)
    }
    seq.mut <- Biostrings::DNAStringSet(sapply(seq.mut, c2s))

    return(input_seq(seq.mut, region = object@region))
  })
