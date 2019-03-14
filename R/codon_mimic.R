#' Mimic a target codon usage bias
#'
#' Mutating the current DNA sequences in the regioned_dna object to mimic a
#' target codon usage pattern.
#'
#' @param object regioned_dna object
#' @param alt target codon usage vector or DNAStringSet object representing
#'   target codon usage
#'
#' @return regioned_dna
#' @exportMethod codon_mimic
#' @include regioned_dna_Class.R
#' @seealso \code{\link{input_seq}}, \code{\link{codon_to}},
#'   \code{\link{codon_random}}, \code{\link{dinu_to}}
#' @examples
#' target <- get_cu(rgd.seq)[2,]
#' new <- codon_mimic(rgd.seq, alt = target)
#' get_cu(new) - get_cu(rgd.seq)
#'
#' target <- Biostrings::DNAStringSet("TTGAAAA-CTC-N--AAG")
#' new <- codon_mimic(rgd.seq, alt = target)
#' get_cu(new) - get_cu(rgd.seq)
#' get_freq(new) - get_freq(rgd.seq)
#' get_rscu(new) - get_rscu(rgd.seq)
#' @name codon_mimic
#' @rdname codon_mimic-methods
setGeneric(
  name = "codon_mimic",
  def = function(object, alt, ...) {
    standardGeneric("codon_mimic")
  }
)


#' @rdname codon_mimic-methods
#' @name codon_mimic
setMethod(
  f = "codon_mimic",
  signature = signature(object = "regioned_dna",
    alt = "vector"),
  definition = function(object, alt) {
    if (length(alt) != 64) {
      stop("The length of the target codon usage must be 64")
    }
    alt <- as.numeric(alt)
    names(alt) <- sort(names(Biostrings::GENETIC_CODE))
    result <- codon.mimic(alt,
      dna.seq = object@dnaseq, region = object@region)
    return(result)
  }
)

#' @name codon_mimic
#' @rdname codon_mimic-methods
setMethod(
  f = "codon_mimic",
  signature = signature(object = "regioned_dna",
    alt = "DNAStringSet"),
  definition = function(object, alt) {
    if (length(alt) > 1) {
      warning("only the first one sequence in the the target was used")
    }
    cu.target <- Biostrings::oligonucleotideFrequency(alt,
      width = 3, step = 3)[1, ]
    result <- codon.mimic(cu.target,
      dna.seq = object@dnaseq,
      region = object@region)
    return(result)
  }
)

codon.mimic <- function(cu.target, dna.seq, region) {
  check.region <- all(is.na(region))
  if (!check.region) {
    seq <- sapply(as.character(dna.seq), function(x) {
      splitseq(s2c(x))
    })
    seq.region <- mapply(function(x, y) {
      return(x[y])
    }, seq, region, SIMPLIFY = FALSE)
    seq.fixed <- mapply(function(x, y) {
      return(x[!y])
    }, seq, region, SIMPLIFY = FALSE)

    cu.fixed <-
      get_cu(Biostrings::DNAStringSet(sapply(seq.fixed, c2s)))
    cu.ori <- get_cu(dna.seq)
    freq.target <- freq(cu.target)

    mut.need <- lapply(seq_along(seq), function(i) {
      count.ori <- codon.count(cu.ori[i, ])
      count.fixed <- codon.count(cu.fixed[i, ])
      mut.usage <- mapply(function(x, y) {
        sum(x) * y
      }, count.ori, freq.target, SIMPLIFY = FALSE)
      mut.need <- mapply(function(x, y) {
        x - y
      },
        mut.usage, count.fixed, SIMPLIFY = FALSE)
      mut.need <-
        sapply(mut.need, function(x) {
          #fix the negative needs
          x.negative <- x[x < 0]
          x.positive <- x[x > 0]
          if (length(x.negative) > 0) {
            return((1 + sum(abs(
              x.negative
            )) / sum(x.positive)) * x.positive)
          } else{
            return(x)
          }
        })
    })

    seq.mut <- sapply(seq_along(seq.region), function(i) {
      seq.tmp <- seq.region[[i]]
      seq.tmp.aa <- seqinr::translate(s2c(c2s(seq.tmp)))
      aa.names <- names(table(seq.tmp.aa))
      mut.need.tmp <- mut.need[[i]]
      for (j in seq_along(aa.names)) {
        for.sample <- as.character(which(seq.tmp.aa == aa.names[j]))
        pos.tmp <- as.numeric(sample(for.sample))
        mut.cd.tmp <-  round(mut.need.tmp[[aa.names[j]]])
        if (!any(is.na(mut.cd.tmp))) {
          suppressWarnings(seq.tmp[pos.tmp] <-
              rep(names(mut.cd.tmp), mut.cd.tmp))
        }
      }
      return(seq.tmp)
    })
    names(seq.mut) <- names(seq.region)

  } else {
    #no restriction of region
    seq.region <- sapply(as.character(dna.seq),
      function(x) {
        splitseq(s2c(x))
      })
    freq.target <- freq(cu.target)

    seq.mut <- sapply(seq_along(seq.region), function(i) {
      seq.tmp <- seq.region[[i]]
      seq.tmp.aa <- seqinr::translate(s2c(c2s(seq.tmp)))
      aa.names <- names(table(seq.tmp.aa))
      for (j in seq_along(aa.names)) {
        for.sample <- as.character(which(seq.tmp.aa == aa.names[j]))
        pos.tmp <- as.numeric(sample(for.sample))
        mut.cd.tmp <-
          round(length(pos.tmp) * freq.target[[aa.names[j]]])
        if (!any(is.na(mut.cd.tmp))) {
          suppressWarnings(seq.tmp[pos.tmp] <-
              rep(names(mut.cd.tmp), mut.cd.tmp))
        }
      }
      return(seq.tmp)
    })
    names(seq.mut) <- names(seq.region)
  }
  # merge region ------------------------------------------------------------

  if (!check.region) {
    seq.mut <- mapply(function(x, y, z) {
      x[y] <- z
      return(x)
    }, seq, region, seq.mut, SIMPLIFY = FALSE)
  }
  seq.mut <- Biostrings::DNAStringSet(sapply(seq.mut, c2s))
  return(new("regioned_dna",
    dnaseq = seq.mut,
    region = region))
}


# internal function -------------------------------------------------------

codon.count <- function(x) {
  base::split(x, seqinr::translate(s2c(c2s(names(
    x
  )))))
}

freq <- function(x) {
  lapply(codon.count(x), function(x) {
    x / sum(x)
  })
}
