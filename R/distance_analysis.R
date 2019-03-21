# dinu_dist ---------------------------------------------------------------

#' Calculating the dinucleotide usage difference between sequences
#'
#' We use a least squares approach to estimate the dinucleotide usage
#' difference between DNA sequences
#'
#' @param seq the input DNA sequnece of \code{DNAStringSet} or
#'   \code{regioned_dna} class.
#' @param ref the reference DNA sequnece of \code{DNAStringSet} or
#'   \code{regioned_dna} class.
#'
#' @return vector
#' @exportMethod dinu_dist
#'
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_cu(rgd.seq)
#'
#' mut.seq <- codon_random(rgd.seq)
#' dinu_dist(mut.seq, rgd.seq)
#' @name dinu_dist
#' @rdname dinu_dist-methods
setGeneric(name = "dinu_dist",
  def = function(seq, ref){standardGeneric("dinu_dist")})

#' @rdname dinu_dist-methods
setMethod(f = "dinu_dist",
  signature = signature("ANY"),
  definition = function(seq, ref){
    seq.du.freq <- apply(get_du(seq), 1, function(x){x/sum(x)})
    ref.du.freq <- apply(get_du(ref), 1, function(x){x/sum(x)})
    check <- all(ncol(seq.du.freq) != ncol(ref.du.freq), ncol(ref.du.freq) != 1)
    if(check){
      stop("the length of the reference sequences should be 1. of the same length
        as the input sequences OR 2. of length 1")
    }
    if(ncol(ref.du.freq) == 1){
      du.diff <- sweep(seq.du.freq, 1, ref.du.freq, "-")
    } else {
      du.diff <- seq.du.freq - ref.du.freq
    }

    du.dist <- apply(du.diff, 2, function(x){
      sqrt(mean(x*x))
    })
    return(du.dist)
  })


# codon_dist --------------------------------------------------------------
#' Calculating the codon usage difference between sequences
#'
#' We use a least squares approach to estimate the codon usage
#' difference between DNA sequences.
#'
#' @param seq the input DNA sequnece of \code{DNAStringSet} or
#'   \code{regioned_dna} class.
#' @param ref the reference DNA sequnece of \code{DNAStringSet} or
#'   \code{regioned_dna} class.
#'
#' @return vector
#' @exportMethod codon_dist
#'
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_cu(rgd.seq)
#'
#' mut.seq <- codon_random(rgd.seq)
#' codon_dist(mut.seq, rgd.seq)
#' mut.seq2 <- codon_random(rgd.seq, keep = TRUE)
#' codon_dist(mut.seq2, rgd.seq)
#' @name codon_dist
#' @rdname codon_dist-methods
setGeneric(name = "codon_dist",
  def = function(seq, ref){standardGeneric("codon_dist")})

#' @rdname codon_dist-methods
setMethod(f = "codon_dist",
  signature = "ANY",
  definition = function(seq, ref){
    seq.cu.freq <- get_freq(seq)
    ref.cu.freq <- get_freq(ref)
    check <- all(nrow(seq.cu.freq) != nrow(ref.cu.freq), nrow(ref.cu.freq) != 1)
    if(check){
      stop("the length of the reference sequences should be 1. of the same length
        as the input sequences OR 2. of length 1")
    }
    if(nrow(ref.cu.freq) == 1){
      cu.diff <- sweep(seq.cu.freq, 2, ref.cu.freq, "-")
    } else {
      cu.diff <- seq.cu.freq - ref.cu.freq
    }

    aa.list <- seqinr::ucoweight("")
    aa.list <- aa.list[names(aa.list) != "*"]

    cu.dist <- sapply(aa.list, function(x){
      cd.name.tmp <- toupper(names(x))
      if(length(cd.name.tmp) == 1){
        cu.diff.sub <- cu.diff[,toupper(colnames(cu.diff)) %in% cd.name.tmp]
        return(cu.diff.sub)
      } else {
        cu.diff.sub <- cu.diff[,toupper(colnames(cu.diff)) %in% cd.name.tmp]
        if(is.null(nrow(cu.diff.sub))){
          rslt <- sqrt(mean(cu.diff.sub*cu.diff.sub))
        } else {
          rslt <- apply(cu.diff.sub, 1, function(x){
            sqrt(mean(x*x))
          })
        }
        return(rslt)
      }
    })

    if(is.null(nrow(cu.dist))){
      return(mean(cu.dist))
    } else {
      return(rowMeans(cu.dist))
    }

  })


