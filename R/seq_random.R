#' Generate \code{n} random DNA sequnces of length \code{m}
#'
#' Generate \code{n} random DNA sequnces of length \code{m}, optional exclude
#' stop codons.
#'
#' @param n the number of the output sequence(s).
#' @param m the length of the ouput sequence(s). Eighter a fixed number or a
#'   vector of different numbers.
#' @param no.stop.codon Default FALSE. If TRUE, the stop codons in the frame 1
#'   would be substituted to another random codon.
#' @param ... ...
#'
#' @return a DNAStringSet object
#' @exportMethod seq_random
#'
#' @examples
#' seq_random(n = 1, m = 99)
#' seq_random(n = 10, m = 30)
#' seq_random(n = 10, m = 1:10)
#' seq.nsc <- seq_random(n = 10, m = 100, no.stop.codon = TRUE)
#' get_cu(seq.nsc)
#' @name seq_random
#' @rdname seq_random-methods
setGeneric(
  name = "seq_random",
  def = function(n = 1,
    m,
    no.stop.codon = FALSE,
    ...) {
    standardGeneric("seq_random")
  }
)

#' @rdname seq_random-methods
setMethod(
  f = "seq_random",
  signature = signature(n = "numeric", m = "numeric"),
  definition = function(n, m, no.stop.codon) {
    check.m.length <- length(m) == 1
    if (check.m.length) {
      tmp <- sapply(seq_len(n), function(x) {
        paste0(sample(
          x = c("a", "t", "g", "c"),
          size = m,
          replace = TRUE
        ), collapse = "")
      })
    } else {
      tmp <- mapply(function(x, y) {
        paste0(sample(
          x = c("a", "t", "g", "c"),
          size = y,
          replace = TRUE
        ), collapse = "")
      }, seq_len(n), m)
    }
    if (no.stop.codon) {
      stop.codons <- c("taa", "tga", "tag")
      ok.codons <- setdiff(seqinr::words(), stop.codons)
      tmp <- sapply(tmp, function(x) {
        sst <- strsplit(x, "")[[1]]
        if (length(sst) %% 3 != 0) {
          mut.length <- floor(length(sst) / 3) * 3
          sst.mut <- sst[seq_len(mut.length)]
          sst.tail <- sst[(mut.length + 1):length(sst)]
        } else{
          sst.mut <- sst
          sst.tail <- ""
        }
        codons <- paste0(sst.mut[c(TRUE, FALSE, FALSE)],
          sst.mut[c(FALSE, TRUE, FALSE)],
          sst.mut[c(FALSE, FALSE, TRUE)])
        stop.id <- which(codons %in% stop.codons)
        if (length(stop.id) > 0) {
          codons[stop.id] <- sapply(seq_along(stop.id), function(x) {
            sample(ok.codons, size = 1)
          })
          return(paste0(paste0(codons, collapse = ""), sst.tail, collapse = ""))
        } else {
          return(x)
        }
      }, USE.NAMES = FALSE)
    }
    Biostrings::DNAStringSet(tmp)
  }
)
