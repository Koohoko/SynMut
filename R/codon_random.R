#' Generate random synonymous mutations
#'
#' Generating random synonymous mutations (in user-defined region), with
#' optionally keeping/not keeping the original codon usage bias.
#'
#' @param object A regioned_dna object.
#' @param n Optional n parameter specifying what proportion of the codons to be
#'   mutate. Default value: \code{1}.
#' @param keep Logical parameter controling whether keeping the codon usage bias
#'   of the original sequence. Default value: \code{FALSE}.
#' @param ... ...
#'
#' @details This method randomly sample synonymous codons for \code{n} propotion
#'   of every mutable codons in the sequences. This process will be likely to
#'   alter the codon usage bias of the original sequences. However the
#'   \code{keep = TRUE} argument help to preserve the codon usage bias. It is
#'   done via the \code{synsequence} function in \code{seqinr} package. The
#'   \code{synsequence} function essentially swaps the position of the
#'   synonymous codons without introducing new codons into the original
#'   sequences.
#'
#' @return A regioned_dna object containing the mutants; Or a DNAStringSet
#'   object if the input is a DNAStringSet object.
#' @seealso \code{\link{input_seq}}, \code{\link{dinu_to}},
#'   \code{\link{codon_to}}, \code{\link{codon_mimic}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' set.seed(2019)
#' get_cu(codon_random(rgd.seq, n = 0.5))
#' get_cu(codon_random(rgd.seq))
#' @name codon_random
#' @rdname codon_random-methods
#'
#' @exportMethod codon_random
#' @importFrom seqinr c2s s2c synsequence splitseq syncodons
#' @importFrom Biostrings DNAStringSet
#' @import methods
#' @include regioned_dna_Class.R input_seq.R
setGeneric(
    name = "codon_random",
    def = function(object,
        n = 1,
        keep = FALSE,
        ...) {
        standardGeneric(f = "codon_random")
    }
)

#' @rdname codon_random-methods
setMethod(
    f = "codon_random",
    signature = "regioned_dna",
    definition = function(object, n, keep) {
        if (!is.logical(keep)) {
            stop("'keep' should be either TRUE or FALSE")
        }
        check.n <- all(n > 0, n <= 1)
        if (is.na(n)) {
            n = 1
        } else if (!check.n) {
            stop("n must be at range (0, 1]")
        } else {
            n = n
        }

        check.region <- all(is.na(object@region))
        seq <- convert_to_seq(object@dnaseq)
        seq.region <- extract_region(object, check.region)

        seq.mut <- mutation_random_main(seq.region, n, keep)

        seq.mut <- region_back(seq.mut, check.region, seq, object)

        return(new(
            "regioned_dna",
            dnaseq = seq.mut,
            region = object@region
        ))
    }
)

#' @rdname codon_random-methods
setMethod(
    f = "codon_random",
    signature = "DNAStringSet",
    definition = function(object, n, keep) {
        object <- c(object, DNAStringSet("ATG"))
        seq <- lapply(as.character(object), function(x) {
            splitseq(s2c(x))
        })

        seq.mut <- mutation_random_main(seq, n, keep)

        seq.mut <- seq.mut[seq_len(length(seq.mut) - 1)]
        seq.mut <- DNAStringSet(unlist(lapply(seq.mut, c2s)))

        return(seq.mut)
    }
)


# helper function ---------------------------------------------------------

mutation_random_main <- function(seq.region, n, keep){
    if (n != 1) {
        id <- lapply(seq.region, function(x) {
            id <- sample(seq_len(length(x)), round(length(x) * n))
        })
        seq.mut <- mapply(function(x, y) {
            x[y]
        }, seq.region, id, SIMPLIFY = FALSE)
    } else {
        seq.mut <- seq.region
    }

    if (keep == FALSE) {
        seq.mut <- lapply(seq.mut, function(x) {
            toupper(vapply(syncodons(x), function(x) {
                sample(x, 1)
            }, character(1)))
        })
    } else {
        seq.mut <- lapply(seq.mut, function(x) {
            splitseq(synsequence(s2c(c2s(x))))
        })
    }

    if (n != 1) {
        seq.mut <- mapply(function(x, y, z) {
            x[y] <- z
            return(x)
        }, seq.region, id, seq.mut, SIMPLIFY = FALSE)
    }

    return(seq.mut)
}

