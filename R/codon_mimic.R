#' Mimic a target codon usage bias
#'
#' Mutating the current DNA sequences in the regioned_dna object to mimic a
#' target codon usage pattern.
#'
#' @param object regioned_dna object
#' @param alt target codon usage vector or DNAStringSet object representing
#'   target codon usage
#' @param numcode The ncbi genetic code number for translation. Default value: \code{1}. Details please refer to \code{?seqinr::translate} ("https://rdrr.io/cran/seqinr/man/translate.html").
#' @param ... ...
#'
#' @details The ideas for \code{codon_mimic} is similar to
#'   \code{\link{codon_to}}: first extract the mutable regions and then do the
#'   mutation. However the codons in the fixed (not mutable) regions will also
#'   alter the final codon usage, thus we have to adjust for the fixed codons
#'   when introducing sysnonymous codons. The ideal deisgn for
#'   \code{codon_mimic} is not unique as the swap between positions of the
#'   synonymous codons will not change the codon usage bias.
#'
#'   Details pleas refer to: https://koohoko.github.io/SynMut/algorithm.html
#'
#' @return regioned_dna
#' @exportMethod codon_mimic
#' @include regioned_dna_Class.R
#' @seealso \code{\link{input_seq}}, \code{\link{codon_to}},
#'   \code{\link{codon_random}}, \code{\link{dinu_to}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' target <- get_cu(rgd.seq)[2,]
#' new <- codon_mimic(rgd.seq, alt = target)
#' get_cu(new) - get_cu(rgd.seq)
#'
#' target <- Biostrings::DNAStringSet("TTGAAAA-CTC-N--AAG")
#' new <- codon_mimic(rgd.seq, alt = target)
#' get_cu(new) - get_cu(rgd.seq)
#' get_freq(new) - get_freq(rgd.seq)
#' get_rscu(new) - get_rscu(rgd.seq)
#' @importFrom Biostrings oligonucleotideFrequency
#' @import methods
#' @import Biostrings
#' @name codon_mimic
#' @rdname codon_mimic-methods
setGeneric(
    name = "codon_mimic",
    def = function(object, alt, numcode = 1, ...) standardGeneric("codon_mimic")
)

#' @rdname codon_mimic-methods
setMethod(
    f = "codon_mimic",
    signature = signature(object = "regioned_dna", alt = "vector"),
    definition = function(object, alt, numcode) {
        if (length(alt) != 64) {
            stop("The length of the target codon usage must be 64")
        }
        alt <- as.numeric(alt)
        names(alt) <- sort(names(GENETIC_CODE))
        result <-
            codon.mimic(alt, dna.seq = object@dnaseq, region = object@region, numcode)
        return(result)
    }
)

#' @rdname codon_mimic-methods
setMethod(
    f = "codon_mimic",
    signature = signature(object = "regioned_dna", alt = "DNAStringSet"),
    definition = function(object, alt, numcode) {
        if (length(alt) > 1) {
            warning("only the first one sequence in the the target was used")
        }
        cu.target <- oligonucleotideFrequency(alt, width = 3, step = 3)[1, ]
        result <- codon.mimic(cu.target,
            dna.seq = object@dnaseq,
            region = object@region,
            numcode)
        return(result)
    }
)

#' @rdname codon_mimic-methods
setMethod(
    f = "codon_mimic",
    signature = signature(object = "DNAStringSet", alt = "DNAStringSet"),
    definition = function(object, alt, numcode) {
        if (length(alt) > 1) {
            warning("only the first one sequence in the the target was used")
        }
        cu.target <- oligonucleotideFrequency(alt, width = 3, step = 3)[1, ]
        object <- input_seq(object)
        result <- codon.mimic(cu.target,
            dna.seq = object@dnaseq,
            region = object@region,
            numcode)
        return(result)
    }
)

codon.mimic <- function(cu.target, dna.seq, region, numcode) {
    check.region <- all(is.na(region))
    seq <- convert_to_seq(dna.seq)
    if (!check.region) {
        #have region
        seq.mut <- codon.mimic.region(cu.target, dna.seq, region, numcode)
    } else {
        #no restriction of region
        seq.mut <- codon.mimic.no.region(cu.target, dna.seq, region, numcode)
    }
    # merge region ------------------------------------------------------------

    if (!check.region) {
        seq.mut <- mapply(function(x, y, z) {
            x[y] <- z
            return(x)
        }, seq, region, seq.mut, SIMPLIFY = FALSE)
    }
    seq.mut <- DNAStringSet(unlist(lapply(seq.mut, c2s)))
    return(new("regioned_dna",
        dnaseq = seq.mut,
        region = region))
}


# helper function -------------------------------------------------------

codon.count <- function(x, numcode) {
    gen_code <- sapply(names(GENETIC_CODE), function(x){
        seqinr::translate(strsplit(x, "")[[1]], numcode = numcode)
    })
    gen_code <- gen_code[order(names(gen_code))]
    base::split(x, gen_code)
}

freq <- function(x, numcode) {
    lapply(codon.count(x, numcode), function(x) {
        x / sum(x)
    })
}

mut.assign.back <- function(seq.region, mut.need, numcode){
    seq.mut <- lapply(seq_along(seq.region), function(i) {
        seq.tmp <- seq.region[[i]]
        seq.tmp.aa <- seqinr::translate(s2c(c2s(seq.tmp)), numcode = numcode)
        aa.names <- names(table(seq.tmp.aa))
        mut.need.tmp <- mut.need[[i]]
        for (j in seq_along(aa.names)) {
            for.sample <- as.character(which(seq.tmp.aa == aa.names[j]))
            pos.tmp <- as.numeric(sample(for.sample))
            mut.cd.tmp <- round(mut.need.tmp[[aa.names[j]]])
            if (!any(is.na(mut.cd.tmp)) & any(mut.cd.tmp>1)) {
                suppressWarnings(seq.tmp[pos.tmp] <-
                        rep(names(mut.cd.tmp), mut.cd.tmp))
            }
        }
        return(seq.tmp)
    })
}

codon.mimic.region <- function(cu.target, dna.seq, region, numcode){
    seq <- convert_to_seq(dna.seq)
    seq.region <- mapply(function(x, y) {
        return(x[y])
    }, seq, region, SIMPLIFY = FALSE)
    seq.fixed <- mapply(function(x, y) {
        return(x[!y])
    }, seq, region, SIMPLIFY = FALSE)

    cu.fixed <- get_cu(DNAStringSet(unlist(lapply(seq.fixed, c2s))))
    cu.ori <- get_cu(dna.seq)
    freq.target <- freq(cu.target, numcode)

    mut.need <- lapply(seq_along(seq), function(i) {
        count.ori <- codon.count(cu.ori[i,], numcode)
        count.fixed <- codon.count(cu.fixed[i,], numcode)
        mut.usage <- mapply(function(x, y) {
            sum(x) * y
        }, count.ori, freq.target, SIMPLIFY = FALSE)
        mut.need <- mapply(function(x, y) {
            x - y
        }, mut.usage, count.fixed, SIMPLIFY = FALSE)
        mut.need <- lapply(mut.need, function(x) {
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
    seq.mut <- mut.assign.back(seq.region, mut.need, numcode)
    names(seq.mut) <- names(seq.region)
    return(seq.mut)
}

codon.mimic.no.region <- function(cu.target, dna.seq, region, numcode){
    seq.region <- convert_to_seq(dna.seq)
    freq.target <- freq(cu.target, numcode)

    seq.mut <- lapply(seq_along(seq.region), function(i) {
        seq.tmp <- seq.region[[i]]
        seq.tmp.aa <- seqinr::translate(s2c(c2s(seq.tmp)), numcode=numcode)
        # print(numcode)
        aa.names <- names(table(seq.tmp.aa))
        for (j in seq_along(aa.names)) {
            for.sample <- as.character(which(seq.tmp.aa == aa.names[j]))
            pos.tmp <- as.numeric(sample(for.sample))
            mut.cd.tmp <- round(length(pos.tmp) * freq.target[[aa.names[j]]])
            if (!any(is.na(mut.cd.tmp))) {
                suppressWarnings(seq.tmp[pos.tmp] <-
                        rep(names(mut.cd.tmp), mut.cd.tmp))
            }
        }
        return(seq.tmp)
    })
    names(seq.mut) <- names(seq.region)

    return(seq.mut)
}
