#' Maximize or minimize the usage of certain dinucleotide.
#'
#' Input string of a dinucleotide to either the "max.dinu = " or "min.codon = "
#' parameter to maximize or minimize the usage of certain codon in the sequence.
#' Using a greedy algorithm with priority given to dinucleotide12 or
#' dinucleotide23.
#'
#' @param object A regioned_dna object.
#' @param max.dinu A string of a dinucleotide.
#' @param min.dinu A string of a dinucleotide.
#' @param keep A logical varibale stating if the codon usage of the original
#'   sequences should be keep. Default: False.
#' @param ... ...
#' @details The detail strategy for this function please refer to:
#' https://koohoko.github.io/SynMut/algorithm.html
#' @return regioned_dna
#' @seealso \code{\link{input_seq}}, \code{\link{codon_to}},
#'   \code{\link{codon_random}}, \code{\link{codon_mimic}}
#' @examples
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#' get_du(dinu_to(rgd.seq, max.dinu = "cg")) - get_du(rgd.seq)
#' get_du(dinu_to(rgd.seq, min.dinu = "AA")) - get_du(rgd.seq)
#' get_du(dinu_to(rgd.seq, max.dinu = "cg", keep = TRUE)) - get_du(rgd.seq)
#' get_cu(dinu_to(rgd.seq, max.dinu = "CG", keep = TRUE)) - get_cu(rgd.seq)
#' @include regioned_dna_Class.R
#' @exportMethod dinu_to
#' @importFrom stringr str_detect
#' @name dinu_to
#' @rdname dinu_to-methods
setGeneric(
    name = "dinu_to",
    def = function(object,
        max.dinu = NA,
        min.dinu = NA,
        keep = FALSE,
        ...) {
        standardGeneric("dinu_to")
    }
)

#' @rdname dinu_to-methods
setMethod(
    f = "dinu_to",
    signature = "regioned_dna",
    definition = function(object, max.dinu, min.dinu, keep) {
        max.dinu <- toupper(max.dinu)
        min.dinu <- toupper(min.dinu)
        dinu.input <- c(max.dinu, min.dinu)
        dinu.input <- dinu.input[!is.na(dinu.input)]
        if (length(dinu.input) == 0) {
            stop("please specify at least one parameter of
                 'max.dinu' or 'min.dinu'")
        }
        if (all(!is.na(c(max.dinu, min.dinu)))) {
            stop("Only one parameter is supported between
                 'max.dinu' and 'min.dinu'")
        }
        if (length(dinu.input) > 2) {
            stop("Only one dinucleotide is allowed for
                 input in 'max.dinu' or 'min.dinu'")
        }
        check.valid.dinu <-
            toupper(dinu.input) %in% toupper(seqinr::words(2))
        if (!check.valid.dinu) {
            stop("Please input valid DNA dinucleotide")
        }

        # extract region ------------------------------------------------------

        region <- object@region
        check.region <- all(is.na(region))
        if (!check.region) {
            # have region
            seq <- lapply(as.character(object@dnaseq), function(x) {
                splitseq(s2c(x))
            })
            seq.region <- mapply(function(x, y) {
                return(x[y])
            }, seq, region, SIMPLIFY = FALSE)
        } else {
            #no region
            seq.region <- lapply(as.character(object@dnaseq),
                function(x) {
                    splitseq(s2c(x))
                })
        }

        # mutation ------------------------------------------------------------

        codon.list <-
            lapply(seqinr::ucoweight(''), function(x) {
                toupper(names(x))
            })
        codon.list.alt <-
            codon.list[vapply(codon.list, function(x) {
                length(x) > 1
            }, logical(1))]

        if (keep == FALSE) {
            seq.mut <-
                dinu_to.not.keep(seq.region, codon.list.alt, max.dinu, min.dinu)
        } else {
            seq.mut <- seq.region
        }

        # synonmous codon interchange -----------------------------------------

        if (!check.region) {
            seq.mut <- mapply(function(x, y, z) {
                x[y] <- z
                return(x)
            }, seq, region, seq.mut, SIMPLIFY = FALSE)
        }

        seq.mut <-
            dinu_to.keep(check.region,
                seq.mut,
                max.dinu,
                min.dinu,
                codon.list.alt,
                region)

        # merge region --------------------------------------------------------

        seq.mut <- Biostrings::DNAStringSet(unlist(lapply(seq.mut, c2s)))
        return(methods::new(
            "regioned_dna",
            dnaseq = seq.mut,
            region = object@region
        ))

    }
)

###########################################################################
# internal_function -------------------------------------------------------
###########################################################################

#1: not keep codon usage, mutate at all 12, 23 and 31 regions
dinu_to.not.keep <-
    function(seq.region,
        codon.list.alt,
        max.dinu,
        min.dinu) {
        ## get optimal codon
        if (!is.na(max.dinu)) {
            #max.dinu
            codon.list.optm <- lapply(codon.list.alt, function(x) {
                filter1 <- x[str_detect(x, pattern = max.dinu)]
                if (length(filter1) == 1) {
                    return(filter1)
                } else if (length(filter1) < 1) {
                    pattern2.1 <- paste0("^", substr(max.dinu, 2, 2))
                    pattern2.2 <-
                        paste0(substr(max.dinu, 1, 1), "$")
                    filter2.1 <-
                        x[str_detect(x, pattern = pattern2.1)]
                    filter2.2 <-
                        x[str_detect(x, pattern = pattern2.2)]
                    filter2 <- c(filter2.1, filter2.2)
                    if (length(filter2) < 1) {
                        return(sample(x, 1))
                    } else{
                        table.sorted <- sort(table(filter2), decreasing = TRUE)
                        names.ordered <- names(table.sorted)
                        choice.good <-
                            names.ordered[table.sorted == table.sorted[1]]
                        return(sample(choice.good, 1))
                    }
                } else {
                    pattern2.1 <- paste0("^", substr(max.dinu, 2, 2))
                    pattern2.2 <-
                        paste0(substr(max.dinu, 1, 1), "$")
                    filter2.1 <-
                        x[str_detect(x, pattern = pattern2.1)]
                    filter2.2 <-
                        x[str_detect(x, pattern = pattern2.2)]
                    filter2 <- c(filter2.1, filter2.2)
                    if (length(filter2) < 1) {
                        return(sample(x, 1))
                    } else{
                        table.sorted <- sort(table(filter2), decreasing = TRUE)
                        names.ordered <- names(table.sorted)
                        choice.good <-
                            names.ordered[table.sorted == table.sorted[1]]
                        return(sample(choice.good, 1))
                    }
                }
            })
        } else {
            #min.dinu
            codon.list.optm <- lapply(codon.list.alt, function(x) {
                filter1 <-
                    x[!str_detect(x, pattern = min.dinu)] #difference
                if (length(filter1) == 1) {
                    return(filter1)
                } else if (length(filter1) < 1) {
                    pattern2.1 <- paste0("^", substr(min.dinu, 2, 2))
                    pattern2.2 <-
                        paste0(substr(min.dinu, 1, 1), "$")
                    filter2.1 <-
                        x[!str_detect(x, pattern = pattern2.1)]
                    filter2.2 <-
                        x[!str_detect(x, pattern = pattern2.2)]
                    filter2 <- c(filter2.1, filter2.2)
                    if (length(filter2) < 1) {
                        return(sample(x, 1))
                    } else{
                        table.sorted <- sort(table(filter2), decreasing = TRUE)
                        names.ordered <- names(table.sorted)
                        choice.good <-
                            names.ordered[table.sorted == table.sorted[1]]
                        return(sample(choice.good, 1))
                    }
                } else {
                    pattern2.1 <- paste0("^", substr(min.dinu, 2, 2))
                    pattern2.2 <-
                        paste0(substr(min.dinu, 1, 1), "$")
                    filter2.1 <-
                        x[!str_detect(x, pattern = pattern2.1)]
                    filter2.2 <-
                        x[!str_detect(x, pattern = pattern2.2)]
                    filter2 <- c(filter2.1, filter2.2)
                    if (length(filter2) < 1) {
                        return(sample(x, 1))
                    } else{
                        table.sorted <- sort(table(filter2), decreasing = TRUE)
                        names.ordered <- names(table.sorted)
                        choice.good <-
                            names.ordered[table.sorted == table.sorted[1]]
                        return(sample(choice.good, 1))
                    }
                }
            })
        }

        ## mutate codons to optimals
        seq.mut <- seq.region
        for (i in seq_along(codon.list.optm)) {
            codons.alt <- codon.list.alt[[i]]
            codons.alt <-
                codons.alt[codons.alt != codon.list.optm[[i]]]
            seq.mut <- lapply(seq.mut, function(x) {
                x[x %in% codons.alt] <- codon.list.optm[[i]]
                return(x)
            })
        }

        return(seq.mut)
    }

#2: keep codon usage, mutate at only 31 regions
dinu_to.keep <-
    function(check.region,
        seq,
        max.dinu,
        min.dinu,
        codon.list.alt,
        region) {
        if (!check.region) {
            #not regioned data
            seq.mut <- dinu_to.keep.no.region(codon.list.alt, seq,
                max.dinu, min.dinu)
        } else {
            #regioned data
            seq.mut <-
                dinu_to.keep.region(codon.list.alt, seq, region,
                    max.dinu, min.dinu)
        }
        return(seq.mut)
    }

#2.1
dinu_to.keep.no.region <- function(codon.list.alt,
    seq,
    max.dinu,
    min.dinu) {
    check.max <- !is.na(max.dinu)
    if (check.max) {
        dinu.target <- max.dinu
    } else {
        dinu.target <- min.dinu
    }
    seq.aa <- lapply(seq, function(x) {
        seqinr::translate(s2c(c2s(x)))
    })
    seq.mut <- seq

    #rearrange all synonymous codons in aa order
    for (i in seq_along(codon.list.alt)) {
        aa.tmp <- names(codon.list.alt)[i]
        codon.tmp <- codon.list.alt[[i]]
        #check whether need to rearrange
        codon.target <-
            codon.tmp[substr(codon.tmp, 3, 3) == substr(dinu.target, 1, 1)]
        check.need <- length(codon.target) > 0
        if (check.need) {
            pos.all <- lapply(seq.aa, function(x) {
                which(x == aa.tmp)
            })
            codon.all <- mapply(function(x, y) {
                x[y]
            }, seq, pos.all, SIMPLIFY = FALSE)
            id.target <- lapply(codon.all, function(x) {
                which(x %in% codon.target)
            })

            if (check.max) {
                id.good <- lapply(seq, function(x) {
                    #max
                    codon2nd <- x[which(x %in% codon.tmp) + 1]
                    which(substr(codon2nd, 1, 1) == substr(dinu.target, 2, 2))
                })
            } else {
                id.good <- lapply(seq, function(x) {
                    #min
                    codon2nd <- x[which(x %in% codon.tmp) + 1]
                    #note this difference between min and max
                    which(substr(codon2nd, 1, 1) != substr(dinu.target, 2, 2))
                })
            }

            codon.all.new <- mapply(function(x, y, z) {
                if (any(length(y) == 0, length(x) == 0)) {
                    return(z)
                } else
                    if (length(x) <= length(y)) {
                        id.new <- as.numeric(sample(as.character(y), length(x)))
                        id.new.other <-
                            which(!(seq_along(z) %in% id.new))
                        codon.select <- z[x]
                        codon.remainder <- z[-x]
                        codon.all.new <- z
                        codon.all.new[id.new] <-
                            sample(codon.select, 1)
                        codon.all.new[id.new.other] <-
                            sample(codon.remainder, 1)
                        return(codon.all.new)
                    } else {
                        id.new <- y
                        id.new.other <-
                            which(!(seq_along(z) %in% id.new))
                        id.choose <-
                            as.numeric(sample(as.character(x), length(y)))
                        codon.select <- z[id.choose]
                        codon.remainder <- z[-id.choose]
                        codon.all.new <- z
                        codon.all.new[id.new] <-
                            sample(codon.select, 1)
                        codon.all.new[id.new.other] <-
                            sample(codon.remainder, 1)
                        return(codon.all.new)
                    }
            },
                id.target,
                id.good,
                codon.all,
                SIMPLIFY = FALSE)

            seq.mut <- mapply(function(x, y, z) {
                z[y] <- x
                return(z)
            },
                codon.all.new,
                pos.all,
                seq.mut,
                SIMPLIFY = FALSE)
        }
    }
    return(seq.mut)
}

#2.2
dinu_to.keep.region <-
    function(codon.list.alt,
        seq,
        region,
        max.dinu,
        min.dinu) {
        check.max <- !is.na(max.dinu)
        if (check.max) {
            dinu.target <- max.dinu
        } else {
            dinu.target <- min.dinu
        }
        seq.aa <- lapply(seq, function(x) {
            seqinr::translate(s2c(c2s(x)))
        })
        seq.mut <- seq

        #rearrange all synonymous codons in aa order
        for (i in seq_along(codon.list.alt)) {
            aa.tmp <- names(codon.list.alt)[i]
            codon.tmp <- codon.list.alt[[i]]
            #check whether need to rearrange
            codon.target <-
                codon.tmp[substr(codon.tmp, 3, 3) == substr(dinu.target, 1, 1)]
            check.need <- length(codon.target) > 0
            if (check.need) {
                #specify region
                pos.all <- mapply(function(x, y) {
                    which(x == aa.tmp & y)
                }, seq.aa, region, SIMPLIFY = FALSE)
                codon.all <- mapply(function(x, y) {
                    x[y]
                }, seq, pos.all, SIMPLIFY = FALSE)
                id.target <- lapply(codon.all, function(x) {
                    which(x %in% codon.target)
                })

                if (check.max) {
                    id.good <- mapply(function(x, y) {
                        #max
                        codon2nd <- x[y + 1]
                        which(substr(codon2nd, 1, 1) == substr(dinu.target,
                            2, 2))
                    }, seq, pos.all, SIMPLIFY = FALSE)
                } else {
                    id.good <- mapply(function(x, y) {
                        #min
                        codon2nd <- x[y + 1]
                        #note this difference between min and max
                        which(substr(codon2nd, 1, 1) != substr(dinu.target,
                            2, 2))
                    }, seq, pos.all, SIMPLIFY = FALSE)
                }

                codon.all.new <- mapply(function(x, y, z) {
                    if (any(length(y) == 0, length(x) == 0)) {
                        return(z)
                    } else
                        if (length(x) <= length(y)) {
                            id.new <- as.numeric(sample(as.character(y),
                                length(x)))
                            id.new.other <-
                                which(!(seq_along(z) %in% id.new))
                            codon.select <- z[x]
                            codon.remainder <- z[-x]
                            codon.all.new <- z
                            codon.all.new[id.new] <-
                                sample(codon.select, 1)
                            codon.all.new[id.new.other] <-
                                sample(codon.remainder, 1)
                            return(codon.all.new)
                        } else {
                            id.new <- y
                            id.new.other <-
                                which(!(seq_along(z) %in% id.new))
                            id.choose <-
                                as.numeric(sample(as.character(x), length(y)))
                            codon.select <- z[id.choose]
                            codon.remainder <- z[-id.choose]
                            codon.all.new <- z
                            codon.all.new[id.new] <-
                                sample(codon.select, 1)
                            codon.all.new[id.new.other] <-
                                sample(codon.remainder, 1)
                            return(codon.all.new)
                        }
                },
                    id.target,
                    id.good,
                    codon.all,
                    SIMPLIFY = FALSE)

                seq.mut <- mapply(function(x, y, z) {
                    z[y] <- x
                    return(z)
                },
                    codon.all.new,
                    pos.all,
                    seq.mut,
                    SIMPLIFY = FALSE)
            }
        }

        return(seq.mut)
    }
