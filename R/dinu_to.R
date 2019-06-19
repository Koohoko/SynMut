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
#' @import methods
#' @import Biostrings
#' @name dinu_to
#' @rdname dinu_to-methods
setGeneric(
    name = "dinu_to",
    def = function(object, max.dinu = NA, min.dinu = NA, keep = FALSE,
                   ...) standardGeneric("dinu_to")
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
        check.valid.dinu <- toupper(dinu.input) %in%
            toupper(mkAllStrings(c("a", "t", "c", "g"), 2))
        if (!check.valid.dinu) {
            stop("Please input valid DNA dinucleotide")
        }

        # extract region ------------------------------------------------------

        region <- object@region
        check.region <- all(is.na(region))
        seq <- convert_to_seq(object@dnaseq)
        seq.region <- extract_region(object, check.region)

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

        # merge region --------------------------------------------------------

        if (!check.region) {
            seq.mut <- mapply(function(x, y, z) {
                x[y] <- z
                return(x)
            }, seq, region, seq.mut, SIMPLIFY = FALSE)
        }

        # synonmous codon interchange -----------------------------------------

        seq.mut <-
            dinu_to.keep(check.region,
                seq.mut,
                max.dinu,
                min.dinu,
                codon.list.alt,
                region)

        # export --------------------------------------------------------------

        seq.mut <- DNAStringSet(unlist(lapply(seq.mut, c2s)))
        return(new(
            "regioned_dna",
            dnaseq = seq.mut,
            region = object@region
        ))

    }
)

###########################################################################
# internal helper function ------------------------------------------------
###########################################################################

#1: not keep codon usage, mutate at all 12, 23 and 31 regions

dinu_to.not.keep <-
    function(seq.region,
        codon.list.alt,
        max.dinu,
        min.dinu) {
        ## get optimal codon
        codon.list.optm <- get_optimal_codon(codon.list.alt, max.dinu, min.dinu)

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

get_optimal_codon <- function(codon.list.alt, max.dinu, min.dinu){
    check.max <- !is.na(max.dinu)
    target.codon <- c(max.dinu, min.dinu)[!is.na(c(max.dinu, min.dinu))]
    codon.list.optm <- lapply(codon.list.alt, function(x) {
        filter1 <- if(check.max){
            x[str_detect(x, pattern = target.codon)]
        } else {
            x[!str_detect(x, pattern = target.codon)]
        }

        if (length(filter1) == 1) {
            return(filter1)
        } else {
            pattern2.1 <- paste0("^", substr(target.codon, 2, 2))
            pattern2.2 <-
                paste0(substr(target.codon, 1, 1), "$")
            if (length(filter1) > 1){
                x <- filter1
            }
            filter2.1 <- if(check.max){
                x[str_detect(x, pattern = pattern2.1)]
            } else {
                x[!str_detect(x, pattern = pattern2.1)]
            }

            filter2.2 <- if(check.max){
                x[str_detect(x, pattern = pattern2.2)]
            } else {
                x[!str_detect(x, pattern = pattern2.2)]
            }

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
    return(codon.list.optm)
}

#2: keep codon usage, mutate at only 31 regions
dinu_to.keep <-
    function(check.region, seq.mut, max.dinu, min.dinu, codon.list.alt,
        region) {
        check.max <- !is.na(max.dinu)
        if (check.max) {
            dinu.target <- max.dinu
        } else {
            dinu.target <- min.dinu
        }
        seq.aa <- lapply(seq.mut, function(x) {
            seqinr::translate(s2c(c2s(x)))
        })

        check.head.g <- lapply(seq.mut, function(x){
            substr(x, 1, 1) == substr(dinu.target, 2, 2)
        })

        if(!check.region){ # have region
            pos.to.mut <- mapply(function(x, y){
                intersect((which(x) - 1), which(y))
            }, check.head.g, region)
        } else { # no region
            pos.to.mut <- lapply(check.head.g, function(x){
                setdiff((which(x) - 1), 0)
            })
        }

        seq.mut <- mapply(function(seq.mut.i, pos.to.mut.i, seq.aa.i, region.i){
            if(length(pos.to.mut.i) < 1){
                return(seq.mut.i)
            } else {
                aa.to.mut <- seq.aa.i[pos.to.mut.i]
                aa.to.mut.unique <- unique(aa.to.mut)

                for(i in seq_along(aa.to.mut.unique)){
                    aa.tmp <- aa.to.mut.unique[i]
                    pos.to.mut.tmp <- pos.to.mut.i[aa.to.mut == aa.tmp]

                    if(!check.region){ # have region
                        pos.for.select <- intersect(which(seq.aa.i == aa.tmp),
                            which(region.i))
                    } else { # no region
                        pos.for.select <- which(seq.aa.i == aa.tmp)
                    }
                    pos.for.select <- setdiff(pos.for.select, pos.to.mut.tmp)

                    if(check.max){ # max dinu
                        pos.to.mut.tmp <-
                            pos.to.mut.tmp[substr(seq.mut.i[pos.to.mut.tmp],
                                3, 3) != substr(dinu.target, 1, 1)]
                        pos.for.select <-
                            pos.for.select[substr(seq.mut.i[pos.for.select],
                                3, 3) == substr(dinu.target, 1, 1)]
                        if(length(pos.to.mut.tmp) > 0){
                            if(length(pos.to.mut.tmp)>=length(pos.for.select)){
                                pos.to.mut.tmp <-as.numeric(sample(as.character(
                                    pos.to.mut.tmp), length(pos.for.select)))
                            } else {
                                pos.for.select <-as.numeric(sample(as.character(
                                    pos.for.select), length(pos.to.mut.tmp)))
                            }
                            codon.to.mut <- seq.mut.i[pos.to.mut.tmp]
                            codon.for.select <- seq.mut.i[pos.for.select]
                            seq.mut.i[pos.to.mut.tmp] <- codon.for.select
                            seq.mut.i[pos.for.select] <- codon.to.mut
                        }

                    } else { # min dinu
                        pos.to.mut.tmp <-
                            pos.to.mut.tmp[substr(seq.mut.i[pos.to.mut.tmp],
                                3, 3) == substr(dinu.target, 1, 1)]
                        pos.for.select <-
                            pos.for.select[substr(seq.mut.i[pos.for.select],
                                3, 3) != substr(dinu.target, 1, 1)]
                        if(length(pos.to.mut.tmp) > 0){
                            if(length(pos.to.mut.tmp)>=length(pos.for.select)){
                                pos.to.mut.tmp <-as.numeric(sample(as.character(
                                    pos.to.mut.tmp), length(pos.for.select)))
                            } else {
                                pos.for.select <-as.numeric(sample(as.character(
                                    pos.for.select), length(pos.to.mut.tmp)))
                            }
                            codon.to.mut <- seq.mut.i[pos.to.mut.tmp]
                            codon.for.select <- seq.mut.i[pos.for.select]
                            seq.mut.i[pos.to.mut.tmp] <- codon.for.select
                            seq.mut.i[pos.for.select] <- codon.to.mut
                        }
                    }
                }
                return(seq.mut.i)
            }
        }, seq.mut, pos.to.mut, seq.aa, region)

        return(seq.mut)
    }
