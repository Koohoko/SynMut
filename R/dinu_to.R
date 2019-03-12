#' Maximize or minimize the usage of certain dinucleotide.
#'
#' Input string of a dinucleotide to either the "max.dinu = " or "min.codon = "
#' parameter to maximize or minimize the usage of certain codon in the sequence.
#' Using a greedy algorithm with priority given to dinucleotide12 or
#' dinucleotide23.
#'
#' @param object A regioned_dna object.
#' @param ...
#' @param max.dinu A string of a dinucleotide.
#' @param min.dinu A string of a dinucleotide.
#'
#' @return regioned_dna
#' @export
#' @seealso \code{\link{input_seq}}, \code{\link{codon_to}},
#'   \code{\link{codon_random}}, \code{\link{codon_mimic}},
#'   \code{\link{codon_mimic_dinu_to}}
#' @examples
#' get_du(dinu_to(rgd.seq, max = "cg")) - get_du(rgd.seq)
#' get_du(dinu_to(rgd.seq, min = "aa")) - get_du(rgd.seq)
#' @include regioned_dna_Class.R
#' @name dinu_to
#' @rdname dinu_to-method
setGeneric(name = "dinu_to",
  def = function(object, max.dinu = NA, min.dinu = NA, ...){
    standardGeneric("dinu_to")
  })

#' @name dinu_to
#' @rdname dinu_to-method
setMethod(f = "dinu_to", signature = "regioned_dna",
  definition = function(object, max.dinu, min.dinu){
    max.dinu <- toupper(max.dinu)
    min.dinu <- toupper(min.dinu)
    dinu.input <- c(max.dinu, min.dinu)
    dinu.input <- dinu.input[!sapply(dinu.input, is.na)]
    if(length(dinu.input) == 0){
      stop("please specify at least one parameter of 'max.dinu' or 'min.dinu'")
    }
    if(all(!is.na(c(max.dinu, min.dinu)))){
      stop("Only one parameter is supported between 'max.dinu' and 'min.dinu'")
    }
    if(length(dinu.input) > 2){
      stop("Only one dinucleotide is allowed for input in 'max.dinu' or 'min.dinu'")
    }
    check.valid.dinu <- toupper(dinu.input) %in% toupper(seqinr::words(2))
    if(!check.valid.dinu){
      stop("Please input valid DNA dinucleotide")
    }

    # extract region ----------------------------------------------------------

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

    # mutation ----------------------------------------------------------------

    codon.list <- sapply(seqinr::ucoweight(''), function(x){toupper(names(x))})
    codon.list.alt <- codon.list[sapply(codon.list, function(x){length(x)>1})]

    ## get optimal dinucleotide -----------------------------------------------

    if(!is.na(max.dinu)){ #max.dinu
      codon.list.optm <- lapply(codon.list.alt, function(x){
        filter1 <- x[stringr::str_detect(x, pattern = max.dinu)]
        if(length(filter1) == 1){
          return(filter1)
        } else if(length(filter1) < 1){
          filter2.1 <- x[stringr::str_detect(x,
            pattern = paste0("^", substr(max.dinu, 2, 2)))]
          filter2.2 <- x[stringr::str_detect(x,
            pattern = paste0(substr(max.dinu, 1, 1), "$"))]
          filter2 <- c(filter2.1, filter2.2)
          if(length(filter2)<1){
            return(sample(x, 1))
          } else{
            return(names(sort(table(filter2), decreasing = T))[1])
          }
        } else {
          filter2.1 <- filter1[stringr::str_detect(filter1,
            pattern = paste0("^", substr(max.dinu, 2, 2)))]
          filter2.2 <- filter1[stringr::str_detect(filter1,
            pattern = paste0(substr(max.dinu, 1, 1), "$"))]
          return(names(sort(table(c(filter2.1, filter2.2)),decreasing = T))[1])
        }
      })
    } else { #min.dinu
      codon.list.optm <- lapply(codon.list.alt, function(x){
        filter1 <- x[!stringr::str_detect(x, pattern = min.dinu)]
        if(length(filter1) == 1){
          return(filter1)
        } else if(length(filter1) < 1){
          filter2.1 <- x[!stringr::str_detect(x,
            pattern = paste0("^", substr(min.dinu, 2, 2)))]
          filter2.2 <- x[!stringr::str_detect(x,
            pattern = paste0(substr(min.dinu, 1, 1), "$"))]
          filter2 <- c(filter2.1, filter2.2)
          if(length(filter2)<1){
            return(sample(x, 1))
          } else{
            return(names(sort(table(filter2), decreasing = T))[1])
          }
        } else {
          filter2.1 <- filter1[!stringr::str_detect(filter1,
            pattern = paste0("^", substr(min.dinu, 2, 2)))]
          filter2.2 <- filter1[!stringr::str_detect(filter1,
            pattern = paste0(substr(min.dinu, 1, 1), "$"))]
          return(names(sort(table(c(filter2.1, filter2.2)),decreasing = T))[1])
        }
      })
    }

    ## mutate codons to optimals ----------------------------------------------

    seq.mut <- seq.region

    for(i in seq_along(codon.list.optm)){
      codons.alt <- codon.list.alt[[i]]
      codons.alt <- codons.alt[codons.alt != codon.list.optm[[i]]]
      seq.mut <- sapply(seq.mut, function(x){
        x[x %in% codons.alt] <- codon.list.optm[[i]]
        return(x)
      })
    }

    # merge region ------------------------------------------------------------

    if(!check.region){
      seq.mut <- mapply(function(x, y, z){
        x[y] <- z
        return(x)
      }, seq, object@region, seq.mut)
    }
    seq.mut <- Biostrings::DNAStringSet(sapply(seq.mut, c2s))
    return(input_seq(seq.mut, region = object@region))

  })