#' codon_usage
#'
#' Constructing \code{codon_usage_data} from DNAStringSet.
#'
#' @param ... filepath or DNAstringSet
#' @return A codon_usage_data-class object
#' @seealso \code{\link{get_cu}}
#' @examples
#' # Creating a codon_usage class directly from system file
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' cu <- codon_usage(filepath)
#'
#' # Creating from exsisting DNAStringSet object
#' seq <- Biostrings::DNAStringSet("ATCGATCGA")
#' cu <- codon_usage(seq)
#' @export
setGeneric(name = "codon_usage",
  def = function(object, ...){standardGeneric("codon_usage")})

#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @rdname codon_usage
#' @export
setMethod(f = "codon_usage", signature = "character",
  definition = function(object){
    dnaseq <- readDNAStringSet(filepath = object)
    cu <- oligonucleotideFrequency(dnaseq, width = 3, step = 3)
    return(new(
      "codon_usage_data",
      dnaseq = dnaseq,
      cu = cu
    ))
  })

#' @rdname codon_usage
#' @export
setMethod(f = "codon_usage", signature = "DNAStringSet",
  definition = function(object){
    cu <- oligonucleotideFrequency(object, width = 3, step = 3)
    return(new(
      "codon_usage_data",
      dnaseq = object,
      cu = cu
    ))
  })

