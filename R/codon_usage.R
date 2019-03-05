#' codon_usage
#'
#' Constructing \code{codon_usage_data} from DNAStringSet.
#'
#' @param \code{object} filepath or DNAstringSet
#' @param \code{region} a data-frame specifying paticular regions that is
#'   allowed to be mutated in the sequences. Both \code{0 / 1} or
#'   \code{TRUE / FALSE} encoding is OK. Please refer to Examples below for
#'   reference.
#' @return A codon_usage_data-class object
#' @seealso \code{\link{get_cu}}, \code{\link{get_region}}
#' @examples
#' # Creating a codon_usage class directly from system file
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' cu <- codon_usage(filepath)
#'
#' #Optionally input with region dataframe
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' filepath.csv <- system.file("extdata", "target_regions.csv", package = "SynMut")
#' region <- read.csv(filepath.csv)
#' cu <- codon_usage(filepath, region)
#'
#' # Creating from exsisting DNAStringSet object
#' seq <- Biostrings::DNAStringSet("ATCGATCGA")
#' cu <- codon_usage(seq)
#' @export
#' @name codon_usage
#' @rdname codon_usage-method
setGeneric(name = "codon_usage",
  def = function(object, region = NA, ...){standardGeneric("codon_usage")})

#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @name codon_usage
#' @rdname codon_usage-method
#' @export
setMethod(f = "codon_usage", signature = signature(object = "character"),
  definition = function(object, region){
    if(all(is.na(region))){
      dnaseq <- readDNAStringSet(filepath = object)
      cu <- oligonucleotideFrequency(dnaseq, width = 3, step = 3)
      return(new(
        "codon_usage_data",
        dnaseq = dnaseq,
        cu = cu,
        region = list(NA)
      ))
    } else {
      region <- as.list(region)
      region <- lapply(region, as.logical)
      dnaseq <- readDNAStringSet(filepath = object)
      cu <- oligonucleotideFrequency(dnaseq, width = 3, step = 3)
      return(new(
        "codon_usage_data",
        dnaseq = dnaseq,
        cu = cu,
        region = region
      ))
    }
  })

#' @name codon_usage
#' @rdname codon_usage-method
#' @export
setMethod(f = "codon_usage", signature = "DNAStringSet",
  definition = function(object, region){
    if(all(is.na(region))){
      cu <- oligonucleotideFrequency(object, width = 3, step = 3)
      return(new(
        "codon_usage_data",
        dnaseq = object,
        cu = cu,
        region = list(NA)
      ))
    } else {
      region <- as.list(region)
      region <- lapply(region, as.logical)
      cu <- oligonucleotideFrequency(object, width = 3, step = 3)
      return(new(
        "codon_usage_data",
        dnaseq = object,
        cu = cu,
        region = region
      ))
    }
  })

