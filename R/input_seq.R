#' Import region / constructing regioned_dna object
#'
#' Constructing \code{regioned_dna} from DNAStringSet. Optionally input a
#' \doce{region} data.frame to define restricted amino-acid region for mutation.
#'
#' @param \code{object} Filepath or DNAstringSet. The input sequences is
#'   suggested to be in open reading frame(ORF).
#' @param \code{region} A data.frame specifying paticular regions (positions in
#'   amino acid sequence) that is allowed to be mutated in the sequences. Both
#'   \code{0 / 1} or \code{TRUE / FALSE} encoding is OK. Please refer to
#'   Examples below for reference.
#' @return A regioned_dna-class object
#' @seealso \code{\link{get_cu}}, \code{\link{get_du}},
#'   \code{\link{get_region}}, \code{\link{get_dnaseq}}
#' @examples
#' # Creating a input_seq class directly from system file
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#'
#' #Optionally input with region dataframe
#' filepath.fasta <- system.file("extdata", "example.fasta", package = "SynMut")
#' filepath.csv <- system.file("extdata", "target_regions.csv", package = "SynMut")
#' region <- read.csv(filepath.csv)
#' rgd.seq <- input_seq(filepath.fasta, region)
#'
#' # Creating from exsisting DNAStringSet object
#' seq <- Biostrings::DNAStringSet("ATCGATCGA")
#' rgd.seq <- input_seq(seq)
#' @export
#' @name input_seq
#' @rdname input_seq-method
setGeneric(name = "input_seq",
  def = function(object, region = NA, ...){standardGeneric("input_seq")})

#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @name input_seq
#' @rdname input_seq-method
#' @export
setMethod(f = "input_seq", signature = signature(object = "character"),
  definition = function(object, region){
    if(all(is.na(region))){
      dnaseq <- readDNAStringSet(filepath = object)
      return(new(
        "regioned_dna",
        dnaseq = dnaseq,
        region = list(NA)
      ))
    } else {
      region <- as.list(region)
      region <- lapply(region, function(x){as.logical(x[!is.na(x)])})
      dnaseq <- readDNAStringSet(filepath = object)
      return(new(
        "regioned_dna",
        dnaseq = dnaseq,
        region = region
      ))
    }
  })

#' @name input_seq
#' @rdname input_seq-method
#' @export
setMethod(f = "input_seq", signature = "DNAStringSet",
  definition = function(object, region){
    if(all(is.na(region))){
      return(new(
        "regioned_dna",
        dnaseq = object,
        region = list(NA)
      ))
    } else {
      region <- as.list(region)
      region <- lapply(region, function(x){as.logical(x[!is.na(x)])})
      return(new(
        "regioned_dna",
        dnaseq = object,
        region = region
      ))
    }
  })

#' @name input_seq
#' @rdname input_seq-method
#' @export
setMethod(f = "input_seq", signature = "DNAString",
  definition = function(object, region){
    if(all(is.na(region))){
      return(new(
        "regioned_dna",
        dnaseq = object,
        region = list(NA)
      ))
    } else {
      region <- as.list(region)
      region <- lapply(region, function(x){as.logical(x[!is.na(x)])})
      return(new(
        "regioned_dna",
        dnaseq = object,
        region = region
      ))
    }
  })