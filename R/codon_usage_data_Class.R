#' \code{codon_usage_data}: An S4 class to record codon usage of a sequence
#'
#' Recording codon usage matrix.
#'
#' @slot cu a matrix recording the codon usage
#' @slot DNAStringSet a DNAStingSet object recording the sequence(s)
#'
#' @importFrom methods setClass
#' @seealso \code{\link{codon_usage}}, \code{\link{get_cu}}
#' @seealso cu
#' @export
#' @author Haogao Gu
setClass(Class = "codon_usage_data",
  slots = c(
    dnaseq = "DNAStringSet",
    cu = "matrix"
    )
  )

setValidity("codon_usage_data",
  function(object){
  check <- nrow(object@cu) == length(object@dnaseq)
  if(check){TRUE} else {stop("dnaseq and cu must have the same length.")}
  })

##viewer
setMethod(f = "show", signature = "codon_usage_data",
  definition = function(object){
    cat("An object of class ", class(object), "\n", sep = "")
    cat("Number of sequences: ", nrow(object@cu), "\n", sep = "")
  })

##Accessor
#' get_cu
#'
#' access the codon usage matrix
#'
#' @param cu codon_usage_data
#' @return matrix
#' @seealso \code{\link{codon_usage}}
#' @examples
#' codon.matrix <- get_cu(cu)
#' @export
setGeneric(name = "get_cu",
  def = function(object, ...){standardGeneric("get_cu")})

#' @rdname get_cu
setMethod(f = "get_cu", signature = "codon_usage_data",
  definition = function(object){
    object@cu
  })



