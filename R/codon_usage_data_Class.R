
# Class definition --------------------------------------------------------

#' codon_usage_data: An S4 class to record codon usage of a sequence
#'
#' Recording codon usage matrix.
#'
#' @slot \code{cu} a matrix recording the codon usage
#' @slot \code{dnaseq} a DNAStingSet object recording the sequence(s)
#' @slot \code{region} a list specifying paticular regions allowed to be mutated
#'
#' @importFrom methods setClass
#' @seealso \code{\link{codon_usage}}, \code{\link{get_cu}}, \code{\link{get_region}}
#' @seealso cu
#' @export
#' @author Haogao Gu
setClass(Class = "codon_usage_data",
  slots = c(
    dnaseq = "DNAStringSet",
    cu = "matrix",
    region = "list"
    )
  )

setValidity("codon_usage_data",
  function(object){
    check.na <- all(is.na(object@region))
    dnaseq.n <- length(object@dnaseq)
    cu.n <- nrow(object@cu)
    region.n <- length(object@region)
    if(check.na){
      check.num <- dnaseq.n == cu.n
    } else {
      check.num <- all(sapply(c(dnaseq.n, cu.n, region.n),
        function(x){
          x == dnaseq.n
        }))
    }
    if(check.num){TRUE} else{"dnaseq, cu and region must have the same length."}
  })

##viewer
setMethod(f = "show", signature = "codon_usage_data",
  definition = function(object){
    cat("An object of class ", class(object), "\n", sep = "")
    cat("Number of sequences: ", nrow(object@cu), "\n", sep = "")
  })


# Accessor ----------------------------------------------------------------

#' get_cu: get codon usage matrix
#'
#' access the codon usage matrix
#'
#' @param \code{cu} codon_usage_data
#' @return matrix
#' @seealso \code{\link{codon_usage}}, \code{\link{get_region}}
#' @examples
#' codon.matrix <- get_cu(cu)
#' @export
#' @name get_cu
#' @rdname get_cu-method
setGeneric(name = "get_cu",
  def = function(object, ...){standardGeneric("get_cu")})

#' @name get_cu
#' @rdname get_cu-method
setMethod(f = "get_cu", signature = "codon_usage_data",
  definition = function(object){
    object@cu
  })

#' get_region: get the variable region
#'
#' access the variable regions
#'
#' @param \code{cu} codon_usage_data
#' @return list
#' @seealso \code{\link{codon_usage}}, \code{\link{get_cu}}
#' @examples
#' codon.matrix <- get_region(cu)
#' @export
#' @name get_region
#' @rdname get_region-method
setGeneric(name = "get_region",
  def = function(object, ...){standardGeneric("get_region")})

#' @name get_region
#' @rdname get_region-method
setMethod(f = "get_region", signature = "codon_usage_data",
  definition = function(object){
    object@region
  })

