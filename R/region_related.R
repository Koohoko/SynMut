# helper function ---------------------------------------------------------

convert_to_seq <- function(dna.seq){
    seq <- lapply(as.character(dna.seq), function(x) {
        seqinr::splitseq(s2c(x))
    })
    return(seq)
}

extract_region <- function(object, check.region){
    if (!check.region) {
        seq <- convert_to_seq(object@dnaseq)
        seq.region <- mapply(function(x, y) {
            return(x[y])
        }, seq, object@region, SIMPLIFY = FALSE)
    } else {
        seq.region <- lapply(as.character(object@dnaseq),
            function(x) {
                splitseq(s2c(x))
            })
    }
    return(seq.region)
}

region_back <- function(seq.mut, check.region, seq, object){
    if (!check.region) {
        seq.mut <- mapply(function(x, y, z) {
            x[y] <- z
            return(x)
        }, seq, object@region, seq.mut, SIMPLIFY = FALSE)
    }
    seq.mut <- Biostrings::DNAStringSet(unlist(lapply(seq.mut, c2s)))
    return(seq.mut)
}
