#' 
#' 
#' @title rDNA - generate random DNA sequences
#' 
#' @description generate random DNA sequences
#' @name rDNA
#' 
#' 
#' @param length the length, in nucleotide positions, of the sequence to simulate
#' @param prob a vector of weights for the four possible nucleotides (A, C, G, and T)
#' 
#' 
#' @details 
#' Uses the function \code{base::sample()} to generate random sequences of nucleotides.
#' 
#' 
#' @return 
#' A matrix of nucleotides consisting of one row (sample) that should be compatible with \code{ape::as.DNAbin()}.
#' 
#' @seealso 
#' \code{base::sample()}
#' \code{ape::as.DNAbin()}
#' 
#' @examples 
#' myDNA <- rDNA()
#' 
#' \dontrun{
#' ape::as.DNAbin(myDNA)
#' }
#' 
#' @export
rDNA <- function(length = 100, prob = c(0.25, 0.25, 0.25, 0.25)){
  myDNA <- matrix(sample(c('A', 'C', 'G', 'T'), size = length, replace = TRUE, prob = prob), nrow = 1)
  return(myDNA)
}

