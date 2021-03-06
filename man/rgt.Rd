% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rgt}
\alias{rgt}
\title{rgt - random genotypes}
\usage{
rgt(nsamp = 4L, nvar = 3L, pphased = as.numeric(c(0.5)),
  pploid = as.numeric(c(0, 1)), pallele = as.numeric(c(0.5, 0.5)),
  verbose = 0L)
}
\arguments{
\item{nsamp}{number of samples (columns) to simulate}

\item{nvar}{number of variants (rows) to simulate}

\item{pphased}{probability each genotype will be phased}

\item{pploid}{probability each genotype will be a ploidy level}

\item{pallele}{probability of each allele}

\item{verbose}{should verbose output be produced (1) or not (0)}
}
\value{
A matrix with variants (loci) in rows and samples in columns.
}
\description{
Generate a matrix of random genotypes.
}
\details{
The function \code{rgt} generates a matrix of random genotypes.
This matrix consists of loci or variants in rows and samples in columns where the number of loci and samples is arbitrary.
Genotypes can be delimited by either the forward slash (/), for unphased genotypes, or the pipe (|), for phased genotypes.
The matrix of genotypes generated by this function is intended to be compatible with a matrix of genotypes extracted from variant cal format data (VCF).
The alleles are encoded in a zero-based numeric manner consistent with the VCF specification.
Here, 0 is the first (typically the reference) allele, 1 is the first alternate allele, 2 is the second alternate allele, and so on.
While these genotypes are intended to conform to the VCF specification it is hoped they will be generalizable to broader applications.


The parameter \strong{nvar} determines the number of loci or variants (rows) to simulate.
This integer value is arbitrary to accomodate a range of locus counts.


The parameter \strong{nsamp} determines the number of samples (columns) to simulate.
This integer value is arbitrary to accomodate a range of sample sizes.


The parameter \strong{pphased} indicates the probability that a genotype will be phased.
This parameter consists of a vector containing one element ranging in value from zero to one and describes the probability an individual genotype will be phased.
For example, \code{rgt(pphased=c(0))} will produce genotypes that have a probability of zero for being phased.
The code \code{rgt(pphased=c(1))} will produce genotypes that have a probability of one for being phased.
The code \code{rgt(pphased=c(0.2))} will produce genotypes that have a probability of 0.2 that they will be phased.
The probability for each genotype returned is calculated independently.


The parameter \strong{pploid} indicates the probability a genotype will be of a particular ploidy or copy number.
This is an n-element vector where each element indicates the probability a genotype will be n-ploid.
For example, \code{rgt(pploid=c(1))} will produce genotypes that have a probability of 1 for being one-ploid (haploid).
The code \code{rgt(pploid=c(0,1))} will produce genotypes that have a probability of 0 for being one-ploid and a probability of 1 for being two-ploid (diploid).
The code \code{rgt(pploid=c(0,0,1))} will produce genotypes that have a probability of 0 for being one-ploid, a probability of 0 for being two-ploid (diploid), and a probability of 1 for being three-ploid (triploid).
The code \code{rgt(pploid=c(0.1,0.5,0.4))} will produce genotypes that have a probability of 0.1 for being one-ploid, 0.5 for being two-ploid, and 0.4 for being three-ploid.
The values in this vector are of arbitrary length but should sum to one.

The parameter \strong{pallele} indicates the probability of each allele to occur in a genotype.
This is an n-element vector where each element indicates the probability of the nth allele occuring.
Here the first element is the probability of the 0th allele.
For example, \code{rgt(pallele = c(1))} will produce genotypes where the 0th allele (0) will occur with a probability of 1.
The code \code{rgt(pallele = c(0,1))} will produce genotypes where the 0th allele (0) will occur with a probability of 0 and the 1st allele will occur with a probability of 1.
The code \code{rgt(pallele = c(0,0,1))} will produce genotypes where the 0th allele (0) will occur with a probability of 0, the 1st allele will occur with a probability of 0, and the 2nd allele will occur with a probability of 1.
The code \code{rgt(pallele = c(0.2,0.7,0.1))} will produce genotypes where the 0th allele (0) will occur with a probability of 0.2, the 1st allele will occur with a probability of 0.7, and the 2nd allele will occur with a probability of 0.1.
The values in this vector are of arbitrary length but should sum to one.
}
\examples{
rgt()


}
