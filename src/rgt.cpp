#include <Rcpp.h>
#include <sstream>
#include <stdlib.h>
//using namespace Rcpp;


std::string simGT(Rcpp::NumericVector pphased,
                  Rcpp::NumericVector pploid,
                  Rcpp::NumericVector pallele,
                  int verbose = 0
){
  int i = 0;
  int j = 0;
  int ploid = 0;
  Rcpp::NumericVector myRand = Rcpp::runif(1);
  
//=int verbose = 0; // FALSE
//int verbose = 1; // TRUE
  
  if(verbose == 1){
    Rcpp::Rcout << "pploid: " << pploid(0);
    for(i=1; i<pploid.size(); i++){
      Rcpp::Rcout << ", " << pploid(i);
    }
    Rcpp::Rcout << "\n";
  }

  
  // Convert ploidy probabilities to thresholds
  Rcpp::NumericVector ploid_threshold(pploid.size() + 1);
  for(i=1; i<ploid_threshold.size(); i++){
    ploid_threshold(i) = ploid_threshold(i-1) + pploid(i - 1);
  }
  
  if(verbose == 1){
    Rcpp::Rcout << "ploid_threshold: " << ploid_threshold(0);
    for(i=1; i<ploid_threshold.size(); i++){
      Rcpp::Rcout << ", " << ploid_threshold(i);
    }
    Rcpp::Rcout << "\n";
  }
  
  // Determine ploid
  for(i=0; i<ploid_threshold.size() - 1; i++){
    if(myRand(0) > ploid_threshold(i) & myRand(0) <= ploid_threshold(i + 1) ){
      ploid = i;
    }
  }
  
  if(verbose == 1){
    Rcpp::Rcout << "Locus has " << ploid + 1 << " copies\n";
  }
  
  // Determine phasing
  std::string delim = "/";
  myRand = Rcpp::runif(1);
  if(myRand(0) < pphased(0)){
    delim = "|";
  }
  
  // Create vector of alleles
  Rcpp::IntegerVector myAlleles(ploid + 1);

  if(verbose == 1){
    Rcpp::Rcout << "pallele is " << pallele(0);
    for(i=1; i<pallele.size(); i++){
      Rcpp::Rcout << ", " << pallele(i);
    }
    Rcpp::Rcout << "\n";
  }

  // Convert allele probabilities to thresholds
  Rcpp::NumericVector allele_threshold(pallele.size() + 1);  
  for(i=1; i<allele_threshold.size(); i++){
//    Rcpp::Rcout << i << " Adding " << allele_threshold( i - 1 ) << " and " << pallele( i - 1 ) << "\n";
    allele_threshold(i) = allele_threshold( i - 1 ) + pallele( i - 1 );
  }

  if(verbose == 1){
    Rcpp::Rcout << "Allelic thresholds " << allele_threshold(0);
    for(i=1; i<allele_threshold.size(); i++){
      Rcpp::Rcout << ", " << allele_threshold(i);
    }
    Rcpp::Rcout << "\n";
  }

  // Assign allelic state
  for(i=0; i<=ploid; i++){
    // Allele copy i
    myRand = Rcpp::runif(1);
    if(verbose == 1){
      Rcpp::Rcout << "Random number is: " << myRand(0) << "\n";
    }
    if(verbose == 1){
      Rcpp::Rcout << "  Allele copy " << i << "\n";
    }
    for(j=0; j<allele_threshold.size() - 1; j++){
      // Determine state of allele i
      if(verbose == 1){
        Rcpp::Rcout << "    Allelic state " << j << ", allele_threshold: " << allele_threshold(j) << " - " << allele_threshold(j + 1) << "\n";
      }
      if( myRand(0) > allele_threshold(j) & myRand(0) <= allele_threshold(j + 1) ){
        myAlleles(i) = j;
        if(verbose == 1){
          Rcpp::Rcout << "      Changed allele: " << j << " to " << myAlleles(i) << "\n";
        }
      }
    }
    if(verbose == 1){
      Rcpp::Rcout << "    Allele called is " << myAlleles(i) << "\n";
    }
  }

  // Concatenate alleles into a string
  std::stringstream sstm;
  sstm << myAlleles(0);
  for(i=1; i<myAlleles.size(); i++){
    sstm << delim << myAlleles(i);
  }

  std::string myGT = sstm.str();
  // std::string myGT = "0/1";
  if(verbose == 1){
    Rcpp::Rcout << "GT is: " << myGT << "\n";
    Rcpp::Rcout << "\n";
  }
  return(myGT);
}



//'
//' @title rgt - random genotypes
//' @description
//' Generate a matrix of random genotypes.
//' @name rgt
//' 
//' @param nsamp number of samples (columns) to simulate
//' @param nvar number of variants (rows) to simulate
//' @param pphased probability each genotype will be phased
//' @param pploid probability each genotype will be a ploidy level
//' @param pallele probability of each allele
//' @param verbose should verbose output be produced (1) or not (0)
//' 
//' 
//' @details
//' The function \code{rgt} generates a matrix of random genotypes.
//' This matrix consists of loci or variants in rows and samples in columns where the number of loci and samples is arbitrary.
//' Genotypes can be delimited by either the forward slash (/), for unphased genotypes, or the pipe (|), for phased genotypes.
//' The matrix of genotypes generated by this function is intended to be compatible with a matrix of genotypes extracted from variant cal format data (VCF).
//' The alleles are encoded in a zero-based numeric manner consistent with the VCF specification.
//' Here, 0 is the first (typically the reference) allele, 1 is the first alternate allele, 2 is the second alternate allele, and so on.
//' While these genotypes are intended to conform to the VCF specification it is hoped they will be generalizable to broader applications.
//' 
//' 
//' The parameter \strong{nvar} determines the number of loci or variants (rows) to simulate.
//' This integer value is arbitrary to accomodate a range of locus counts.
//' 
//' 
//' The parameter \strong{nsamp} determines the number of samples (columns) to simulate.
//' This integer value is arbitrary to accomodate a range of sample sizes.
//' 
//' 
//' The parameter \strong{pphased} indicates the probability that a genotype will be phased.
//' This parameter consists of a vector containing one element ranging in value from zero to one and describes the probability an individual genotype will be phased.
//' For example, \code{rgt(pphased=c(0))} will produce genotypes that have a probability of zero for being phased.
//' The code \code{rgt(pphased=c(1))} will produce genotypes that have a probability of one for being phased.
//' The code \code{rgt(pphased=c(0.2))} will produce genotypes that have a probability of 0.2 that they will be phased.
//' The probability for each genotype returned is calculated independently.
//' 
//' 
//' The parameter \strong{pploid} indicates the probability a genotype will be of a particular ploidy or copy number.
//' This is an n-element vector where each element indicates the probability a genotype will be n-ploid.
//' For example, \code{rgt(pploid=c(1))} will produce genotypes that have a probability of 1 for being one-ploid (haploid).
//' The code \code{rgt(pploid=c(0,1))} will produce genotypes that have a probability of 0 for being one-ploid and a probability of 1 for being two-ploid (diploid).
//' The code \code{rgt(pploid=c(0,0,1))} will produce genotypes that have a probability of 0 for being one-ploid, a probability of 0 for being two-ploid (diploid), and a probability of 1 for being three-ploid (triploid).
//' The code \code{rgt(pploid=c(0.1,0.5,0.4))} will produce genotypes that have a probability of 0.1 for being one-ploid, 0.5 for being two-ploid, and 0.4 for being three-ploid.
//' The values in this vector are of arbitrary length but should sum to one.
//' 
//' The parameter \strong{pallele} indicates the probability of each allele to occur in a genotype.
//' This is an n-element vector where each element indicates the probability of the nth allele occuring.
//' Here the first element is the probability of the 0th allele.
//' For example, \code{rgt(pallele = c(1))} will produce genotypes where the 0th allele (0) will occur with a probability of 1.
//' The code \code{rgt(pallele = c(0,1))} will produce genotypes where the 0th allele (0) will occur with a probability of 0 and the 1st allele will occur with a probability of 1.
//' The code \code{rgt(pallele = c(0,0,1))} will produce genotypes where the 0th allele (0) will occur with a probability of 0, the 1st allele will occur with a probability of 0, and the 2nd allele will occur with a probability of 1.
//' The code \code{rgt(pallele = c(0.2,0.7,0.1))} will produce genotypes where the 0th allele (0) will occur with a probability of 0.2, the 1st allele will occur with a probability of 0.7, and the 2nd allele will occur with a probability of 0.1.
//' The values in this vector are of arbitrary length but should sum to one.
//' 
//' 
//' @return
//' A matrix with variants (loci) in rows and samples in columns.
//' 
//' 
//' @examples
//' rgt()
//' 
//'
//' @export
// [[Rcpp::export]]
Rcpp::CharacterMatrix rgt(
    int nsamp = 4,
    int nvar = 3,
    Rcpp::NumericVector pphased = Rcpp::NumericVector::create(0.5),
    Rcpp::NumericVector pploid = Rcpp::NumericVector::create(0,1),
    Rcpp::NumericVector pallele = Rcpp::NumericVector::create(0.5,0.5),
    int verbose = 0
){
  // Initialize return matrix.
  Rcpp::CharacterMatrix gtMat(nvar, nsamp);
  Rcpp::CharacterVector myCols(nsamp);
  Rcpp::CharacterVector myRows(nvar);

  int i = 0;
  int j = 0;
  Rcpp::NumericVector mySum(1);
  
  // Check 0 <= pphased <= 1.
  if( pphased[0] < 0 | pphased[0] > 1 ){
    Rcpp::Rcerr << "pphased: " << pphased[0] << "\n";
    Rcpp::Rcerr << "pphased must be between 0 and 1.\n";
    return( Rcpp::CharacterMatrix(1) );
  }
  
  
  // Check pploid sum == 1.
  mySum[0] = 0;
  for(i=0; i<pploid.size(); i++){
    mySum[0] = mySum[0] + pploid(i);
  }

  if( mySum[0] < 0.999 | mySum[0] > 1.001 ){
    Rcpp::Rcerr << "pploid deos not sum to one.\n";
    Rcpp::Rcerr << pploid(0);
    for(i=1;i<pploid.size();i++){
      Rcpp::Rcerr << ", " << pploid(i);
    }
    Rcpp::Rcerr << " == " << mySum[0];
    Rcpp::Rcerr << "\n";
    return( Rcpp::CharacterMatrix(1) );
  }

  
  // Check pallele sum == 1.
  mySum[0] = 0;
  for(i=0; i<pallele.size(); i++){
    mySum[0] = mySum[0] + pallele(i);
  }

  if( mySum[0] < 0.999 | mySum[0] > 1.001 ){
    Rcpp::Rcerr << "pallele deos not sum to one.\n";
    Rcpp::Rcerr << pallele(0);
    for(i=1;i<pallele.size();i++){
      Rcpp::Rcerr << ", " << pallele(i);
    }
    Rcpp::Rcerr << " == " << mySum[0];
    Rcpp::Rcerr << "\n";
    return( Rcpp::CharacterMatrix(1) );
  }

  
  // Loop through matrix
  for(i=0; i<nsamp;i++){
    std::stringstream sstm;
    sstm << "Sample" << i;
    myCols(i) = sstm.str();
    
    for(j=0; j<nvar; j++){
      std::stringstream sstm;
      sstm << "Var" << j;
      myRows(j) = sstm.str();
      
      // Simulate a single genotype
      gtMat(j,i) = simGT(pphased, pploid, pallele, verbose);
    }
  }
  
  // Add col and row names
  colnames(gtMat) = myCols;
  rownames(gtMat) = myRows;
  return(gtMat);
}

