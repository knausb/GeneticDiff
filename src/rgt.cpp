#include <Rcpp.h>
#include <sstream>
#include <stdlib.h>
//using namespace Rcpp;


std::string simGT(Rcpp::NumericVector pphased,
                  Rcpp::NumericVector pploid,
                  Rcpp::NumericVector pallele
){
  int i = 0;
  int j = 0;
  int ploid = 0;
  Rcpp::NumericVector myRand = Rcpp::runif(1);
  
//
int verbose = 0; // FALSE
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
//' 
//' 
//' @details
//' The function \code{rgt} generates a matrix of random genotypes.
//' The parameter `nsamp` determines the number of samples (columns) to simulate.
//' The parameter `nvar` determines the number of variants (rows) to simulate.
//' 
//' @return
//' A matrix with samples in columns and variants in rows.
//' 
//'
//' @export
// [[Rcpp::export]]
Rcpp::CharacterMatrix rgt(
    int nsamp = 4,
    int nvar = 3,
    Rcpp::NumericVector pphased = Rcpp::NumericVector::create(0.5),
    Rcpp::NumericVector pploid = Rcpp::NumericVector::create(0,1),
    Rcpp::NumericVector pallele = Rcpp::NumericVector::create(0.5,0.5)
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
      gtMat(j,i) = simGT(pphased, pploid, pallele);
    }
  }
  
  // Add col and row names
  colnames(gtMat) = myCols;
  rownames(gtMat) = myRows;
  return(gtMat);
}

