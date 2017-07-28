#include <Rcpp.h>
#include <sstream>
#include <stdlib.h>
//using namespace Rcpp;


std::string simGT(int pphased,
                  Rcpp::IntegerVector pploid,
                  Rcpp::IntegerVector pallele
){
  
  // Determine phasing
  std::string delim = "/";
//  unsigned seed = 1;
//  srand(seed);
//  int myRand = rand() % 100;
//  Rcpp::Rcout << myRand << "\n";
  
//  if(myRand < pphased){
//    delim = "|";
//  }


  return("0/1");
}



//'
//' @title rgt
//' @description
//' Generate a matrix of random genotypes.
//' 
//' @param nsamp number of samples to simulate
//' @param nvar number of variants to simulate
//' @param pphased probability each genotype will be phased
//' @param pploid probability each genotype will be a ploidy level
//' @param pallele probability of each allele
//' 
//' 
//' @details
//' Generate a matrix of random genotypes.
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
    int pphased = 50,
    Rcpp::IntegerVector pploid = Rcpp::IntegerVector::create(0,100),
    Rcpp::IntegerVector pallele = Rcpp::IntegerVector::create(50,50)
){
  // Initialize return matrix.
  Rcpp::CharacterMatrix gtMat(nvar, nsamp);
  Rcpp::CharacterVector myCols(nsamp);
  Rcpp::CharacterVector myRows(nvar);

  int i = 0;
  int j = 0;
  
  // Check pploid sum < 100.
  int mySum = 0;
  for(i=0; i<pploid.size(); i++){
    mySum = mySum + pploid(i);
  }
  if(mySum > 100){
    Rcpp::Rcerr << "pploid sum is greater than 100\n";
    return( Rcpp::CharacterMatrix(1) );
  }
  
  // Check pallele sum < 100.
  mySum = 0;
  for(i=0; i<pallele.size(); i++){
    mySum = mySum + pallele(i);
  }
  if(mySum > 100){
    Rcpp::Rcerr << "pallele sum is greater than 100\n";
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

