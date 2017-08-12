#include <Rcpp.h>
// using namespace Rcpp;


//' @title genotype_split - split a genotype into alleles
//' @description
//' Split a genotype into alleles based on the delimiters forward slash (/)  and pipe (|).
//' @name genotype_split
//' 
//' @param myGT a string containing a genotype
//' 
//' @details
//' The funciton \code{genotype_split} splits a genotype into a vector of alleles.
//' The alleles are expected to be delimited by either a forward slash ('/') or a pipe ('|').
//' This format is expected to conform with the VCF specification but may handle other formats as well.
//' 
//' 
//' @examples
//' genotype_split("0/1/1|2")
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::StringVector genotype_split(std::string myGT) {
//Rcpp::StringVector genotype_split(Rcpp::StringVector myGT) {
  
  int i = 0;
  int myStart = 0;
  
  int nalleles = 0;
  // If first position is not a delimeter, increment.
  if(myGT[0] != '/' | myGT[0] != '|'){
    nalleles = nalleles + 1;
  }
  
  // Determine the number of alleles.
  for(i=0; i<myGT.size(); i++){
//    Rcpp::Rcout << myGT[i] << "\n";
    if(myGT[i] == '/' | myGT[i] == '|'){
      nalleles = nalleles + 1;
    }
  }
//  Rcpp::Rcout << "nalleles: " << nalleles << "\n";
  
  // Initialize return.
  Rcpp::StringVector myAlleles(nalleles);
  nalleles = 0;

  // Parse genotype into alleles.
  for(i=0; i<myGT.size(); i++){
    if(myGT[i] == '/' | myGT[i] == '|'){
//      Rcpp::Rcout << "  nalleles: " << nalleles << ", myStart: " << myStart << ", i: " << i << "\n";
      myAlleles[nalleles] = myGT.substr(myStart, i - myStart);
      nalleles = nalleles + 1;
      myStart = i + 1;
    }
  }
  // Last allele.
  myAlleles[nalleles] = myGT.substr(myStart, i - myStart);
  
  return myAlleles;
}


