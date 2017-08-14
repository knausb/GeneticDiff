#include <Rcpp.h>
//using namespace Rcpp;


//'
//' @title count_alleles
//' 
//' @description
//' Summarize genotypes and alleles.
//' 
//' @name count_alleles
//' 
//' @param GT a matrix of genotypes represented by strings
//' 
//' 
//' @details
//' The function \code{rgt} generates a matrix of random genotypes.
//' 
//' 
//' @return
//' A DataFramw.
//' 
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame count_alleles(Rcpp::StringMatrix GT){
  Rcpp::IntegerVector nAllele(GT.nrow());
  Rcpp::IntegerVector nHet(GT.nrow());
  Rcpp::IntegerVector nHom(GT.nrow());
  Rcpp::IntegerVector allele_counts(GT.nrow());
  
  int i = 0;
  int j = 0;
  Rcpp::StringVector myAlleles;
  
  for(i=0; i<GT.nrow(); i++){
    // Row counter
    for(j=0; j<GT.ncol(); j++){
      // Column counter
//      myAlleles = GeneticDiff::genotype_split(GT[i,j]);
    }
  }
  
  
  
  return Rcpp::DataFrame::create(Rcpp::_["nAllele"] = nAllele,
                                 Rcpp::_["nHet"] = nHet,
                                 Rcpp::_["nHom"] = nHom, 
                                 Rcpp::_["allele_counts"] = allele_counts
                                );
}

