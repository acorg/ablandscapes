
#include <RcppArmadillo.h>
#include "titer_likelihood.h"

//' Calculate the negative log likelihood of a linear model
//' 
//' This is the base function used by the optimizer to calculate the negative
//' log likelihood of a given set of linear model parameters
//' 
//' @param par A vector of parameters, intercept followed by coefficients for each
//'   coordinate dimension, or in the case of getting likelihood for a given height 
//'   simply the coefficients for each coordinate dimension
//' @param max_titers Numeric vector of the upper bounds of the measured titers
//' @param min_titers Numeric vector of the lower bounds of the measured titers
//' @param ag_coords Matrix of antigenic coordinates relative to the landscape
//'   coordinates being modelled
//' @param ag_weights A vector of weights to apply to each antigen, according to
//'   their distance from the point being modelled
//' @param error_sd The expected standard deviation of titer error
//' 
//' @name negll_titer_lm
//' 

//' @rdname negll_titer_lm
// [[Rcpp::export]]
double negll_titer_lm(
    const arma::vec &par,
    const arma::mat &max_titers,
    const arma::mat &min_titers,
    const arma::mat &ag_coords,
    const arma::vec &ag_weights,
    const double &error_sd
  ) {
  
  double overall_prob = 0;
  double predicted_titer;
  double log_prediction_prob;
  
  for(arma::uword i = 0; i < ag_coords.n_rows; ++i) {
    if(ag_weights(i) > 0) {
      predicted_titer = par(0);
      for(arma::uword j = 0; j < ag_coords.n_cols; ++j) {
        predicted_titer += par(j+1) * ag_coords(i,j);
      }
      for(arma::uword x = 0; x < max_titers.n_rows; ++x) {
        if(std::isfinite(max_titers(x,i))) {
          log_prediction_prob = titer_prediction_negll(
            max_titers(x,i), 
            min_titers(x,i), 
            predicted_titer, 
            error_sd
          );
          overall_prob -= log_prediction_prob * ag_weights(i);
        }
      }
    }
  }
  
  return(overall_prob);
  
}


//' @rdname negll_titer_lm
// [[Rcpp::export]]
double negll_lndscp_height(
    arma::vec par,
    const double &lndscp_height,
    const arma::mat &max_titers,
    const arma::mat &min_titers,
    const arma::mat &ag_coords,
    const arma::vec &ag_weights,
    const double &error_sd
  ) {
  
  // Set the intercept as the height provided
  par.insert_rows(0, lndscp_height);
  
  // Pass on the the normal lm function
  return(
    negll_titer_lm(
      par,
      max_titers,
      min_titers,
      ag_coords,
      ag_weights,
      error_sd
    )
  );
  
}
