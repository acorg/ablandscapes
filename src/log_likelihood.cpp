#include <Rcpp.h>
using namespace Rcpp;

//' Calculate log likelihood of single fitted HI titer
//' 
//' This is the base function for performing a maximum likelihood calculation for a 
//' single fitted HI titer given upper and lower limits of the measured value.
//'
//' @param max_titer The upper bound of the measured titer.
//' @param min_titer The lower bound of the measured titer.
//' @param pred_titer The predicted titer.
//' @param error_sd The standard deviation of the error.
//'
//' @details This function simply calculates to log likelihood of a predicted measurement 
//' given the upper and lower bounds of the measurement, the main assumption being that the 
//' associated error is normally distributed.
//'
//' @return Returns the log-likelihood of the measured titer given the measured titer 
//' bounds and error standard deviation supplied.
// [[Rcpp::export]]
double calc_titer_loglik(double max_titer,
                         double min_titer,
                         double pred_titer,
                         double error_sd)
{
  return(R::logspace_sub(R::pnorm5(max_titer, pred_titer, error_sd,1,1),
                         R::pnorm5(min_titer, pred_titer, error_sd,1,1)));
}


//' Calculate the total negative log-likelihood of a predicted titer set
//' 
//' This is a base function to sum the total _negative_ log likelihood of a predicted 
//' and measured set of HI titers.
//' 
//' @param max_titers Numeric vector of the upper bounds of the measured titers.
//' @param min_titers Numeric vector of the lower bounds of the measured titers.
//' @param pred_titers Numeric vector of predicted titers.
//' @param titer_weights Numeric vector of titer weights.
//' @param hi_error_sd The standard deviation of the error.
//' 
//' @details This function simply calculates to log likelihood of a series of predicted 
//' measurements given the upper and lower bounds of the measurements, the main assumption 
//' being that the associated error is normally distributed.
//' @export
// [[Rcpp::export]]
double calc_titer_set_negll(NumericVector max_titers,
                            NumericVector min_titers,
                            NumericVector pred_titers,
                            NumericVector titer_weights,
                            double hi_error_sd)
{
  
  double total_negll = 0;
  for(int i = 0; i < min_titers.length(); ++i) {
    total_negll -= calc_titer_loglik(max_titers[i],
                                     min_titers[i],
                                     pred_titers[i],
                                     hi_error_sd)*titer_weights[i];
  }
  return(total_negll);
  
}


//' Calculate the total negative log-likelihood of a mean titer
//' 
//' This is a base function to sum the total _negative_ log likelihood of a mean titer.
//' 
//' @param max_titers Numeric vector of the upper bounds of the measured titers.
//' @param min_titers Numeric vector of the lower bounds of the measured titers.
//' @param predicted_mean The predicted mean titer.
//' @param titer_sd The expected standard deviation of titers.
//' 
//' @details This function calculates the total negative log-likelihood of a predicted mean 
//' titer given a set of titers. The main assumption is that both measurement error and 
//' titer variation are normally distributed. Note that the argument \code{titer_sd} is the 
//' total expected standard deviation of the titer set, i.e. measurement error plus titer 
//' variation.
//' @export
// [[Rcpp::export]]
double calc_mean_titer_negll(NumericVector max_titers,
                             NumericVector min_titers,
                             double predicted_mean,
                             double titer_sd) {
  
  double total_negll = 0;
  for(int i = 0; i < min_titers.length(); ++i) {
    total_negll -= calc_titer_loglik(max_titers[i],
                                     min_titers[i],
                                     predicted_mean,
                                     titer_sd);
  }
  return(total_negll);
  
}

