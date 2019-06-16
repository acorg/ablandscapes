#include <Rcpp.h>
using namespace Rcpp;

// This is my modified tobit-like model
//' @export
// [[Rcpp::export]]
double get_lndscp_fit_negll(double max_titre,
                            double min_titre,
                            double pred_titre,
                            double error_sd)
{
  return(R::logspace_sub(R::pnorm5(max_titre,pred_titre,error_sd,1,1),
                         R::pnorm5(min_titre,pred_titre,error_sd,1,1)));
}

//' @export
// [[Rcpp::export]]
NumericVector get_lndscp_prediction_set_negll(NumericVector min_titres,
                                              NumericVector max_titres,
                                              NumericVector pred_titres,
                                              double error_sd)
{
  
  int num_titres = pred_titres.length();
  NumericVector all_neglls(num_titres);
  for(int i = 0; i < num_titres; ++i) {
    all_neglls[i] = -get_lndscp_fit_negll(max_titres[i],
                                          min_titres[i],
                                                    pred_titres[i],
                                                               error_sd);
  }
  return(all_neglls);
  
}

//' @export
// [[Rcpp::export]]
double get_negll_hi_lm(NumericVector par,
                       NumericMatrix max_titres,
                       NumericMatrix min_titres,
                       NumericMatrix ag_coords,
                       NumericVector ag_weights,
                       double error_sd)
{
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(int i = 0; i < ag_coords.nrow(); ++i) {
    if(ag_weights[i] > 0) {
      predicted_titre = par[0];
      for(int j = 0; j < ag_coords.ncol(); ++j) {
        predicted_titre += par[j+1] * ag_coords(i,j);
      }
      for(int x = 0; x < max_titres.nrow(); ++x) {
        if(!NumericVector::is_na(max_titres(x,i))) {
          log_prediction_prob = get_lndscp_fit_negll(max_titres(x,i), min_titres(x,i), predicted_titre, error_sd);
          overall_prob -= log_prediction_prob * ag_weights[i];
        }
      }
    }
  }
  
  return(overall_prob);
  
}


//' @export
// [[Rcpp::export]]
double get_negll_hi_height(NumericVector par,
                           double lndscp_height,
                           NumericMatrix max_titres,
                           NumericMatrix min_titres,
                           NumericMatrix ag_coords,
                           NumericVector ag_weights,
                           double error_sd)
{
  
  par.insert(0, lndscp_height);
  return(get_negll_hi_lm(par,
                         max_titres,
                         min_titres,
                         ag_coords,
                         ag_weights,
                         error_sd));
  
}




//' @export
// [[Rcpp::export]]
double get_negll_lm(NumericVector par,
                    NumericVector pred_titres,
                    NumericVector max_titres,
                    NumericVector min_titres,
                    double error_sd)
{
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(int i = 0; i < pred_titres.length(); ++i) {
    predicted_titre = par[0] + par[1]*pred_titres[i];
    if(!NumericVector::is_na(max_titres[i])) {
      log_prediction_prob = get_lndscp_fit_negll(max_titres[i], min_titres[i], predicted_titre, error_sd);
      overall_prob -= log_prediction_prob;
    }
  }
  
  return(overall_prob);
  
}



//' @export
// [[Rcpp::export]]
double get_negll_lm_intercept(NumericVector par,
                              NumericVector pred_titres,
                              NumericVector max_titres,
                              NumericVector min_titres,
                              double lm_gradient,
                              double error_sd)
{
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(int i = 0; i < pred_titres.length(); ++i) {
    predicted_titre = par[0] + lm_gradient*pred_titres[i];
    if(!NumericVector::is_na(max_titres[i])) {
      log_prediction_prob = get_lndscp_fit_negll(max_titres[i], min_titres[i], predicted_titre, error_sd);
      overall_prob -= log_prediction_prob;
    }
  }
  
  return(overall_prob);
  
}



//' @export
// [[Rcpp::export]]
double get_negll_lm_gradient(NumericVector par,
                             NumericVector pred_titres,
                             NumericVector max_titres,
                             NumericVector min_titres,
                             double lm_intercept,
                             double error_sd)
{
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(int i = 0; i < pred_titres.length(); ++i) {
    predicted_titre = lm_intercept + par[0]*pred_titres[i];
    if(!NumericVector::is_na(max_titres[i])) {
      log_prediction_prob = get_lndscp_fit_negll(max_titres[i], min_titres[i], predicted_titre, error_sd);
      overall_prob -= log_prediction_prob;
    }
  }
  
  return(overall_prob);
  
}




