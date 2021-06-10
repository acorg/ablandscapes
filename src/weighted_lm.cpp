
#include <RcppArmadillo.h>

// [[Rcpp::export]]
double get_lndscp_fit_negll(
    const double &max_titre,
    const double &min_titre,
    const double &pred_titre,
    const double &error_sd
  ) {
  return(
    R::logspace_sub(
      R::pnorm5(max_titre,pred_titre,error_sd,1,1),
      R::pnorm5(min_titre,pred_titre,error_sd,1,1)
    )
  );
}

// [[Rcpp::export]]
arma::vec get_lndscp_prediction_set_negll(
    const arma::vec &min_titres,
    const arma::vec &max_titres,
    const arma::vec &pred_titres,
    const double &error_sd
  ) {
  
  int num_titres = pred_titres.n_elem;
  arma::vec all_neglls(num_titres);
  for(int i = 0; i < num_titres; ++i) {
    all_neglls(i) = -get_lndscp_fit_negll(
      max_titres(i),
      min_titres(i),
      pred_titres(i),
      error_sd
    );
  }
  return(all_neglls);
  
}

// [[Rcpp::export]]
double get_negll_hi_lm(
    const arma::vec &par,
    const arma::mat &max_titres,
    const arma::mat &min_titres,
    const arma::mat &ag_coords,
    const arma::vec &ag_weights,
    const double &error_sd
  ) {
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(arma::uword i = 0; i < ag_coords.n_rows; ++i) {
    if(ag_weights(i) > 0) {
      predicted_titre = par(0);
      for(arma::uword j = 0; j < ag_coords.n_cols; ++j) {
        predicted_titre += par(j+1) * ag_coords(i,j);
      }
      for(arma::uword x = 0; x < max_titres.n_rows; ++x) {
        if(std::isfinite(max_titres(x,i))) {
          log_prediction_prob = get_lndscp_fit_negll(
            max_titres(x,i), 
            min_titres(x,i), 
            predicted_titre, 
            error_sd
          );
          overall_prob -= log_prediction_prob * ag_weights(i);
        }
      }
    }
  }
  
  return(overall_prob);
  
}


// [[Rcpp::export]]
double get_negll_hi_height(
    arma::vec par,
    const double &lndscp_height,
    const arma::mat &max_titres,
    const arma::mat &min_titres,
    const arma::mat &ag_coords,
    const arma::vec &ag_weights,
    const double &error_sd
  ) {
  
  par.insert_rows(0, lndscp_height);
  return(
    get_negll_hi_lm(
      par,
      max_titres,
      min_titres,
      ag_coords,
      ag_weights,
      error_sd
    )
  );
  
}

// [[Rcpp::export]]
double get_negll_lm(
    const arma::vec &par,
    const arma::vec &pred_titres,
    const arma::vec &max_titres,
    const arma::vec &min_titres,
    const double &error_sd
  ) {
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(arma::uword i = 0; i < pred_titres.n_elem; ++i) {
    predicted_titre = par(0) + par(1)*pred_titres(i);
    if(std::isfinite(max_titres(i))) {
      log_prediction_prob = get_lndscp_fit_negll(
        max_titres(i), 
        min_titres(i), 
        predicted_titre, 
        error_sd
      );
      overall_prob -= log_prediction_prob;
    }
  }
  
  return(overall_prob);
  
}

// [[Rcpp::export]]
double get_negll_lm_intercept(
    const arma::vec &par,
    const arma::vec &pred_titres,
    const arma::vec &max_titres,
    const arma::vec &min_titres,
    const double &lm_gradient,
    const double &error_sd
  ) {
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(arma::uword i = 0; i < pred_titres.n_elem; ++i) {
    predicted_titre = par(0) + lm_gradient*pred_titres(i);
    if(std::isfinite(max_titres(i))) {
      log_prediction_prob = get_lndscp_fit_negll(
        max_titres(i), 
        min_titres(i), 
        predicted_titre, 
        error_sd
      );
      overall_prob -= log_prediction_prob;
    }
  }
  
  return(overall_prob);
  
}

// [[Rcpp::export]]
double get_negll_lm_gradient(
    const arma::vec &par,
    const arma::vec &pred_titres,
    const arma::vec &max_titres,
    const arma::vec &min_titres,
    const double &lm_intercept,
    const double &error_sd
  ) {
  
  double overall_prob = 0;
  double predicted_titre;
  double log_prediction_prob;
  
  for(arma::uword i = 0; i < pred_titres.n_elem; ++i) {
    predicted_titre = lm_intercept + par(0)*pred_titres(i);
    if(std::isfinite(max_titres(i))) {
      log_prediction_prob = get_lndscp_fit_negll(
        max_titres(i), 
        min_titres(i), 
        predicted_titre, 
        error_sd
      );
      overall_prob -= log_prediction_prob;
    }
  }
  
  return(overall_prob);
  
}
