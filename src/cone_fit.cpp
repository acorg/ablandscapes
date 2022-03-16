
#include <RcppArmadillo.h>
#include "titer_likelihood.h"

double euc_dist(
    const arma::rowvec &coord1,
    const arma::rowvec &coord2
) {
  
  return(
    std::sqrt(
      std::pow(coord1(0) - coord2(0), 2) +
      std::pow(coord1(1) - coord2(1), 2)
    )
  );
  
}

// [[Rcpp::export]]
double negll_cone_pars(
    const arma::vec &par,
    const arma::vec &cone_heights,
    double cone_slope,
    arma::mat cone_coords,
    const arma::mat &ag_coords,
    const arma::mat &max_titers,
    const arma::mat &min_titers,
    const double &error_sd,
    const bool optimise_cone_slope,
    const bool optimise_cone_coords
) {
  
  // Set parameters
  if (optimise_cone_slope) {
    
    cone_slope = par(0);
    if (optimise_cone_coords) {
      for (arma::uword sr = 0; sr < max_titers.n_rows; sr++) {
        cone_coords(sr, 0) = par(sr*2 + 1);
        cone_coords(sr, 1) = par(sr*2 + 2);
      }
    }
    
  } else {
    
    if (optimise_cone_coords) {
      for (arma::uword sr = 0; sr < max_titers.n_rows; sr++) {
        cone_coords(sr, 0) = par(sr*2);
        cone_coords(sr, 1) = par(sr*2 + 1);
      }
    }
    
  }
  
  // Calculate likelihood
  double negll = 0;
  for (arma::uword sr = 0; sr < max_titers.n_rows; sr++) {
    for (arma::uword ag = 0; ag < ag_coords.n_rows; ag++) {
      
      double dist = euc_dist(cone_coords.row(sr), ag_coords.row(ag));
      double predicted_titer = cone_heights(sr) - dist*cone_slope;
      if(std::isfinite(max_titers(sr, ag))) {
        negll -= titer_prediction_negll(
          max_titers(sr,ag), 
          min_titers(sr,ag), 
          predicted_titer, 
          error_sd
        );
      }
      
    }
  }
  return(negll);
  
}
