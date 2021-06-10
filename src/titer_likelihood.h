
#include <RcppArmadillo.h>

#ifndef ablandscapes__titer_likelihood__h
#define ablandscapes__titer_likelihood__h

double titer_prediction_negll(
    const double &max_titer,
    const double &min_titer,
    const double &pred_titer,
    const double &error_sd
);

#endif
