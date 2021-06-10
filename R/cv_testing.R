
# # Generate CV training and test sets
# lndscp.cv.set <- function(
#   titer_table,
#   train.prop,
#   num.runs
# ){
#   
#   
#   
# }

# Generate CV predictions from a test set
lndscp.cv.predict <- function(
){
  
  # Create test bandwidth and test degree set
  test.pars <- expand.grid(
    test.bandwidths,
    test.degrees
  )[,c(2,1)]
  colnames(test.pars) <- c("degree", "bandwidth")
  
  # Setup for cross-validation runs
  cvruns <- lapply(seq_len(num.runs), function(runnum){
    
    # Generate your training set
    train_table <- t(apply(titer_table, 1, function(titers){
      
      ntrain <- floor(train.prop*sum(!is.na(titers)))
      train_ags <- rep(FALSE, length(titers))
      train_ags[sample(
        x       = which(!is.na(titers)),
        size    = ntrain,
        replace = FALSE
      )] <- TRUE
      train_ags
      
    }))
    
    # Work out your testing set
    test_table <- !is.na(titer_table) & !train_table
    
    # Work out the chull table
    chull_table <- t(vapply(seq_len(nrow(test_table)), function(i){
      
      # Get test and train
      test_ags  <- test_table[i,]
      train_ags <- train_table[i,]
      
      # Check if test within train chull
      within_chull <- rep(NA, length(test_ags))
      within_chull[test_ags] <- apply(
        ag_coords[test_ags,,drop=FALSE], 1, function(point_coords){
          lndscp_checkChull(
            point_coords = point_coords,
            ag_coords    = ag_coords[train_ags,,drop=FALSE]
          )
        }
      )
      
      # Return chull results
      within_chull
      
    }, logical(ncol(test_table))))
    
    # Return the run setups
    list(
      train_table = train_table,
      test_table  = test_table
    )
    
  })
  
}

#' Perform cross-validation testing
#' @export
lndscp.cv.test <- function(
  titer_table,
  ag_coords,
  test.bandwidths,
  test.degrees,
  train.prop,
  num.runs,
  error.sd,
  control = list(),
  verbose = TRUE
){
  
  # Replace * with NA
  titer_table[titer_table == "*"] <- NA
  
  # Setup an empty table
  na_table <- matrix(
    nrow = nrow(titer_table),
    ncol = ncol(titer_table)
  )
  
  # Create test bandwidth and test degree set
  test.pars <- expand.grid(
    test.bandwidths,
    test.degrees
  )[,c(2,1)]
  colnames(test.pars) <- c("degree", "bandwidth")
  
  # Setup for cross-validation runs
  cvruns <- lapply(seq_len(num.runs), function(runnum){
    
    # Generate your training set
    train_table <- t(apply(titer_table, 1, function(titers){
      
      ntrain <- floor(train.prop*sum(!is.na(titers)))
      train_ags <- rep(FALSE, length(titers))
      train_ags[sample(
        x       = which(!is.na(titers)),
        size    = ntrain,
        replace = FALSE
      )] <- TRUE
      train_ags
      
    }))
    
    # Work out your testing set
    test_table <- !is.na(titer_table) & !train_table
    
    # Work out the chull table
    chull_table <- t(vapply(seq_len(nrow(test_table)), function(i){
      
      # Get test and train
      test_ags  <- test_table[i,]
      train_ags <- train_table[i,]
      
      # Check if test within train chull
      within_chull <- rep(NA, length(test_ags))
      within_chull[test_ags] <- apply(
        ag_coords[test_ags,,drop=FALSE], 1, function(point_coords){
          lndscp_checkChull(
            point_coords = point_coords,
            ag_coords    = ag_coords[train_ags,,drop=FALSE]
          )
        }
      )
      
      # Return chull results
      within_chull
      
    }, logical(ncol(test_table))))
    
    # Return the run setups
    list(
      train_table = train_table,
      test_table  = test_table,
      chull_table = chull_table
    )
    
  })
  
  # Setup par prediction tables
  par_predictions <- apply(test.pars, 1, function(pars){
    
    # Track the run number
    if(verbose) message(paste(pars, collapse = ":"))
    
    # Predict for each cv run
    lapply(seq_len(num.runs), function(runnum){
      
      # Track the run number
      if(verbose) message(".", appendLF = FALSE)
      
      # Setup the table of predictions
      prediction_table <- na_table
      
      # Predict for each table row
      for(rownum in seq_len(nrow(titer_table))){
        
        # Get titers and test and training ags
        titers    <- titer_table[rownum,]
        train_ags <- cvruns[[runnum]]$train_table[rownum,]
        test_ags  <- cvruns[[runnum]]$test_table[rownum,]
        
        # Fit training set
        lndscp <- ablandscape.fit(
          titers    = titers[train_ags],
          coords    = ag_coords[train_ags,,drop=FALSE],
          bandwidth = pars[2],
          degree    = pars[1],
          error.sd  = error.sd,
          control   = control
        )
        
        # Predict the test set
        prediction_table[rownum, test_ags] <- predict(
          lndscp,
          coords     = ag_coords[test_ags,,drop=FALSE],
          crop2chull = FALSE
        )
        
      }
      
      # Track last run
      if(runnum == num.runs && verbose) message("done.")
      
      # Return the prediction table
      prediction_table
      
    })
    
  })
  
  # Signal that the fitting is done
  if(verbose) message("done.")
  
  # Return the cv object
  cv_results <- list(
    titer_table          = titer_table,
    ag_coords            = ag_coords,
    test.pars            = test.pars,
    train.prop           = train.prop,
    cv.runs              = cvruns,
    test.par.predictions = par_predictions
  )
  class(cv_results) <- c("lndscp.cv", "list")
  cv_results
  
}


#' Print a summary of cv testing results
#' @export
summary.lndscp.cv <- function(cv_results,
                              crop2chull = TRUE,
                              trim       = 0){
  
  # Get log titers
  log_titers <- titer_to_logtiter(cv_results$titer_table)
  
  # Define lessthans, morethans etc
  lessthans  <- grepl("<", cv_results$titer_table)
  morethans  <- grepl(">", cv_results$titer_table)
  natiters   <- is.na(cv_results$titer_table)
  measurable <- !lessthans & !morethans & !natiters
  
  # For each par set get the breakdown of errors
  par_rmses <- lapply(cv_results$test.par.predictions, function(par.predictions){
    
    measurable_rmses <- rep(NA, length(cv_results$cv.runs))
    lessthan_rmses   <- rep(NA, length(cv_results$cv.runs))
    
    for(runnum in seq_along(cv_results$cv.runs)){
      
      prediction_table <- par.predictions[[runnum]]
      chull_table      <- cv_results$cv.runs[[runnum]]$chull_table
      
      prediction_table[prediction_table < -1] <- -1
      if(crop2chull){ prediction_table[!chull_table] <- NA }
      
      measurable_rmses[runnum] <- sqrt(mean((prediction_table[measurable] - log_titers[measurable])^2, na.rm = TRUE, trim = trim))
      lessthan_rmses[runnum]   <- sqrt(mean((prediction_table[lessthans] - log_titers[lessthans])^2, na.rm = TRUE, trim = trim))
      
    }
    
    list(
      measurable = measurable_rmses,
      lessthan   = lessthan_rmses
    )
    
  })
  
  # Assemble results table
  measurable_rmse <- vapply(par_rmses, function(x){ x$measurable }, numeric(length(cv_results$cv.runs)))
  lessthan_rmse   <- vapply(par_rmses, function(x){ x$lessthan   }, numeric(length(cv_results$cv.runs)))
  
  # Compile results
  output <- data.frame(
    degree    = cv_results$test.pars$degree,
    bandwidth = cv_results$test.pars$bandwidth,
    measurable_rmse_mean  = colMeans(measurable_rmse),
    lessthan_rmse_mean    = colMeans(lessthan_rmse),
    measurable_rmse_var   = apply(measurable_rmse, 2, var),
    lessthan_rmse_var     = apply(lessthan_rmse, 2, var),
    measurable_rmse_nruns = colSums(!is.na(measurable_rmse)),
    lessthan_rmse_nruns   = colSums(!is.na(lessthan_rmse))
  )
  
  # Add a cv summary class
  class(output) <- c("cv.summary", class(output))
  output
  
}


#' Plotting CV summary results
#' @export
plot.cv.summary <- function(cv.summary){
  
  for(titer_type in c("Measurable")){
    
    if(titer_type == "Measurable"){
      rmse_results <- cv.summary$measurable_rmse_mean
    }
    
    if(titer_type == "Non-detectable"){
      rmse_results <- cv.summary$lessthan_rmse_mean
    }
    
    # Plot the summary results
    par(mar = c(5,5,4,4))
    plot.new()
    plot.window(
      xlim = range(cv.summary$bandwidth),
      ylim = min(rmse_results) + c(0, 0.1)
    )
    box()
    
    # Title the plot
    title(
      main = paste(titer_type, "titers"),
      xlab = "Test bandwidth",
      ylab = "Prediction RMSE"
    )
    
    # Draw axis
    axis(side = 1,
         at   = unique(cv.summary$bandwidth))
    
    axis(side = 2,
         las = 1)
    
    # Plot results by degree
    degrees     <- unique(cv.summary$degree)
    degree_cols <- colorRampPalette(c("blue", "#ff3399"))(length(degrees))
    
    # Draw legend
    legend(
      "topright",
      legend = paste("Degree", degrees),
      fill   = degree_cols,
      bty    = "n"
    )
    
    for(i in seq_along(degrees)){
      
      degree     <- degrees[i]
      bandwidths <- cv.summary$bandwidth[cv.summary$degree == degree]
      rmses      <- rmse_results[cv.summary$degree == degree]
      
      # Plot the RMSEs for each test bandwidth
      points(
        x = bandwidths,
        y = rmses,
        col = degree_cols[i],
        pch = 16,
        type = "b"
      )
      
      # Draw a horizontal line at the minimum RMSE
      min_rmse <- min(rmses, na.rm = T)
      abline(h   = min_rmse,
             col = degree_cols[i],
             lty = 3)
      
      # Mark the minimum on the right axis side
      axis(side   = 4,
           at     = min_rmse,
           labels = round(min_rmse, 3),
           col    = degree_cols[i],
           las    = 1)
      
      # Circle the lowest rmse value
      points(
        x = bandwidths[which.min(rmses)],
        y = rmses[which.min(rmses)],
        col = degree_cols[i],
        cex = 1.8
      )
      
    }
    
  }
  
}

