
titer_types <- function(titers){
  
  titer_types  <- titers
  titer_types[]                              <- "measured"
  titer_types[substr(titers, 1, 1) == "<"]   <- "lessthan"
  titer_types[substr(titers, 1, 1) == ">"]   <- "morethan"
  titer_types[titers == "*" | is.na(titers)] <- "omitted"
  titer_types
  
}

titer_to_logtiter <- function(titers){
  
  # Get titer types
  titer_types <- titer_types(titers)
  
  # Get log titers
  threshold_titers <- titer_types == "lessthan" | titer_types == "morethan"
  log_titers                    <- titers
  log_titers[titer_types == "omitted"] <- NA
  log_titers[threshold_titers]  <- substr(log_titers[threshold_titers], 2, nchar(log_titers[threshold_titers]))
  mode(log_titers)              <- "numeric"
  log_titers                    <- log2(log_titers/10)
  log_titers[titer_types == "lessthan"] <- log_titers[titer_types == "lessthan"] - 1
  log_titers[titer_types == "morethan"] <- log_titers[titer_types == "morethan"] + 1
  log_titers
  
}

logtiter_to_titer <- function(logtiters, titer_types, round_titers = TRUE){
  
  # Convert back to raw titers
  titers <- logtiters
  titers[titer_types == "lessthan"] <- titers[titer_types == "lessthan"] + 1
  titers[titer_types == "morethan"] <- titers[titer_types == "morethan"] + 1
  titers <- 2^titers*10
  
  if(round_titers){
    titers <- round(titers)
  }
  
  titers[titer_types == "lessthan"] <- paste0("<", titers[titer_types == "lessthan"])
  titers[titer_types == "morethan"] <- paste0(">", titers[titer_types == "morethan"])
  titers[titer_types == "omitted"]  <- "*"
  titers
  
}

