

# function to calculate CV
calc_CV <- function(index) {
  sd(index, na.rm = TRUE) / mean(index, na.rm = TRUE)
}


# function to calculate CD at log scale
calc_CD <- function(index) {
  index <- index[!is.na(index)]
  
  # Ensure positive values (required for log)
  if (any(index <= 0)) {
    stop("Index contains zero or negative values; CD (log-ratio) is undefined.")
  }
  
  mean(abs(diff(log(index))))
}




# function to detect trend difference between two time series based on Mann-Kendall test
mk_trend_difference <- function(index_A, index_B) {
  
  # Basic checks
  if (length(index_A) != length(index_B)) {
    stop("index_A and index_B must have the same length.")
  }
  
  if (any(index_A <= 0, na.rm = TRUE) || any(index_B <= 0, na.rm = TRUE)) {
    stop("Indices must be positive to compute log-ratios.")
  }
  
  # Remove paired NA values
  keep <- complete.cases(index_A, index_B)
  index_A <- index_A[keep]
  index_B <- index_B[keep]
  
  # Log-ratio series
  log_ratio <- log(index_A / index_B)
  
  # Mann–Kendall test
  library(Kendall)
  mk <- MannKendall(log_ratio)
  
  # Return key results
  list(
    tau = mk$tau,
    p_value = mk$sl,
    interpretation = ifelse(
      mk$sl < 0.05,
      ifelse(mk$tau > 0,
             "Series A increases faster than Series B",
             "Series B increases faster than Series A"),
      "No significant difference in monotonic trends"
    )
  )
}