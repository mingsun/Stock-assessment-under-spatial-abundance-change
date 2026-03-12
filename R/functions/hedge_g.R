# hedge's g function

calc_hedge_s <- function(group.1, group.2) {
  
  # treatment group is 1
  # group.1 <- filter(catch_by_tow, OWF == "INSIDE")
  n.1 = nrow(group.1)
  mean.1 = mean(group.1$NUMBER_adj, na.rm = TRUE)
  sd.1 =  sd(group.1$NUMBER_adj, na.rm = TRUE)
  
  # control group is 2
  # group.2 <- filter(catch_by_tow, OWF == "OUTSIDE")
  n.2 = nrow(group.2)
  mean.2 = mean(group.2$NUMBER_adj, na.rm = TRUE)
  sd.2 =  sd(group.2$NUMBER_adj, na.rm = TRUE)
  
  # pooled
  df = n.1 + n.2 -2
  sd.pooled <- sqrt(((n.1 - 1) * sd.1^2 + (n.2 - 1) * sd.2^2)/df)
  J = 1-3/4/df
  var.pooled = ((n.1+n.2)/(n.1 * n.2) + (sd.1-sd.2)^2/(2 * (n.1+n.2))) * J^2
  
  g = (mean.1 - mean.2)/sd.pooled * J
  
  return(g)
}
  
