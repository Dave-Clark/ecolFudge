getGroupDists <- function(adonisResult, data, distanceMatrix){
  factorVector <- attr(adonisResult$terms, "term.label")
  factorCombs <- as.data.frame(expand.grid(data[, factorVector],
    data[, factorVector]))
  factorCombs$ind <- ifelse(factorCombs$Var1 == factorCombs$Var2,
    "Within", "Between")
  factMatr <- matrix(factorCombs$ind, ncol = nrow(data), nrow = nrow(data))
  distanceMatrix <- as.dist(distanceMatrix, diag = F)
  grpDists <- data.frame(distance = as.vector(distanceMatrix),
    comp = as.vector(factMatr[lower.tri(factMatr, diag = F)]))
  return(grpDists)
}
