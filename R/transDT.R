transDT <- function(dt, transCol, rowID){
  newRowNames <- colnames(dt)
  newColNames <- dt[, transCol, with = F]
  transposedDt <- transpose(dt[, !colnames(dt) %in% transCol, with = F])
  colnames(transposedDt) <- unlist(newColNames)
  transposedDt[, rowID] <- newRowNames[newRowNames != transCol]
  return(transposedDt)
}
