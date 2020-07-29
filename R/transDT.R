#' A function for quickly transposing data.table objects
#'
#' This function allows you transpose data.table objects very rapidly, whilst maintaining control of row and column names.
#' @param dt The data.table object you wish to transpose.
#' @param transCol The name of the column that you wish to pivot on. Values in this column will become the new colnames.
#' @param rowID This will be the name of the new rownames column
#' @keywords data.table
#' @export
#' @examples
#' transDT()

transDT <- function(dt, transCol, rowID){
  newRowNames <- colnames(dt)
  newColNames <- dt[, transCol, with = F]
  transposedDt <- transpose(dt[, !colnames(dt) %in% transCol, with = F])
  colnames(transposedDt) <- unlist(newColNames)
  transposedDt[, rowID] <- newRowNames[newRowNames != transCol]
  return(transposedDt)
}
