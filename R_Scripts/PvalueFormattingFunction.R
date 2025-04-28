format_pvalue <- function(p, digits = 3) {
  # Vectorize the function to handle multiple p-values
  if (length(p) > 1) {
    return(sapply(p, format_pvalue, digits = digits))
  }
  
  if (is.na(p)) return("NA")
  
  # For very small p-values
  if (p < .001) {
    return("*p* < .001")
  }
  
  # For all other p-values
  formatted <- sprintf("*p* = %s", format(round(p, digits), nsmall = digits))
  # Remove the leading 0 for all p-values
  formatted <- gsub("0\\.", ".", formatted)
  return(formatted)
}