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
  
  # For p-values less than .05 but greater than .001
  if (p < .05) {
    # Format p-value without scientific notation
    formatted <- sprintf("*p* = %s", format(round(p, digits), nsmall = digits))
    # Remove the leading 0 to match APA style
    formatted <- gsub("0\\.", ".", formatted)
    return(formatted)
  }
  
  # For larger p-values, maintain 3 decimal places
  formatted <- sprintf("*p* = %s", format(round(p, digits), nsmall = digits))
  return(formatted)
}