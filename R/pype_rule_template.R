




#' Generate code template for cellpype rules
#'
#' @param n_rules Number of lines (rules) to generate
#' 
#' @description This function \link[cellpypes]{rule} code snippet with neat 
#' text alignment to the console.
#' Paste this into your script and start changing the rules.
#'
#' @return Prints rules to the consoles.
#' @export
#'
#' @examples
#' pype_code_template()
pype_code_template <- function(n_rules=3) {
  assertthat::is.count(n_rules)
  rule_call <-  "  rule(  'class_name',     'marker_gene_name',   '>',  1.00 )             %>%"
  for(i in 1:n_rules) cat(rule_call, fill=1)
  cat("  plot_last()")
}