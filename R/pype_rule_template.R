




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
  rlang::is_integerish(n_rules)
  stopifnot("n_rules must be scalar."=length(n_rules)==1)
  stopifnot("n_rules must be a positive integer."=n_rules>0)
  is_integerish <- is.integer(n_rules) || (is.numeric(n_rules) && 
                                             all(n_rules == trunc(n_rules)) && 
                      !is.na(n_rules))
  stopifnot("n_rules must be a positive integer."=is_integerish)
  rule_call <-  "  rule(  'class_name',     'marker_gene_name',   '>',  1.00 )             %>%"
  for(i in 1:n_rules) cat(rule_call, fill=1)
  cat("  plot_last()")
}