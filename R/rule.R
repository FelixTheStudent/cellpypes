


#' Add 
#'
#' @param obj 
#' @param class 
#' @param feature 
#' @param operator 
#' @param t 
#' @param parent 
#'
#' @return
#' @export
#'
#' @examples
rule <- function(obj,
                 class,
                 feature,
                 operator = ">",
                 threshold,
                 parent = NULL) {
  if( any(is.null(obj), is.null(class), is.null(feature), is.null(threshold) )) { 
    stop("rule does not accept NULL input for arguments: obj, class, feature and threshold.", call. = F) }
  stopifnot(all(
    is_scalar_character(class), is_scalar_character(feature),
    is_scalar_character(operator), is_scalar_character(parent)|is.null(parent),
    is_scalar_double(threshold) ) )
  if(any(parent==class)) stop("Class and parent cannot be the same.")

  
  
  
  return(obj)
}


# examples:
#   obj %>% rule("T", "CD3E", ">", .1)
#   obj %>% rule("T", "percent.mito", "<", .25)

# write `rule` function by completing each step in the following:
(rule_cases <- expand.grid(list(existing_class=c(F,T), existing_feature=c(T,F), NULL_parent = c(T,F))))
cbind(rule_cases, action = c(
 "new class can't have existing feature, throw error", 
 1,
 1,
 1, 
 "new class can't have existing feature, throw error", 
 1, 
 "create class, create rule, parent=parent",
 "modify class, create rule, modify parent"
 )) 

# compute pooled counts if(existing_feature=F)


