

#' Check if obj$classes looks as expected.
#' is_class returns FALSE for example in these cases:
#'   is_classes(NULL)
#'   is_classes(data.frame())
#'   is_classes(data.frame(class=c("T","T"), parent=c("..root..","..root..")))
#' @param classes 
#'
#' @return
#'
#' @examples
is_classes <- function(classes) {
    inherits(classes, "data.frame"  ) &&
    ncol(classes) > 0 &&  # handles this: `is_classes(data.frame())`
    all(colnames(classes)[1:2] == c("class", "parent")) &&
    all(!duplicated(classes$class))
  }

#' Check if obj$rules looks as expected.
#'
#' @param rules 
#'
#' @return
#'
#' @examples
is_rules <- function(rules) {
    inherits(rules, "data.frame"  ) &&
    ncol(rules) > 0 &&  # handles this: `is_rules(data.frame())`
    all(colnames(rules)[1:4] == c("class","feature","operator","threshold")) &&
    all(!duplicated(paste0(rules$class, rules$feature)))
  }



#' Add a cell type rule.
#'
#' @param obj 
#' @param class 
#' @param feature 
#' @param operator One of the following: \code{c("<", ">", "<=", ">=")}.
#' @param t 
#' @param parent The parent class (e.g. "T" or "T cell").
#' If NULL, new classes get "..root.." and
#' existing classes keep their current parent.
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
  # check inputs:
  if( any(is.null(obj), is.null(class), is.null(feature), is.null(threshold) )) { 
    stop("rule does not accept NULL input for arguments: obj, class, feature and threshold.", call. = F) }
  stopifnot(all(
    rlang::is_scalar_character(class),    rlang::is_scalar_character(feature),
    rlang::is_scalar_character(operator), rlang::is_scalar_character(parent)|is.null(parent),
    rlang::is_scalar_double(threshold) ) )
  stopifnot(operator %in% c("<", ">", "<=", ">="))
  if(any(parent==class)) stop("Class and parent cannot be the same.")
  if(any(parent=="..root..")) stop("Cellpypes internally uses '..root..'. Call your parent something else!")

  # check that obj has everything this rule needs:
  stopifnot(feature %in% colnames(obj$raw))
  stopifnot(is.null(obj$classes) || is_classes(obj$classes) )
  stopifnot(is.null(obj$rules) || is_rules(obj$rules) )
  stopifnot(is.null(obj$classes) == is.null(obj$rules))
  
  # is_classes also checks if obj$class is NULL
  existing_class <- ifelse(is_classes(obj$class) && class %in% obj$class,
                           TRUE,
                           FALSE)
  existing_feature <- ifelse(is_rules(obj$rules) && class %in% obj$rules$class &&
                               feature %in% obj$rules$feature,
                             TRUE,
                             FALSE)
  
  # design choices:
  #  * parent=NULL means "use existing parent"
  #  * top classes have parent="..root..":
  if (is.null(parent)) {
    parent <- ifelse(existing_class,
                     obj$classes$parent[obj$classes==class],
                     "..root..")
  }
  
  class_dat <- data.frame(class=class, parent=parent)
  rule_dat  <- data.frame(class=class, feature=feature, 
                          operator=operator, threshold=threshold)
  # save to object
  if (existing_class) {
    obj$classes[obj$classes$class==class,] <- class_dat
  } else {
    obj$classes <- rbind(class_dat, obj$classes)
  }
  if (existing_feature) {
    obj$rules[obj$rules$class==class && obj$rules$feature==feature,] <- rule_dat
  } else {
    obj$rules <- rbind(rule_dat, rules) # newest rule is on top
  }
  
  
  
  return(obj)
}


# examples:
#   obj %>% rule("T", "CD3E", ">", .1)
#   obj %>% rule("T", "percent.mito", "<", .25)

# write `rule` function by completing each step in the following:
rule_cases <- expand.grid(list(existing_class=c(F,T), existing_feature=c(T,F), NULL_parent = c(T,F)))
cbind(rule_cases, action = c(
 "new class can't have existing feature, throw error", 
 "modify class, modify rule, parent=parent",
 "create class, create rule, parent=parent",
 1, 
 "new class can't have existing feature, throw error", 
 1, 
 "create class, create rule, create parent",
 "modify class, create rule, modify parent"
 )) 

# compute pooled counts if(existing_feature=F)


