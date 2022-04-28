

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
  # Helpful messages to the user:
  if(!all(classes$parent %in% c("..root..",classes$class))) {
    stop("A class has parent that does not exist -- double-check your rules!")
  }
    inherits(classes, "data.frame"  ) &&
    ncol(classes) > 0 &&  # handles this: `is_classes(data.frame())`
    all(colnames(classes)[1:2] == c("class", "parent")) &&
    all(!duplicated(classes$class)) &&
    all(classes$parent %in% c("..root..",classes$class))
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



check_obj <- function(obj) {
  # check that object has everything we expect from obj.
  stopifnot(is.null(obj$classes) || is_classes(obj$classes) )
  stopifnot(is.null(obj$rules) || is_rules(obj$rules) )
  stopifnot(is.null(obj$classes) == is.null(obj$rules))
  
  # If rules exist they make sense:
  if (is_rules(obj$rules)) {
    stopifnot( all(obj$rules$feature %in% rownames(obj$raw)) )
  }
  
  # If totalUMI exist they make sense:
  if (!is.null(obj$totalUMI)) {
    s <- obj$totalUMI
    stopifnot( length(s) == ncol(obj$raw) )
    stopifnot( all.equal(s, as.integer(s), check.attributes=FALSE) )
  }
}



#' Add a cell type rule.
#'
#' @param obj 
#' @param class 
#' @param feature 
#' @param operator One of the following: \code{c("<", ">", "<=", ">=")}.
#' @param threshold 
#' @param use_CP10K Is the threshold provided as counts per 10 thousand UMIs
#' (CP10K)? If TRUE, cellpypes multiplies the threshold with 1e-4.
#' @param parent The parent class (e.g. "T" or "T cell").
#' If NULL, new classes get "..root.." and
#' existing classes keep their current parent.
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom rlang is_scalar_character is_scalar_double
rule <- function(obj,
                 class,
                 feature,
                 operator = ">",
                 threshold,
                 parent = NULL,
                 use_CP10K=TRUE) {
  # check inputs:
  check_obj(obj)
  if( any(is.null(obj), is.null(class), is.null(feature), is.null(threshold) )) { 
    stop("rule does not accept NULL input for arguments: obj, class, feature and threshold.", call. = F) }
  stopifnot(all(
    is_scalar_character(class),    is_scalar_character(feature),
    is_scalar_character(operator), is_scalar_character(parent)|is.null(parent),
    is_scalar_double(threshold) ) )
  stopifnot(operator %in% c("<", ">", "<=", ">="))
  if(any(parent==class)) stop("Class and parent cannot be the same.")
  if(any(parent=="..root..")) stop("Cellpypes internally uses '..root..'. Call your parent something else!")

  
  # is_classes also checks if obj$classes is NULL
  existing_class <- ifelse(
    is_classes(obj$classes) && class %in% obj$classes$class,
    TRUE,
    FALSE)
  existing_rule <- ifelse(
    is_rules(obj$rules) && 
    any( obj$rules$class==class & obj$rules$feature == feature ),
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
  if(use_CP10K) {threshold <- threshold * 1e-4}
  rule_dat  <- data.frame(class=class, feature=feature, 
                          operator=operator, threshold=threshold)
  # save to object
  if (existing_class) {
    obj$classes[obj$classes$class==class,] <- class_dat
  } else {
    obj$classes <- rbind(obj$classes, class_dat) # newest class comes last
  }
  # Updating a rule also moves it to the position of most recent rule:
  if (existing_rule) {
    modified_rule <- obj$rules$class==class & obj$rules$feature==feature
    obj$rules <- obj$rules[!modified_rule,] 
  } 
  obj$rules <- rbind(obj$rules, rule_dat) # newest rule comes last
  
  
  
  return(obj)
}
# backlog:
#  * compute pooled counts if feature already exists in rules$feature
#  * add examples
#  * perhaps allow using meta-data such as percent.mito or totalUMI in rules. 
#    This is how it could look like:
#   obj %>% rule("T", "CD3E", ">", .1)
#   obj %>% rule("T", "percent.mito", "<", .25)


# I have written the `rule` function by considering each case in the following:
#       (note from 2 weeks later: existing_feature now is existing_rule, not
#        sure if the table below still makes sense to look at)
# rule_cases <- expand.grid(list(existing_class=c(F,T), existing_feature=c(T,F), NULL_parent = c(T,F)))
# cbind(rule_cases, action = c(
#  "new class can't have existing feature, throw error", 
#  "modify class, modify rule, parent=parent",
#  "create class, create rule, parent=parent",
#  1, 
#  "new class can't have existing feature, throw error", 
#  1, 
#  "create class, create rule, create parent",
#  "modify class, create rule, modify parent"
#  )) 



