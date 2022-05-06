

#' Check if obj$classes looks as expected.
#' is_class returns FALSE for example in these cases:
#'   is_classes(NULL)
#'   is_classes(data.frame())
#'   is_classes(data.frame(class=c("T","T"), parent=c("..root..","..root..")))
#' @param classes The obj$classes you want to check.
#'
#' @return logical scalar.
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
#' @param rules The obj$rules slot of a cellpypes object.
#'
#' @return logical scalar
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
#' @template param_obj
#' @param class Character scalar with the class name. Typically, 
#' cellpypes classes
#' are literature cell types ("T cell") or any subpopulation of interest
#' ("CD3E+TNF+LAG3-"). 
#' @param feature Character scalar naming the gene you'd like to threshold. 
#' Must be a row name in \code{obj$raw}.
#' @template param_operator
#' @param threshold Numeric scalar with the expression threshold separating positive
#' from negative cells.
#' Experiment with this value, until expression and selected cells agree well
#' in UMAP (see examples on 
#' \href{https://github.com/FelixTheStudent/cellpypes}{gitHub}).
#' @param use_CP10K If TRUE, \code{threshold} is taken to be 
#' counts per 10 thousand UMI counts, a measure for RNA molecule fractions. 
#' We recommend CP10K for human intuition (1 CP10K is roughly 1 UMI per cell),
#' but the results are the exact same whether you use  
#' \code{threshold=1,CP10K=TRUE} or
#' \code{threshold=1e-4,CP10K=FALSE}.
#' @param parent Character scalar with the parent class 
#' (e.g. "T cell" for "Cytotoxic T cells").
#' Only has to be specified once per class (else most recent one is taken),
#' and defaults to
#' "..root.." if NULL is passed in all rules.
#' 
#' @description This is the heart of cellpypes and best used by piping from
#' one rule into the next
#' with `magrittr::%>%`. Check out examples at
#' \href{https://github.com/FelixTheStudent/cellpypes}{gitHub})!
#' 
#' @details Calling `rule` is computationally cheap because it only stores
#' the cell type rule while all computations
#' happen in \link[cellpypes]{classify}.
#' If you have classes with multiple rules, the most recent \code{parent} and
#' \code{feature}-\code{threshold} combination counts.
#' It is ok to mix rules with and without \code{use_CP10K=TRUE}.
#' 
#' 
#' @seealso 
#' To have nicely formatted code in the end, copy the output of
#' \code{pype_code_template()} to your script and start editing.
#' 
#' 
#' @return \code{obj} is returned, but with the rule and class stored in
#' \code{obj$rules} and \code{obj$classes}, to be used by
#' \link[cellpypes]{classify}.
#' 
#' @template cellpypes_obj
#' 
#' @export
#'
#' @examples
#' # T cells are CD3E+:
#' obj <- rule(simulated_umis, "T", "CD3E", ">", .1)
#' # T cells are MS4A1-:
#' obj <- rule(obj, "T", "MS4A1", "<", 1)
#' # Tregs are a subset of T cells:
#' obj <- rule(obj, "Treg", "FOXP3", ">", .1, parent="T") 
#' 
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
  stopifnot(    "class must be a single string"= is.character(class) && length(class) == 1)
  stopifnot(  "feature must be a single string"= is.character(feature) && length(feature) == 1)
  stopifnot( "operator must be a single string"= is.character(operator) && length(operator) == 1)
  stopifnot("threshold must be a single number"= is.numeric(threshold) && length(threshold) == 1)
  stopifnot(   "parent must be NULL or a single string"= is.null(parent) | is.character(parent) && length(parent) == 1)
  stopifnot("operator has to be one of these: <, >"=operator %in%
              c("<", ">"))
  if(any(parent==class)) stop("Class and parent cannot be the same.")
  if(any(parent=="..root..")) stop("Cellpypes internally uses '..root..'. Call your parent something else!")
  stopifnot("feature must be in rownames(obj$raw)"= feature %in% rownames(obj$raw) )
  
  
  
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



