# Functions that return classification outcomes (class assignments)



# classes encodes a phylogenic tree of classes/cell types.
# These are the functions to work with that:
#  tree_ancestry     function(tree, class)   returns class and its ancestry
#  tree_leafs        function(tree)          returns names of all leaf nodes   
#  tree_plot         function(obj)           ggplot object or something :)


#' Finds leaf nodes, i.e. classes without children
#'
#' @param classes. Put obj$classes here, i.e. a tree of class definitions 
#' created with the rule function.
#'
#' @return
#' @export
#'
#' @examples
tree_leaf_nodes <- function(classes) {
  stopifnot(is_classes(classes))
  
  res <- classes
  res$hasChild <- FALSE
  res$hasChild[res$class %in% res$parent] <- TRUE
 
  res$class[!res$hasChild] 
}




#' Find parent, parent's parent and so on for a class using recursive programming
#'
#' @param classes The class definitions of a cellpypes object, i.e. obj$classes.
#' @param class A character vector of length one with the class.
#'
#' @return
#' @export
#'
tree_ancestry <- function(classes, class) {
  sel <- classes$class==class
  current_parent <- classes[sel, "parent"]
  current_class  <- classes[sel, "class"] 
  
  if (current_parent == "..root..") {
    return(current_class)
  } else {
    return(c(tree_ancestry(classes, current_parent), current_class))
  }
}







#' Classify cells with boolean logic
#'
#' @param obj Cellpypes object.
#' @param classes Character vector with one or more class names.
#' @param drop If TRUE (default) and only a single class is supplied,
#' a boolean vector is returned. A boolean matrix otherwise.
#'
#' @return A boolean vector or matrix.
#' @export
#'
#' @examples
# class_boolean <- function(obj, classes, drop=TRUE) {
#   stopifnot(is.character(classes))
#   stopifnot(is_classes(obj$classes))
#   stopifnot(is_rules(obj$rules))
#   
# 
#   boolean_matrix <- sapply(classes, function(class) { # loop through classes
#     ancestry <- tree_ancestry(obj$classes, class)
#     class_rules <- obj$rules[obj$rules$class %in% ancestry,]
#     # Reduce can be thought of here as looping through class_rules:
#     base::Reduce(function(x,y) x&y, 
#                 lapply(1:nrow(class_rules), function(i){
#                   evaluate_rule(obj      = obj,
#                                 class    = class_rules[i, "class"],
#                                 feature  = class_rules[i, "feature"],
#                                 operator = class_rules[i, "operator"],
#                                 threshold= class_rules[i, "threshold"]
#                                 ) })
#                 )
#   })
#   
#   if (ncol(boolean_matrix)==1 && drop) {
#     drop(boolean_matrix)
#   } else{
#     boolean_matrix
#   }
#      
# }






# obj example to develop functions here:
# obj <- simulated_umis %>%
#   rule("B", "MS4A1", ">", .1e-4)                     %>%
#   rule("T", "CD3E", ">", .1e-4)                      %>%
#   rule("Treg", "FOXP3", ">", .1e-4, parent="T")      %>%
#   rule("Tother", "FOXP3", "<", .1e-4, parent="T")    %>%
#   rule("actTreg", "ICOS", ">", .1e-4, parent="Treg")





#' Title
#'
#' @param obj 
#' @param classes 
#' @param replace_overlap_with 
#' @param replace_unassigned_with 
#' @param return_logical_matrix 
#'
#' @return I think the order of factor levels is the same as classes, and
#' if is.null(classes) then it's the same as in obj$classes$class, i.e. the
#' order in which classes were created by the user. I am not sure, however,
#' but realize this is important (see e.g. f1_all.Rmd, building an f1-data.frame
#' relies on it).
#' @export
#'
#' @examples
classify <- function(
  obj,
  classes=NULL,
    # If NULL, uses all childless classes (leafs).
    # unique(obj$classes$class) returns both leafs and parents.
  replace_overlap_with="Unassigned", # alternatives: 'common_parent', NA or any scalar character
  replace_unassigned_with="Unassigned", # "common_parent", NA or any scalar character
  return_logical_matrix =FALSE # ignore overlap/unassigned rules and just output 
  # a logical matrix. If a single class is supplied, the matrix has exactly one
  # column and the user can pipe it into "drop" to convert it to a vector.
) {
  # checks specific for classify: classes have to be present in obj$classes
  # other sanity checks:
  check_obj(obj) 
  

  if (is.null(classes)) {
    # Convention: use tree leaves if classes=NULL 
    classes <- tree_leaf_nodes(obj$classes)
    # when using leaves, rules from all classes are relevant for classification:
    relevant_classes <- obj$classes$class
  } else {
    stopifnot(all(classes %in% obj$classes$class))
    # To omit irrelevant computations below I note down relevant classes here:
    relevant_classes <- lapply(classes, tree_ancestry, classes=obj$classes)
    relevant_classes <- unique(do.call(c, relevant_classes))
  }

  
  
  
  # rules_eval:
  #   It is important to evaluate once for each rule (rather than once for 
  #   each class), since in the latter case the computation is repeated for each
  #   (grand)child node.
  # rules_info:
  #   Its columns correspond directly to the columns in rules_eval by design,
  #   so we can use rules_info to subset rules_eval below.
  rules_info <-obj$rules[obj$rules$class %in% relevant_classes, ]
  rules_eval <- mapply(
    FUN = function(feature, operator, threshold) {
      K <- pool_across_neighbors(obj$raw[, feature], 
                                 obj$neighbors)
      S <- pool_across_neighbors(Matrix::rowSums(obj$raw),
                                 obj$neighbors)
      cdf <- ppois(K, S*threshold)
      switch(operator,
             ">" = cdf > .99,
             ">=" =cdf > .01,
             "<"  =cdf < .01,
             "<=" =cdf < .99)
      },
    rules_info$feature,
    rules_info$operator,
    rules_info$threshold)
  # mapply sets (non-unique) features as colnames, preclude confusion with NULL:
  colnames(rules_eval) <- NULL

  res <- sapply(relevant_classes, function(class) {
    # We can use info to subset eval because their columns directly correspond:
    is_class_rule <- rules_info$class %in% tree_ancestry(obj$classes, class)
    sum(is_class_rule) == rowSums(rules_eval[, is_class_rule, drop=FALSE])
  })
  
    
  # massage according to replace_overlap_with, replace_unassigned_with, etc.
  if (return_logical_matrix) {
    return(res[, classes, drop=FALSE])
    }
 
  if (replace_overlap_with!="Unassigned" || replace_unassigned_with!="Unassigned")
    stop("Not implemented yet, sorry.")
  # Note to myself: for "common_parent" res[, classes] is not enough.
 
   
  class_res <- res[, classes]
  class_factor <- rep("Unassigned", nrow(class_res))
  for (i in 1:ncol(class_res)) {
    class_factor[ class_res[,i] ] <- classes[i]
  }
  class_factor[base::rowSums(class_res)>1] <- replace_overlap_with
  
  # Unassigned class may or may not be present, under different names:
  #  1. convert to factor at very end (adding new class to factor is hard)
  #  2. It is a design choice whether to add classes with 0 cells or not:
  class_factor <- factor( class_factor,
                          levels=c(classes, "Unassigned"))
  return( class_factor )
}




# skeleton:
# obj %>%
#   rule(  "A",       "gene1",          ">",   1e-3              )             %>%
#   rule(  "A",       "gene1",          ">",   1e-3, parent = "A")             %>%
#   rule(  "A",       "gene1",          ">",   1e-3, parent = "A")             



# Alternative parametrization, which I think is stupid (delete soon):
# classes = c("Treg", "Ttox", "B"). 
# 
# overlap = c("common_parent", NA or any character string: "Unassigned"/"multiplelabels"/...)
# resort_to_parent=TRUE/FALSE # for cells not assigned to classes 
# # (leafs or user-defined)
# output_logical=TRUE/FALSE  # ignore overlap/resort_to_parent and just output a logical
# # matrix. If a single class is supplied, pipe it into "drop".


