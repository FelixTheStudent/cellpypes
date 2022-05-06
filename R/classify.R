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
#' @param invert. If TRUE, return classes that are NOT leafs instead of leaf nodes.
#'
#' @return Character vector with the tree's leafs.
#' @keywords internal
#'
tree_leaf_nodes <- function(classes, invert=FALSE) {
  stopifnot(is_classes(classes))
  
  res <- classes
  res$hasChild <- FALSE
  res$hasChild[res$class %in% res$parent] <- TRUE
  
  if(invert) return(res$class[res$hasChild] )
  return( res$class[!res$hasChild] )
}




#' Find parent, parent's parent and so on for a class using recursive programming
#'
#' @param classes The class definitions of a cellpypes object, i.e. obj$classes.
#' @param class A character vector of length one with the class.
#'
#' @return Character vector with the ancestry of a class.
#' @keywords internal
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

#' Find child, child's child and so on for class(es) using recursive programming
#'
#' @param classes The class definitions of a cellpypes object, i.e. obj$classes.
#' @param class A character vector with one or multiple classes. 
#' @param leafs Has to be the output of tree_leaf_nodes(classes). It's passed
#' as argument so that tree_leaf_nodes is not executed in each recursion.
#'
#' @return Character vector with the descendants of a class.
#' @keywords internal
#'
tree_descendants <- function(classes, class, leafs) {

  current_children <- classes[classes$parent %in% class,"class"]
  is_leaf <- current_children %in% leafs
  if(all(current_children %in% leafs)) {
    return(current_children)
  } else {
    return(c(current_children[is_leaf],
             current_children[!is_leaf],
             tree_descendants(classes, current_children[!is_leaf], leafs)))
    
  }
  

}








#





# obj example to develop functions here:
# obj <- simulated_umis %>%
#   rule("B", "MS4A1", ">", .1e-4)                     %>%
#   rule("T", "CD3E", ">", .1e-4)                      %>%
#   rule("Treg", "FOXP3", ">", .1e-4, parent="T")      %>%
#   rule("Tother", "FOXP3", "<", .1e-4, parent="T")    %>%
#   rule("actTreg", "ICOS", ">", .1e-4, parent="Treg")





#' Classify cells on previously defined rules
#' @template param_obj
#' @template classify_params
#' 
#'
#' @return A factor with cell type labels.
#' @template cellpypes_obj
#' @template handling_overlap
#' @export
#'
#' @examples
#' classify(rule(simulated_umis, "Tcell", "CD3E", ">", 1))
classify <- function(
  obj,
  classes=NULL,
    # If NULL, uses all childless classes (leafs).
    # unique(obj$classes$class) returns both leafs and parents.
  replace_overlap_with="Unassigned", # alternatives: 'common_parent', NA or any scalar character
  return_logical_matrix =FALSE # ignore overlap/unassigned rules and just output 
  # a logical matrix. If a single class is supplied, the matrix has exactly one
  # column and the user can pipe it into "drop" to convert it to a vector.
) {
  # checks specific for classify: classes have to be present in obj$classes
  # other sanity checks:
  check_obj(obj) 
  
  # factors give unexpected behaviour:
  if(is.factor(classes)) stop("Argument classes should be character, not factor. Use as.character!")

  if (is.null(classes)) {
    # Convention: use tree leaves if classes=NULL 
    classes <- tree_leaf_nodes(obj$classes)
    # when using leaves, rules from all classes are relevant for classification:
    relevant_classes <- obj$classes$class
  } else {
    if("Unassigned" %in% classes) stop("Don't pass 'Unassigned' to the 'classes' argument.")
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
  rules_eval <- mapply(FUN=evaluate_rule,
                       MoreArgs=list(obj=obj),
                       feature=rules_info$feature,
                       operator=rules_info$operator,
                       threshold=rules_info$threshold)
  # mapply sets (non-unique) features as colnames, preclude confusion with NULL:
  colnames(rules_eval) <- NULL

  res <- sapply(relevant_classes, function(class) {
    # We can use info to subset eval because their columns directly correspond:
    is_class_rule <- rules_info$class %in% tree_ancestry(obj$classes, class)
    sum(is_class_rule) == rowSums(rules_eval[, is_class_rule, drop=FALSE])
  })
  
    
  if (return_logical_matrix) {
    return(res[, classes, drop=FALSE])
    }
 
  # massage according to replace_overlap_with, etc.
  
  if (replace_overlap_with=="common_parent" || is.na(replace_overlap_with))
    stop("Not implemented yet, sorry.")
  # Note to myself: for "common_parent" res[, classes] is not enough.
 
   
  class_res <- res[, classes, drop=FALSE]
  class_factor <- rep("Unassigned", nrow(class_res))
  the_leafs <- tree_leaf_nodes(obj$classes) 
  for (i in 1:ncol(class_res)) {
    # If a class has overlap with one of its ancestor classes, 
    # I set the ancestor to FALSE here.
    # This way, classify replaces overlap EXCEPT for overlap with ones ancestor.
    descendants <- tree_descendants(obj$classes, classes[i], leafs=the_leafs)
    # Not all descendants are necessarily in class_res:
    descendants <- descendants[descendants %in% colnames(class_res)]
    cell_is_descendant <- base::rowSums(class_res[,descendants, drop=F]) > 0
    class_res[cell_is_descendant,i] <- FALSE
    # assign class:
    class_factor[ class_res[,i] ] <- classes[i]
  }
  # replace overlap, except for overlap with ones ancestors.
  class_factor[base::rowSums(class_res)>1] <- replace_overlap_with
  
  # Unassigned class may or may not be present, under different names:
  #  1. convert to factor at very end (adding new class to factor is hard)
  #  2. It is a design choice whether to add classes with 0 cells or not:
  class_factor <- factor( class_factor,
                          levels=unique(c(classes, "Unassigned", replace_overlap_with)))
  return( class_factor )
}





# I think the order of factor levels is the same as classes, and
# if is.null(classes) then it's the same as in obj$classes$class, i.e. the
# order in which classes were created by the user. I am not sure, however,
# but realize this is important (see e.g. f1_all.Rmd, building an f1-data.frame
# relies on it).



