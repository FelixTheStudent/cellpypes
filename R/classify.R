# Functions that return classification outcomes (class assignments)



# classes encodes a phylogenic tree of classes/cell types.
# These are the functions to work with that:
#  tree_ancestry     function(tree, class)   returns class and its ancestry
#  tree_leafs        function(tree)          returns names of all leaf nodes   
#  tree_plot         function(obj)           ggplot object or something :)




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
class_boolean <- function(obj, classes, drop=TRUE) {
  stopifnot(is.character(classes))
  stopifnot(is_classes(obj$classes))
  stopifnot(is_rules(obj$rules))
  

  boolean_matrix <- sapply(classes, function(class) { # loop through classes
    ancestry <- tree_ancestry(obj$classes, class)
    class_rules <- obj$rules[obj$rules$class %in% ancestry,]
    # Reduce can be thought of here as looping through class_rules:
    base::Reduce(function(x,y) x&y, 
                lapply(1:nrow(class_rules), function(i){
                  evaluate_rule(obj      = obj,
                                class    = class_rules[i, "class"],
                                feature  = class_rules[i, "feature"],
                                operator = class_rules[i, "operator"],
                                threshold= class_rules[i, "threshold"]
                                ) })
                )
  })
  
  if (ncol(boolean_matrix)==1 && drop) {
    drop(boolean_matrix)
  } else{
    boolean_matrix
  }
     
}






# obj example to develop functions here:
# obj <- simulated_umis %>%
#   rule("B", "MS4A1", ">", .1e-4)                     %>%
#   rule("T", "CD3E", ">", .1e-4)                      %>%
#   rule("Treg", "FOXP3", ">", .1e-4, parent="T")      %>%
#   rule("Tother", "FOXP3", "<", .1e-4, parent="T")    %>%
#   rule("actTreg", "ICOS", ">", .1e-4, parent="Treg")





classify <- function(
  obj,
  classes,
    # If NULL, uses all childless classes (leafs).
    # unique(obj$classes$class) returns both leafs and parents.
  replace_overlap_with="Unassigned", # alternatives: 'common_parent', NA or any scalar character
  replace_unassigned_with="Unassigned", # "common_parent", NA or any scalar character
  return_logical_matrix =FALSE # ignore overlap/resort_to_parent and just output a logical
# matrix. If a single class is supplied, pipe it into "drop".
) {
  # checks specific for classify: classes have to be present in obj$classes
  # other sanity checks:
  check_obj(obj) 
  


  # all operations
  relevant_classes <- lapply(classes, tree_ancestry, classes=obj$classes)
  relevant_classes <- unique(do.call(c, relevant_classes))
  
  
  # evaluate all rules in obj$rules to obtain a
  # boolean matrix (see class_boolean function). 
  # Doing this step first ensures this computation is not repeated for
  # ancestor classes (parents of multiple leafs).
  # Its columns correspond directly to class names in obj$rule$class,
  # so you can use obj$rule$class to extract all rules for class "T", for example.
  
  # select classes you want to return. Could be all leaves, or user-defined.
  
  # for a given class, subset columns in boolean matrix with
  #  obj$rule$class %in% tree_ancestry(obj, class) and
  # combine to obtain logical vector per class.
  
  # massage according to replace_overlap_with, replace_unassigned_with and
  # return_logical
  
  # return
}

# Useful to develop code:
# obj <- simulated_umis %>%
#   rule("T", "CD3E", ">", .1e-3) %>%
#   rule("Treg", "FOXP3", ">", .01e-3, parent="T") %>%
#   rule("Treg_act", "ICOS", ">", .01e-3, parent="Treg") %>%
#   rule("B", "MS4A1", ">", .1e-3)
# obj %>% classify(classes=c("Treg_act", "B"))

# Alternative parametrization, which I think is stupid (delete soon):
# classes = c("Treg", "Ttox", "B"). 
# 
# overlap = c("common_parent", NA or any character string: "Unassigned"/"multiplelabels"/...)
# resort_to_parent=TRUE/FALSE # for cells not assigned to classes 
# # (leafs or user-defined)
# output_logical=TRUE/FALSE  # ignore overlap/resort_to_parent and just output a logical
# # matrix. If a single class is supplied, pipe it into "drop".


