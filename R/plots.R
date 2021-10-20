




#' Plot the last modified rule or class
#'
#' @param obj A cellpypes object.
#' @param show_feat If TRUE (default), a second panel shows the feature plot of
#' the relevant gene.
#' @param what Either "rule" or "class".
#'
#' @return The plot is drawn as side-effect with `print`. The function
#' returns obj invisibly, so that you can use it in a pipe.
#' 
#' 
#' @export
#' 
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed xlab ylab ggtitle
#' @importFrom ggplot2 scale_color_manual theme_bw
#'
#' @examples
plot_last <- function(obj, show_feat=TRUE, what="rule") {
  check_obj(obj)
  if(what=="rule") {
    last_rule <- obj$rules[nrow(obj$rules),]
    boolean=evaluate_rule(obj,
                          last_rule$class, last_rule$feature,
                          last_rule$operator,   last_rule$threshold)
    plot_title <- paste0("Rule: ", last_rule$feature, last_rule$operator, 
                         last_rule$threshold)
  } else if (what=="class") {
    last_class <- obj$rules[nrow(obj$rules), "class"]
    boolean=drop(classify(obj, classes=last_class, return_logical_matrix = T))
    plot_title <- paste0("Class: ", last_class)
  } else { stop("plot_last argument 'what' should either be rule or class.") }
  
  p <- ggplot(data=data.frame(V1=obj$embed[,1, drop=T], # tbl makes drop necessary
                         V2=obj$embed[,2, drop=T],
                         last=boolean),
         aes(V1, V2, col=last))+coord_fixed()+
    geom_point()  +
    xlab( colnames(obj$embed)[1] ) + 
    ylab( colnames(obj$embed)[2] ) +
    ggtitle(plot_title) + 
    scale_color_manual(values = c("TRUE"="#44AA99", # cartoColors (colorblind friendly)
                                  "FALSE"="#888888")) +
    theme_bw()
  
  # Functions are either transforming (rule) or side-effects (plot_last). 
  # Side-effect functions return the obj so that you can use them in pipes.
  if(show_feat) print(p+feat(obj, last_rule$feature)) else print(p)
  return( invisible(obj) ) # enables this: obj %>% plot_last() %>% rule(...)
}





#' Call and visualize 'classify' function
#'
#' @param obj 
#' @param ... Same parameters as \link[cellpypes]{classify}.
#'
#' @return
#' @export
#'
#' @examples
plot_classes <- function(obj, ...) {
  check_obj(obj)
  
  labels <- classify(obj, ...)
  # manually set Unassigned color:
  colors <- scales::hue_pal()(length(unique(labels))-1) 
  names(colors) <- unique(labels)[unique(labels)!="Unassigned"]
  colors <- c(colors, Unassigned="#888888")
  # do the plot
  plot_dat = bind_cols(
    dim1=obj$embed[,1],
    dim2=obj$embed[,2],
    class=labels)
  axis_names <- if (is.null(colnames(obj$embed))) {
    c("Dim1", "Dim2")
  } else {
    colnames(obj$embed)
  }
  ggplot(plot_dat, aes(dim1, dim2, col=class))+
    coord_fixed()+
    geom_point(size=.4) +
    scale_color_manual(values=colors) + 
    xlab(axis_names[1]) +
    ylab(axis_names[2]) +
    theme_bw()
  
}






#' Title
#'
#' @param obj 
#' @param feature_name 
#'
#' @return
#'
#' @examples
feat_data <- function(obj, feature_name) {
  if (is.null(obj$totalUMI)) { 
    obj$totalUMI <- Matrix::rowSums(obj$raw)
  } 
  df <- data.frame(obj$embed,
                   obj$raw[,feature_name],
                   obj$totalUMI)
  colnames(df) <- c("X1","X2","k","s")
  
  return(df)
}






#' Feature plots: Color gene expression in 2D embeddings
#'
#' @param obj 
#' @param features 
#' @param ... 
#'
#' @return
#' @export
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed xlab ylab ggtitle theme
#' @importFrom ggplot2 element_blank element_rect
#' @examples
feat <- function(obj, features, ...) {
  check_obj(obj)
  
  axis_names <- if (is.null(colnames(obj$embed))) {
    c("Dim1", "Dim2")
  } else {
    colnames(obj$embed)
  }
  
  # user might enter duplicated or nonsensical feature names:
  features <- unique(features)
  does_exist<- features %in% colnames(obj$raw)
  if(sum(does_exist)==0) stop("None of the supplied features were found in your object.") 
  if(sum(does_exist) < length(does_exist)) {
    cat("These features were not found and will be ignored: ",
        "\n",
        paste(features[!does_exist], sep=", "),
        "\n")
  }
  features <- features[does_exist]
  
  l = lapply(features, function(gene) {
    dat <- feat_data(obj, gene)
    dat$expr <- 1e4*dat$k/dat$s # UMI per ten thousand (per10k)
    # Replacing zeros with a tenth of minimum works well I found:
    if(all(dat$k==0)) {
      dat$expr <- .01
    } else {
      dat$expr <- ifelse(dat$expr > 0,
                         dat$expr,
                         min(dat$expr[dat$expr>0])/10)
    }
    ggplot(dat,
           aes(x = X1, y = X2, col = expr)) +
      coord_fixed() +
      geom_point(size=.4) +
      ggtitle(gene) + xlab(axis_names[1]) + ylab(axis_names[2]) +
      viridis::scale_color_viridis(
        name = "per10k",
        trans="log2",
        breaks=scUtils::closed_breaks_log2,
        labels= function(br) scUtils::closed_labels(
          br, 
          min_is_zero = any(dat$k==0))
      ) +
      # declutter plots:
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.background = element_rect(fill="white", color = "black"),
            panel.grid = element_blank())
      
  })
  

  cowplot::plot_grid(plotlist = l, ...) 
}
