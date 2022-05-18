




#' Plot the last modified rule or class
#'
#' @template param_obj
#' @param show_feat If TRUE (default), a second panel shows the feature plot of
#' the relevant gene.
#' @param what Either "rule" or "class".
#' @template param_fast
#' @param legend_rel_width Relative width compared to the other two plots
#' (only relevant if \code{show_feat=TRUE}). 
#'
#' @return Returns a ggplot2 object with the plot.
#' 
#' @template cellpypes_obj
#' 
#' @export
#' 
#' @importFrom ggplot2 ggplot aes_string geom_point coord_fixed xlab ylab ggtitle
#' @importFrom ggplot2 scale_color_manual theme_bw element_text margin
#' @importFrom cowplot plot_grid get_legend
#'
#' @examples
#' plot_last(rule(simulated_umis, "T", "CD3E",">", 1))
plot_last <- function(obj, show_feat=TRUE, what="rule", fast=NULL,
                      legend_rel_width=0.3) {
  check_obj(obj)
  if(is.null(fast)) fast <- ifelse(ncol(obj$raw)>10e3, TRUE, FALSE)
  stopifnot(is.logical(fast) && length(fast)==1)
  
  
  if(what=="rule") {
    last_rule <- obj$rules[nrow(obj$rules),]
    boolean=evaluate_rule(obj      = obj,
                          feature  = last_rule$feature,
                          operator = last_rule$operator, 
                          threshold=last_rule$threshold)
    plot_title <- paste0(last_rule$feature, " ", last_rule$operator, " ",
                         1e4*last_rule$threshold, " CP10K")
  } else if (what=="class") {
    last_class <- obj$rules[nrow(obj$rules), "class"]
    boolean=drop(classify(obj, classes=last_class, return_logical_matrix = T))
    plot_title <- paste0("Class: ", last_class)
  } else { stop("plot_last argument 'what' should either be rule or class.") }
  
  p <- ggplot(data=data.frame(V1=obj$embed[,1, drop=T], # tbl makes drop necessary
                         V2=obj$embed[,2, drop=T],
                         last=boolean),
         aes_string("V1", "V2", col="last"))+coord_fixed()+
    xlab( colnames(obj$embed)[1] ) + 
    ylab( colnames(obj$embed)[2] ) +
    ggtitle(plot_title) + 
    scale_color_manual(name="Rule",
                       values = c("TRUE"="#44AA99", # cartoColors (colorblind friendly)
                                  "FALSE"="#888888")) +
    theme_bw()+
    theme(plot.title = element_text(color="#44AA99"), legend.position = "none")
  if(fast) {p <- p+scattermore::geom_scattermore()} else {p <- p + geom_point()}
    
 
  # For saving etc. it is convenient to return the plot directly. 
  # I wanted plotting as side-effect, so with return(invisible(obj)),
  # to enable pipes like rule %>% plot_last %>% rule, but I learned
  # from Sveta that a T-pipe can do this anyways.
  if(show_feat&what=="rule") {
    pfeat <- feat(obj, last_rule$feature, fast=fast)
    legend <- cowplot::get_legend(pfeat +
                                    # create some space to the right of the legend
                                    theme(legend.box.margin = margin(0, 12, 0, 0))
    )
    assembled <- cowplot::plot_grid(
      p, pfeat+theme(legend.position = "none"), legend,
      ncol=3, rel_widths = c(1,1,legend_rel_width)
    )
    return(assembled) 
    # with patchwork:   p+feat(obj, last_rule$feature)
  }
  return(p)
  
}





#' Call and visualize 'classify' function
#'
#' @template param_obj
#' @param point_size Dot size used by \link[ggplot2]{geom_point}.
#' @param point_size_legend Dot size displayed in legend. 
#' Legend colors are easier to read with larger points.
#' @param base_size The base_size of \link[ggplot2]{theme_bw}, i.e. 
#' how large text is displayed. Default: 15.
#' @template classify_params
#' @template param_fast
#'
#' @return A ggplot2 object.
#' @template cellpypes_obj
#' @template handling_overlap
#' 
#' @export
#' 
#' @importFrom dplyr bind_cols
#'
#' @examples
#' plot_classes(rule(simulated_umis, "T", "CD3E",">", 1))
plot_classes <- function(obj,
                         classes=NULL,
                         replace_overlap_with="Unassigned", 
                         return_logical_matrix =FALSE,
                         fast = NULL,
                         point_size=.4,
                         point_size_legend=2,
                         base_size=15) {
  check_obj(obj)
  if(is.null(fast)) fast <- ifelse(ncol(obj$raw)>10e3, TRUE, FALSE)
  stopifnot(is.logical(fast) && length(fast)==1)
  
  labels <- classify(obj,
                     classes=classes, 
                     replace_overlap_with=replace_overlap_with, 
                     return_logical_matrix=return_logical_matrix) 
  if(is.logical(labels)) stop("Please set return_logical_matrix to FALSE.")
  # manually set Unassigned color:
  colors <- scales::hue_pal()(length(levels(labels))-1) 
  names(colors) <- levels(labels)[levels(labels)!="Unassigned"]
  colors <- c(colors, Unassigned="#888888")
  # do the plot
  plot_dat = bind_cols(
    dim1=obj$embed[,1, drop=TRUE],
    dim2=obj$embed[,2, drop=TRUE],
    class=labels)
  axis_names <- if (is.null(colnames(obj$embed))) {
    c("Dim1", "Dim2")
  } else {
    colnames(obj$embed)
  }
  p <- ggplot(plot_dat, aes_string("dim1", "dim2", col="class"))+
    coord_fixed()+
    scale_color_manual(values=colors) + 
    xlab(axis_names[1]) +
    ylab(axis_names[2]) +
    theme_bw(base_size = base_size) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = point_size_legend)))
  
  if(fast) { 
    p+scattermore::geom_scattermore(pointsize = point_size)
  } else {
    p + geom_point(size=point_size) 
  }
    
}











#' Feature plots: Color gene expression in 2D embeddings
#' 
#' Highlight gene expression in UMAP embeddings, for example.
#'
#' @template param_obj
#' @param features A vector of genes (features) to colour by.
#' @template param_fast
#' @param verbose feat ignores gene names not present in your object and 
#' warns you about them by default. `verbose`=FALSE will suppress the warning
#' (not recommended in interactive use).
#' @param ... Arguments passed to cowplot's \link[cowplot]{plot_grid} function,
#' for example ncol or rel_widths.
#' 
#'
#' @return A ggplot object (assembled by cowplot).
#' @template cellpypes_obj
#' 
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point coord_fixed xlab ylab ggtitle theme
#' @importFrom ggplot2 element_blank element_rect element_text
#' @examples
#' feat(simulated_umis, "CD3E")
feat <- function(obj, features, fast=NULL, verbose=TRUE, ...) {
  check_obj(obj)

  if(is.null(fast)) fast <- ifelse(ncol(obj$raw)>10e3, TRUE, FALSE)
  # User may pass features individuall, not as vector, by accident.
  # This happens to myself once per day, so I want an intelligible error message:
  if(fast %in% rownames(obj$raw)) stop(paste0(
    "Make sure to pass features as vector, not individually.\n",
    "This error appears because the fast is not logical but instead a gene in your object." )
    )
  stopifnot(is.logical(fast) && length(fast)==1)
  
  # old code from before I had the argument 'fast':
  #plotgrid_args <- list(...)
  #if(any(plotgrid_args %in% rownames(obj$raw))) {
  #  warning(paste0("Make sure to pass features as vector, not individually.\n",
  #  "This warning appears because arguments in ... are genes in your object." )
  #  )}
 
  axis_names <- if (is.null(colnames(obj$embed))) {
    c("Dim1", "Dim2")
  } else {
    colnames(obj$embed)
  }
  
  # user might enter duplicated or nonsensical feature names:
  features <- unique(features)
  does_exist<- features %in% rownames(obj$raw)
  if(sum(does_exist)==0) stop("None of the supplied features were found in your object.") 
  if(verbose && (sum(does_exist) < length(does_exist))) {
    cat("These features were not found and will be ignored: ",
        "\n",
        paste(features[!does_exist], sep=", "),
        "\n")
  }
  features <- features[does_exist]
  

  
  
  if (is.null(obj$totalUMI) & inherits(obj$raw, "Matrix")) { 
    obj$totalUMI <- Matrix::colSums(obj$raw)
  }  
  if (is.null(obj$totalUMI) & inherits(obj$raw, "matrix")) { 
    obj$totalUMI <- base::colSums(obj$raw)
  }  
  
  l = lapply(features, function(gene) {

    dat <- data.frame(obj$embed,
                     obj$raw[gene,],
                     obj$totalUMI)
    colnames(dat) <- c("X1","X2","k","s")
    
    dat$expr <- 1e4*dat$k/dat$s # UMI per ten thousand (per10k)
    # Replacing zeros with a tenth of minimum works well I found:
    if(all(dat$k==0)) {
      dat$expr <- .01
    } else {
      dat$expr <- ifelse(dat$expr > 0,
                         dat$expr,
                         min(dat$expr[dat$expr>0])/10)
    }
    p <- ggplot(dat,
           aes_string(x = "X1", y = "X2", col = "expr")) +
      coord_fixed() +
      ggtitle(gene) + xlab(axis_names[1]) + ylab(axis_names[2]) +
      viridis::scale_color_viridis(
        name = "CP10K",
        trans="log2",
        breaks=scUtils::closed_breaks_log2,
        labels= function(br) scUtils::closed_labels(
          br, 
          min_is_zero = any(dat$k==0))
      ) +
      theme_bw()
    
    if(fast) {
      p+scattermore::geom_scattermore()
    } else {
      p + geom_point(size=.4)
    }
      # declutter plots:
      # theme(axis.ticks = element_blank(),
      #       axis.text = element_blank(),
      #       panel.background = element_rect(fill="white", color = "black"),
      #       panel.grid = element_blank())
      
  })
  
  if(length(l) == 1) return(l[[1]])   # single plot need not be list

  cowplot::plot_grid(plotlist = l, ...) 
}
