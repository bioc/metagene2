#' Produce a metagene plot
#' 
#' @param df a \code{data.frame} obtained with the \code{get_data_frame}
#' function. Must have the following columns: "region", "design", "bin", 
#' "value", "qinf" and "qsup".
#'
#' @param facet_by A formula to be used for facetting the metagene plot. This
#'                 formula can include any design metadata, or region_metadata 
#.                 columns that were part of the \code{split_by} argument.
#'                 \code{NA} can be used to keep the previous value.
#'                 Default: \code{NA}.
#' @param group_by A string representing a single column from design_metadata or
#'                 region_metadata which will be used to group observations 
#'                 together into lines and which will be used to generate the 
#'                 color scale. \code{NA} can be used to keep the previous 
#'                 value. Default: \code{NA}.
#'
#' @return A `ggplot` object.
#'
#' @import ggplot2
#' @examples
#' mg <- get_demo_metagene()
#' df <- mg$add_metadata()
#' p <- metagene2:::plot_metagene(df)
plot_metagene <- function(df, facet_by=NULL, group_by=NULL) {
    df$design <- as.factor(df$design)

    expected_cols <- c("bin", "value", "qinf", "qsup", "group")
    assert_subset<-df[,which(colnames(df) %in% expected_cols)]
    expected_class <- c("integer", rep("numeric", 3), "factor")
    names(expected_class) = expected_cols
    stopifnot(all(expected_cols %in% colnames(assert_subset)))
    actual_classes = vapply(assert_subset, class, character(1))
    actual_classes = actual_classes[expected_cols]
    stopifnot(all(actual_classes == expected_class))

    if(is.null(group_by)) {
        group_by="group"
    }
    
    p <- ggplot(df, aes(x=bin, y=value, ymin=qinf, ymax=qsup)) +
        geom_ribbon(aes_string(fill=group_by), alpha=0.3) +
        geom_line(aes_string(color=group_by), size=1) +
        theme(panel.grid.major = element_line()) +
        theme(panel.grid.minor = element_line()) +
        theme(panel.background = element_blank()) +
        theme(panel.background = element_rect()) +
        theme_bw(base_size = 20)
        
    if(!is.null(facet_by)) {
        p <- p + facet_grid(facet_by)
    }
    
    return(p)
}