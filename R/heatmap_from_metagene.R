# Take the nested list of (design group) -> (region_group) coverages
# returned by split_coverages_by_regions and melts it into a data-frame
# where each region/bin/design_group is an individual row.
combine_coverages = function(cov_matrices, region_order) {
    df_list = list()
    for(bam in names(cov_matrices)) {
        for(region_name in names(cov_matrices[[bam]])) {
            input_matrix = cov_matrices[[bam]][[region_name]][region_order[[region_name]],]
            label = paste(bam, region_name)
            df_list[[label]] = combine_list_coverages(input_matrix, bam, region_name)
        }
    }
    
    do.call(rbind, df_list)
}

# Melts a single binned coverage matrix and melts it, adding a bam and region 
# group column.
#' @importFrom reshape2 melt
combine_list_coverages = function(x, bam, region) {
    res = reshape2::melt(x, varnames=c("region", "bin"), value.name="coverage")
    res$bam_name=bam
    res$region_name=region
    
    res
}

#' Returns an "as-is" ordering of regions.
#'
#' This function creates an ordering of regions to be used with the 
#' \code{\link{metagene2_heatmap}} function. The regions are not actually
#' reordered, but returned as-is.
#'
#' @param metagene The metagene object whose grouped regions should be ordered.
#' @return A list, with as many elements as there are region groups in the
#'         metagene object. Each element of that list is an ordering of the
#'         regions of that group based on their original ordering in the 
#'         metagene2 object.
#' @examples
#'   demo_metagene = get_demo_metagene()
#'   as_is_region_order(demo_metagene)
#' @export
as_is_region_order <- function(metagene) {
    nrow_order = function(x) { seq_len(nrow(x)) }
    lapply(metagene$split_coverages_by_regions()[[1]], nrow_order)
}

# Parallel sum, an analog to pmax.
psum = function(...) {
    mat = do.call(cbind, list(...))
    apply(mat, 1, sum)
}

#' Determines ordering of regions as a function of coverage.
#'
#' This function creates an ordering of regions within region groups based on 
#' ascending or descending mean coverage. This is used with the 
#' \code{\link{metagene2_heatmap}} function.
#'
#' @param metagene The metagene object whose grouped regions should be ordered.
#' @param design_groups A vector of design groups to be used for determining the 
#'                      ordering. If \code{NULL}, all design groups are used.
#' @param decreasing If TRUE, regions are ordered from the highest mean coverage
#'                   to the lowest mean coverage, and vice versa.
#' @return A list, with as many elements as there are region groups in the
#'         metagene object. Each element of that list is an ordering of the
#'         regions of that group based on their mean coverage.
#' @examples
#'   demo_metagene = get_demo_metagene()
#'   coverage_order(demo_metagene)
#' @importFrom purrr pmap
#' @export
coverage_order <- function(metagene, design_groups=NULL, decreasing=TRUE) {
    # Get the coverages, split by bam and region groups.
    cov_matrices = metagene$split_coverages_by_regions()
    
    # Get the coverage means for all individual regions.
    means = lapply(cov_matrices, function(x) { lapply(x, rowMeans) })
    
    # If a subset of design groups was not specified, choose all design groups.
    if(is.null(design_groups)) {
        design_groups = names(means)
    }
    
    # Sum the mean coverages over specified design groups and order the results.
    lapply(purrr::pmap(means[design_groups], psum), order, decreasing=decreasing)
}

#' Plots a heatmap of coverages from a metagene2 object.
#'
#' This function creates an ordering of regions within region groups based on 
#' ascending or descending mean coverage. This is used with the 
#' \code{\link{metagene2_heatmap}} function.
#'
#' @param metagene The metagene object to be plotted as a heatmap.
#' @param region_order A named list with as many elements as there are region 
#'                     groups, with each element containing an ordering for the 
#'                     regions within that group. 
#'                     The \code{\link{as_is_region_order}} and 
#'                     \code{\link{coverage_order}} functions can be used to
#'                     generate a valid ordering. By default, 
#'                     \code{\link{as_is_region_order}} is used.
#' @param scale_trans A character string giving the transformation that should
#'                    be applied to the coverage values. Common values are
#'                    "identity" and "log1p". See the ggplot2 documentation for
#'                    scale_continuous for more details.
#' @return A ggplot object containing a heatmap representation of the metagene2 
#'         object.
#' @examples
#'   demo_metagene = get_demo_metagene()
#'   metagene2_heatmap(demo_metagene)
#' @import ggplot2
#' @export
metagene2_heatmap <- function(metagene, region_order = as_is_region_order(metagene), scale_trans="identity") {
    stopifnot(is(metagene, "metagene2"))
    stopifnot(is(region_order, "list"))
    # Make sure the names of the region groups and ordering match.
    if(any(sort(names(region_order)) != sort(names(metagene$get_split_regions())))) {
        stop("The names of the provided ordering does not match those of the ",
             "region groups.")
    }
    
    # Make sure the lengths of the region groups and ordering match.
    for(i in names(region_order)) {
        if(length(region_order[[i]]) != length(metagene$get_split_regions()[[i]])) {
           stop("The length of the region ordering does not match the length of regions.")
        }
    }
    
    # Combine the coverages into a data-frame suitable for ggplot.
    heatmap_df = combine_coverages(metagene$split_coverages_by_regions(),
                                   region_order)

    # Build the plot.
    ggplot(heatmap_df) +
        geom_raster(mapping=aes_string(x="bin", y="region", fill="coverage")) + 
        facet_grid(region_name~bam_name) +
        scale_y_continuous(trans = "reverse") +
        scale_fill_gradient(low="#FFFFFF", high="#18448c", trans=scale_trans) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              strip.background = element_blank(),
              panel.border = element_rect(fill=NA, colour = "black"))
}
