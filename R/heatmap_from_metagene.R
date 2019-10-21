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

combine_list_coverages = function(x, bam, region) {
    res = melt(x, varnames=c("region", "bin"), value.name="coverage")
    res$bam_name=bam
    res$region_name=region
    
    res
}

as_is_region_order <- function(metagene) {
    lapply(metagene$split_coverages_by_regions()[[1]], function(x) { seq_len(nrow(x)) })
}

psum = function(...) {
    mat = do.call(cbind, list(...))
    apply(mat, 1, sum)
}

coverage_order <- function(metagene, bam_names=NULL, decreasing=TRUE) {
    cov_matrices = metagene$split_coverages_by_regions()
    means = lapply(cov_matrices, function(x) { lapply(x, rowMeans) })
    if(is.null(bam_names)) {
        bam_names = names(means)
    }
    
    lapply(purrr::pmap(means, psum), order, decreasing=decreasing)
}

# region_order:
#   as_is: Don't touch
#   
metagene2_heatmap <- function(metagene, region_order = as_is_region_order(metagene)) {
    heatmap_df = combine_coverages(metagene$split_coverages_by_regions(),
                                   region_order)

    ggplot(heatmap_df) +
        geom_tile(mapping=aes(x=bin, y=region, fill=coverage)) + 
        facet_grid(region_name~bam_name) +
        scale_y_continuous(trans = "reverse") +
        scale_fill_gradient(low="#FFFFFF", high="#18448c") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              strip.background = element_blank(),
              panel.border = element_rect(fill=NA, colour = "black"))
}
    
    
# test = get_demo_metagene()
# 
# metagene2_heatmap(test)
# metagene2_heatmap(test, coverage_order(test))
# metagene2_heatmap(test, coverage_order(test, "align2_rep2"))



