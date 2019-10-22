if(FALSE) {
    library( "RUnit" )
    library( "metagene2" )
}

test.metagene_heatmap <- function() {
    test = get_demo_metagene()
    
    checkTrue(is(metagene2_heatmap(test), "ggplot"))
    checkTrue(is(metagene2_heatmap(test, coverage_order(test)), "ggplot"))
    checkTrue(is(metagene2_heatmap(test, coverage_order(test, "align2_rep2")), "ggplot"))
}

