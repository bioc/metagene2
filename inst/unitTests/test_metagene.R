## Test functions present in the metagene.R file

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene2" )
    library( "data.table" )
    library( "dplyr" )
}

### }}}

bam_files <- get_demo_bam_files()
named_bam_files <- bam_files
names(named_bam_files) <- letters[1:(length(named_bam_files))]
not_indexed_bam_file <- metagene2:::get_not_indexed_bam_file()
regions <- metagene2:::get_demo_regions()
design <- data.frame(Samples = c("align1_rep1.bam", "align1_rep2.bam",
                    "align2_rep1.bam", "align2_rep2.bam", "ctrl.bam"),
                    align1 = c(1,1,0,0,2), align2 = c(0,0,1,1,2))
design$Samples <- paste0(system.file("extdata", package = "metagene2"), "/",
                        design$Samples)
set.seed(1)
demo_mg <- metagene2$new(regions = get_demo_regions(),
                        bam_files = get_demo_bam_files())
                        
full_mg = demo_mg$clone(deep=TRUE)
full_mg$produce_metagene()

util_test_invalid_param_constructor_value <- function(param_name, param_value, error_value) {
    arg_list_new = list(bam_files=get_demo_bam_files(), 
                        regions=get_demo_regions())
    arg_list_new[[param_name]] = param_value
    
    obs <- tryCatch(do.call(metagene2$new, arg_list_new),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
}

util_test_invalid_param_value <- function(param_name, param_value, step_function, error_value) {
    arg_list_new = list(bam_files=get_demo_bam_files(), 
                        regions=get_demo_regions())
    arg_list_new[[param_name]] = param_value
    
    obs <- tryCatch(do.call(metagene2$new, arg_list_new),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
    
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(do.call(mg[[step_function]], arg_list_new[param_name]),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
    
    mg <- demo_mg$clone(deep=TRUE)
    obs <- tryCatch(do.call(mg$produce_metagene, arg_list_new[param_name]),
                    error = conditionMessage)
    checkIdentical(obs, error_value)      
}
                        
#region <- regions[1]
#bam_file <- bam_files[1]
#demo_mg_min <- metagene2$new(regions = region, bam_files = bam_file)

# Load fake SAMs for value tests.
fake_bam_files = file.path(system.file("extdata", package = "metagene2"),
                           c("fake_align1.bam", "fake_align2.bam", "fake_align3.bam"))
# Define a design that will test various combinations of bam files.
fake_bam_design = data.frame(BAM=fake_bam_files,
                             fake_align1=c(1, 0, 0),
                             fake_align2=c(0, 1, 0),
                             fake_align3=c(0, 0, 1),
                             fake_align12=c(1, 1, 0),
                             fake_align13=c(1, 0, 1),
                             fake_align23=c(0, 1, 1),
                             fake_align123=c(1, 1, 1),
                             stringsAsFactors=FALSE)

# Define the region over which single-region tests will be run.
# fake_align1 has 1 read over the whole region.
# fake_align2 has 2 reads over the whole region.
# fake_align3 has 4 reads over the first half of the region, 8 over the second half.
fake_bam_region_1_chr="chr1"
fake_bam_region_1_pos_start=1000000
fake_bam_region_1_pos_end = 1004999
fake_bam_region_1_pos_range=fake_bam_region_1_pos_start:fake_bam_region_1_pos_end
fake_bam_region_1_test_region_unique = GRanges(paste0(fake_bam_region_1_chr, ":", fake_bam_region_1_pos_start, "-", fake_bam_region_1_pos_end))

# Load the expected genomic coverages generated along with the bam files.
fake_bam_expected_coverages = list()
fake_bam_expected_rpm = list()
for(i in fake_bam_design$BAM) {
    cov_file = gsub(".bam", ".sam", i)
    load(paste0(cov_file, ".coverage.RData"))
    fake_bam_expected_coverages[[i]] = cov_rle
    
    load(paste0(cov_file, ".coverage_rpm.RData"))
    fake_bam_expected_rpm[[i]] = rpm_rle
}

###################################################
## Test the metagene2$new() function (initialize)
###################################################

## Invalid verbose value
test.metagene_invalid_verbose <- function() {
    util_test_invalid_param_constructor_value("verbose", "ZOMBIES", "verbose must be a logicial value (TRUE or FALSE)")
}

## Invalid force_seqlevels value
test.metagene_invalid_force_seqlevels <- function() {
    util_test_invalid_param_constructor_value("force_seqlevels", "ZOMBIES", "force_seqlevels must be a logicial value (TRUE or FALSE)")
}

## Negative padding_size value
test.metagene_invalid_padding <- function() {
    util_test_invalid_param_constructor_value("padding_size", "ZOMBIES", "padding_size must be a non-negative integer")
    util_test_invalid_param_constructor_value("padding_size", -1, "padding_size must be a non-negative integer")
    util_test_invalid_param_constructor_value("padding_size", 1.2, "padding_size must be a non-negative integer")
}


## Negative core value
test.metagene_invalid_core <- function() {
    util_test_invalid_param_constructor_value("cores", "ZOMBIES", "cores must be a positive integer or a BiocParallelParam instance.")
    util_test_invalid_param_constructor_value("cores", -1, "cores must be a positive integer or a BiocParallelParam instance.")
    util_test_invalid_param_constructor_value("cores", 1.2, "cores must be a positive integer or a BiocParallelParam instance.")
    util_test_invalid_param_constructor_value("cores", 0, "cores must be a positive integer or a BiocParallelParam instance.")    
}

## Non-character vector bam_files value
util_test_invalid_bam_file <- function(value, error_value) {
    obs <- tryCatch(metagene2:::metagene2$new(regions = get_demo_regions(),
                                            bam_files = value),
                    error = conditionMessage)
    checkIdentical(obs, error_value)
}

test.metagene_invalid_bam_files <- function() {
    util_test_invalid_bam_file(c(2,4,3), "bam_files must be a vector of BAM filenames.")
    util_test_invalid_bam_file(list(a = "a.txt", b = "b.txt"), "bam_files must be a vector of BAM filenames.")
    util_test_invalid_bam_file(not_indexed_bam_file, "All BAM files must be indexed")
    util_test_invalid_bam_file(c(bam_files, not_indexed_bam_file), "All BAM files must be indexed")
}

util_test_invalid_region <- function(invalid_region, error, ...) {
    obs <- tryCatch(do.call(metagene2$new, c(list(...), 
                                             list(bam_files=get_demo_bam_files(),
                                                  regions = invalid_region))),
                    error = conditionMessage)
    checkIdentical(obs, error)
}

util_test_valid_region <- function(valid_region, ...) {
    obs <- do.call(metagene2$new, c(list(...), 
                                    list(bam_files=get_demo_bam_files(),
                                         regions = valid_region)))
    util_test_valid_metagene(obs)
}

util_test_valid_metagene <- function(mg) {
    checkIdentical(class(mg), c("metagene2", "R6"))
}

# regions is the wrong class
test.metagene_initialize_invalid_region_value <- function() {
    region <- array(data = NA, dim = c(2,2,2))
    util_test_invalid_region(region, paste0("regions must be either a vector of BED filenames, a ",
                                       "GRanges object or a GrangesList object"))
}

test.metagene_regions_seqlevels <- function() {
    # GRanges have seqlevels (indicating all possible seqnames)
    # and seqnames themselves. Extra seqlevels are no concerns.
    region_with_extra_seq_level = get_demo_regions()[[1]]
    GenomeInfoDb::seqlevels(region_with_extra_seq_level) <- c(GenomeInfoDb::seqlevels(region_with_extra_seq_level),
                                                              "extra_seqlevels")
    
    util_test_invalid_region(region_with_extra_seq_level, "Some seqlevels of regions are absent in bam_file")
    util_test_valid_region(region_with_extra_seq_level, force_seqlevels = TRUE)

    # Extra seqnames means we must drop certain regions. We only do so
    # if force_seqlevels=TRUE, otherwise we throw an error.
    region_with_extra_seq = region_with_extra_seq_level
    seqnames(region_with_extra_seq)[1] = "extra_seqlevels"

    util_test_invalid_region(region_with_extra_seq, "Some seqlevels of regions are absent in bam_file")
    util_test_valid_region(region_with_extra_seq, force_seqlevels = TRUE)

    # Sometimes there can be no seqnames left after removing
    # those with unknown levels. This happens often in chromosome names
    # are mismatched ("chr1" vs "1")
    region_no_common_seq = region_with_extra_seq_level
    seqnames(region_no_common_seq) = factor(rep("extra_seqlevels", 50), levels=c("chr1", "extra_seqlevels"))
    
    util_test_invalid_region(region_no_common_seq, "No seqlevels matching between regions and bam file", force_seqlevels=TRUE)
}


# Valid regions narrowPeak
test.metagene_initialize_valid_narrowpeak <- function() {
    region <- metagene2:::get_narrowpeak_file()
    mg <- metagene2$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
                   qValue = "numeric", peak = "integer")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    exp$region_name="list1"
    checkIdentical(obs, exp)
}

# Valid regions broadPeak
test.metagene_initialize_valid_broadpeak <- function() {
    region <- metagene2:::get_broadpeak_file()
    mg <- metagene2$new(regions = region, bam_files = bam_files[1])
    obs <- mg$get_regions()
    extraCols <- c(signalValue = "numeric", pValue = "numeric",
                   qValue = "numeric")
    exp <- rtracklayer::import(region, format = "BED", extraCols = extraCols)
    exp$region_name="list1"
    checkIdentical(obs, exp)
}

# Valid named bam files
test.metagene_initialize_valid_bam_files <- function() {
    mg <- metagene2$new(regions = get_demo_regions(), bam_files = get_demo_bam_files())
    
    # Make sure all bam_files were kept.
    obs <- mg$get_params()[["bam_files"]]
    exp <- get_demo_bam_files()
    checkTrue(all(obs == exp))
    
    # Make sure we have coverage for all bam files.
    obs <- names(mg$get_raw_coverages())
    exp <- get_demo_bam_files()
    checkTrue(all(obs == exp))
    
    # NAmed bam files should have no impact (Historically, they changed the coverage names.)
    named_bam_files = get_demo_bam_files()
    names(named_bam_files) = letters[seq_along(named_bam_files)]
    mg <- metagene2$new(regions = get_demo_regions(), bam_files = named_bam_files)
    
    # Make sure we have coverage for all bam files under the correct names.
    obs <- names(mg$get_raw_coverages())
    exp <- unname(named_bam_files)
    checkIdentical(obs, exp)    
}

##################################################
# Test the metagene2$get_params() function
##################################################

## Valid usage
test.metagene_get_params_valid_usage <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    params <- mg$get_params()
    checkIdentical(unname(params[["bam_files"]]), get_demo_bam_files())
    checkIdentical(params[["padding_size"]], 0)
    checkIdentical(params[["verbose"]], FALSE)
    checkIdentical(params[["force_seqlevels"]], FALSE)
}

##################################################
# Test the metagene2$get_design() function
##################################################

## Valid usage
test.metagene_get_design_valid_usage <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    mg$group_coverages(design=get_demo_design())
    design <- mg$get_design()
    checkIdentical(mg$get_design()[,-1], get_demo_design()[,-1])
    checkTrue(all(mg$get_design()[,1] == mg$get_params()[["bam_files"]]))
}

##################################################
# Test the metagene2$get_regions() function
##################################################

## Valid usage default
test.metagene_get_regions_valid_usage_default <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    regions <- mg$get_regions()
    checkIdentical(length(regions), length(unlist(get_demo_regions())))
    region_1_in_regions = length(get_demo_regions()[[1]])
    region_1_in_metagene = sum(regions$region_name==names(get_demo_regions())[1])
    checkIdentical(region_1_in_regions, region_1_in_metagene)
}

##################################################
# Test the metagene2$get_matrice() function
##################################################

# Tests to write:
#  Test single region
#  Test bin_size = 1
#  Test bin_count < width(regions)

test.metagene_group_coverages_valid_usage_default = function(){
    mg <- demo_mg$clone(deep=TRUE)
    grouped_coverages = mg$group_coverages(design=get_demo_design())

    checkIdentical(length(grouped_coverages), ncol(get_demo_design()) - 1L)
}

# Make sure bam_files order does not change results.
# See issue #7: https://github.com/ArnaudDroitLab/metagene2/issues/7
test.metagene_group_coverages_bam_files_order = function(){
    mg  = metagene2$new(regions=get_demo_regions(), bam_files=get_demo_bam_files())
    mg2 = metagene2$new(regions=get_demo_regions(), bam_files=rev(get_demo_bam_files()))
    checkIdentical(mg$bin_coverages()[["align1_rep1"]][1:10,1:10],
                   mg2$bin_coverages()[["align1_rep1"]][1:10,1:10])
}


test.metagene_bin_coverages_valid_usage_default_design = function(){
    mg <- demo_mg$clone(deep=TRUE)
    mg$group_coverages()
    binned_coverages = mg$bin_coverages()

    design = mg$get_design()
    checkIdentical(length(binned_coverages), ncol(design) - 1L)
    checkIdentical(names(binned_coverages), colnames(design)[-1])
    for(i in 1:length(binned_coverages)) {
        checkIdentical(nrow(binned_coverages[[i]]), length(mg$get_regions()))
        checkTrue(ncol(binned_coverages[[i]]) == mg$get_params()[["bin_count"]])
    }
}

test.metagene_bin_coverages_valid_usage_custom_design = function(){
    mg <- demo_mg$clone(deep=TRUE)
    mg$group_coverages(design=get_demo_design())
    binned_coverages = mg$bin_coverages()

    checkIdentical(length(binned_coverages), ncol(get_demo_design()) - 1L)
    checkIdentical(names(binned_coverages), colnames(get_demo_design())[-1])
    for(i in 1:length(binned_coverages)) {
        checkIdentical(nrow(binned_coverages[[i]]), length(mg$get_regions()))
        checkTrue(ncol(binned_coverages[[i]]) == mg$get_params()[["bin_count"]])
    }
}

test.metagene_split_coverages_valid_usage_default = function(){
    # We'll populate a list with split_coverages results which should all be identical,
    # then run validations against those list elements.
    split_coverages_list = list()
    
    # Case 1: Use default split_by, which will default to region_name.
    #         So we'll pass a GRangesList upon initialization to get 3 region groups.
    split_regions_default = GRangesList(
        Just4=GRanges(data.frame(seqnames="chr1", start=c(1000001, 1000101), end=c(1000100, 1000200))),
        Just8=GRanges(data.frame(seqnames="chr1", start=c(1002501, 1002601), end=c(1002600, 1002700))),
        Mixed=GRanges(data.frame(seqnames="chr1", start=c(1000001, 1002501), end=c(1000100, 1002600))))

    mg_default <- metagene2$new(bam_files=fake_bam_design$BAM, regions=split_regions_default)
    split_coverages_list[["default"]] =  mg_default$split_coverages_by_regions()
        
    # Case 2: Use custom split_by. So we'll pass a GRanges where the regions have a "Custom" attribute,
    # which we'll split on.
    split_regions_custom = unlist(split_regions_default)
    split_regions_custom$Custom = c("Just4", "Just4", "Just8", "Just8", "Mixed", "Mixed")
    
    mg_custom <- metagene2$new(bam_files=fake_bam_design$BAM, regions=split_regions_custom, split_by="Custom")
    split_coverages_list[["custom"]] =  mg_custom$split_coverages_by_regions()
    
    # Finally, validate both split_coverages_by_regions results are good.
    for(split_coverages in split_coverages_list) {
        # Make sure that we've got matrix list per bam file/design group.
        checkTrue(length(split_coverages) == length(fake_bam_design$BAM))
        
        # Make sure each matrix list has has many matrices as there are region groups.
        for(i in 1:length(split_coverages)) {
            for(j in 1:length(split_regions_default)) {
                checkTrue(nrow(split_coverages[[i]][[j]]) == length(split_regions_default[[j]]))
            }
        }
        
        # Make sure numeric values within the matrices are as expected.
        checkTrue(all(split_coverages[["fake_align_3"]]$Just4 == 4))
        checkTrue(all(split_coverages[["fake_align_3"]]$Just8 == 8))
        checkTrue(all(split_coverages[["fake_align_3"]]$Mixed[,1] == 4))
        checkTrue(all(split_coverages[["fake_align_3"]]$Mixed[,2] == 8))
    }
}


test.metagene_raw_coverage_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=fake_bam_region_1_test_region_unique)
    
    # Test genome-wide raw and normalized coverages for
    # strand_mode=FALSE
    raw_coverages = mg$get_raw_coverages()
    norm_coverages = mg$get_normalized_coverages()
    for(i in fake_bam_design$BAM) {
        obs = raw_coverages[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        exp = fake_bam_expected_coverages[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        checkIdentical(obs, exp)
                       
        obs = norm_coverages[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        exp = fake_bam_expected_rpm[[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        checkTrue(all(abs(obs - exp) < 10e-8))
    }
}

test.metagene_group_coverage_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=fake_bam_region_1_test_region_unique)
     
    # Test group coverages.
    group_coverages = mg$group_coverages(design=fake_bam_design, simplify=FALSE)
    
    # Grouped coverages for unit designs should be the same as the plain ones.
    expected_names = file.path(system.file("extdata", package = "metagene2"), 
                               paste0(c("fake_align1", "fake_align2", "fake_align3"), ".bam"))
    names(expected_names) = c("fake_align1", "fake_align2", "fake_align3")
    
    for(i in c("fake_align1", "fake_align2", "fake_align3")) {
        obs = group_coverages[["*"]][[i]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        exp = fake_bam_expected_coverages[[expected_names[i]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
        checkIdentical(obs, exp)
    }

    # Test combined group coverages.
    obs = group_coverages[["*"]][["fake_align12"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[[expected_names["fake_align1"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[[expected_names["fake_align2"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)
    
    obs = group_coverages[["*"]][["fake_align23"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[[expected_names["fake_align2"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[[expected_names["fake_align3"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)
    
    obs = group_coverages[["*"]][["fake_align13"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[[expected_names["fake_align1"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[[expected_names["fake_align3"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)

    obs = group_coverages[["*"]][["fake_align123"]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    exp = fake_bam_expected_coverages[[expected_names["fake_align1"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[[expected_names["fake_align2"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range] +
          fake_bam_expected_coverages[[expected_names["fake_align3"]]][[fake_bam_region_1_chr]][fake_bam_region_1_pos_range]
    checkIdentical(obs, exp)
}

test.metagene_bin_coverage_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=fake_bam_region_1_test_region_unique, design=fake_bam_design)
    
    # Test binned coverages for bins 
    # where all values are the same.
    bin_coverages = mg$bin_coverages()
    checkTrue(all(bin_coverages[["fake_align1"]][1,]==1))
    checkTrue(all(bin_coverages[["fake_align2"]][1,]==2))
    checkTrue(all(bin_coverages[["fake_align3"]][1,1:50]==4))
    checkTrue(all(bin_coverages[["fake_align3"]][1,51:100]==8))
    checkTrue(all(bin_coverages[["fake_align123"]][1,1:50]==7))
    checkTrue(all(bin_coverages[["fake_align123"]][1,51:100]==11))
    
    # Test binned coverages where some values are means.
    bin_coverages = mg$bin_coverages(bin_count=5)
    checkTrue(all(bin_coverages[["fake_align1"]][1,]==1))
    checkTrue(all(bin_coverages[["fake_align2"]][1,]==2))
    checkTrue(bin_coverages[["fake_align3"]][1,1]==4)
    checkTrue(bin_coverages[["fake_align3"]][1,2]==4)
    checkTrue(bin_coverages[["fake_align3"]][1,3]==6)
    checkTrue(bin_coverages[["fake_align3"]][1,4]==8)
    checkTrue(bin_coverages[["fake_align3"]][1,5]==8)
}

test.metagene_calculate_ci_values_unique_region = function(){
    # Create the metagene object.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, regions=fake_bam_region_1_test_region_unique, design=fake_bam_design, bin_count=5)
 
    full_df = mg$add_metadata()

    # Only one region, all confidence intervals should be NA.
    checkTrue(all(is.na(full_df$qinf) & is.na(full_df$qsup)))
    
    # There should be five bins.
    checkTrue(max(full_df$bin)==5 && length(unique(full_df$bin))==5)
    
    # Bin values should match
    checkTrue((full_df %>% filter(bin==3 & design=="fake_align3") %>% pull(value)) == 6)
    
    # All groups should be represented.
    checkTrue(all(full_df$design %in% colnames(fake_bam_design)[-1]))
    
    test_plot = mg$produce_metagene()
    checkTrue("ggproto" %in% class(test_plot$scales))
}


##################################################
# Test the metagene2$get_data_frame() function
##################################################

## Valid usage default
test.metagene_get_data_frame_valid_usage_default <- function() {
 df <- full_mg$get_data_frame()
 regions <- get_demo_regions()
 bam_files <- get_demo_bam_files()
 checkTrue(is.data.frame(df))
 checkTrue(ncol(df) == 8)
 checkTrue(nrow(df) == length(regions) * length(bam_files) * full_mg$get_params()$bin_count)
}

### Valid usage subset
test.metagene_get_data_frame_valid_usage_subset <- function() {
 single_region = names(get_demo_regions())[1]
 three_designs = colnames(full_mg$get_design())[2:4]

 df <- full_mg$get_data_frame(region_names = single_region, design_names =three_designs )

 checkTrue(is.data.frame(df))
 checkTrue(ncol(df) == 8)
 checkTrue(nrow(df) == length(single_region) * length(three_designs) * full_mg$get_params()$bin_count)
}

## Valid usage get_data_frame return by copy of data_frame
test.metagene_get_data_frame_check_copy_of_data_frame <- function() {
 mg <- full_mg$clone(deep=TRUE)
 
 df1 <- mg$get_data_frame()
 #modification of table by reference
 df1$c <- 1:1000

 #Is table copied and unchanged ? 
 df2 <- mg$get_data_frame()
 checkIdentical(ncol(df1) == ncol(df2), FALSE)
}

## Valid usage no data_frame produced
test.metagene_get_data_frame_valid_usage_no_data_frame <- function() {
 mg <- demo_mg$clone(deep=TRUE)
 df <- mg$get_data_frame()
 checkTrue(is.null(df))
 df_subset <- mg$get_data_frame(get_demo_regions()[1],
                                get_demo_bam_files()[1:2])
 checkTrue(is.null(df_subset))
}

##################################################
# Test the metagene2$get_plot() function
##################################################

## Valid case no graph
test.metagene_get_plot_valid_case_no_graph <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    plot <- mg$get_plot()
    checkTrue(is.null(plot))
}

## Valid case graph
test.metagene_get_plot_valid_case_graph <- function() {
    plot <- full_mg$get_plot()
    checkTrue(all(class(plot) == c("gg", "ggplot")))
}

##################################################
# Test the metagene2$get_raw_coverages() function
##################################################

exp_raw <- GenomicAlignments::readGAlignments(bam_files[1])
exp_raw <- GenomicAlignments::coverage(exp_raw)

## Default filenames
test.metagene_get_raw_coverages <- function() {
    obs <- full_mg$get_raw_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_raw[[i]]),
                        logical(1))))
}

##################################################
# Test the metagene2$get_normalized_coverages() function
##################################################

count <- Rsamtools::countBam(bam_files[1])$records
weight <- 1 / (count / 1000000)
exp_norm <- exp_raw * weight

## Default filenames
test.metagene_get_normalized_coverages <- function() {
    mg <- demo_mg$clone(deep=TRUE)
    obs <- mg$get_normalized_coverages()[[1]]
    checkTrue(all(vapply(1:length(obs),
                        function(i) all(obs[[i]]==exp_norm[[i]]),
                        logical(1))))
}

# Invalid bin_count class
test.metagene_invalid_bin_count <- function() {
  util_test_invalid_param_value("bin_count", "a", "bin_coverages", "bin_count must be a positive integer")
  util_test_invalid_param_value("bin_count", -1, "bin_coverages", "bin_count must be a positive integer")
  util_test_invalid_param_value("bin_count", 1.2, "bin_coverages", "bin_count must be a positive integer")  
}

# Invalid normalization class
test.metagene_invalid_normalization <- function() {
 util_test_invalid_param_value("normalization", 1234, "group_coverages", 'normalization must be NULL, "RPM" or "NCIS".')
 util_test_invalid_param_value("normalization", "CSI", "group_coverages", 'normalization must be NULL, "RPM" or "NCIS".')
}
##################################################
# Test the metagene2$produce_data_frame() function
##################################################

### Invalid alpha value
test.metagene_invalid_alpha <- function(){
    util_test_invalid_param_value("alpha", 'test', "calculate_ci", "is.numeric(alpha) is not TRUE")
    util_test_invalid_param_value("alpha", -0.8, "calculate_ci", "alpha >= 0 & alpha <= 1 is not TRUE")    
}

### Invalid sample_count class
test.metagene_invalid_sample_count <- function(){
    util_test_invalid_param_value("sample_count", 'test', "calculate_ci", "is.numeric(sample_count) is not TRUE")
    util_test_invalid_param_value("sample_count", -10, "calculate_ci", "sample_count >= 0 is not TRUE")    
}

test.metagene_replace_metadata <- function() {
    regions_gr <- unlist(get_demo_regions())
    demo_metadata = data.frame(BedName=names(regions_gr),
                               EvenStart=ifelse((start(regions_gr) %% 2) == 0, "Even", "Odd"),
                               Strand=strand(regions_gr))

    mg <- metagene2$new(regions = get_demo_regions(),
                       region_metadata=demo_metadata,
                       bam_files = get_demo_bam_files(),
                       assay='chipseq')
    
    mg$produce_metagene(split_by=c("EvenStart", "Strand"))
    test_meta = mg$get_regions_metadata()
    
    # Add column and replace. This should work fine.
    test_meta$Foo = c("Foo", "Bar")
    mg$replace_region_metadata(test_meta)

    # Remove region_name column. It should be added back.
    old_region_name = mg$get_regions_metadata()$region_name
    test_meta = test_meta[setdiff(colnames(test_meta), "region_name")]
    mg$replace_region_metadata(test_meta)
    
    checkTrue(all(mg$get_regions_metadata()$region_name==old_region_name))
    
    # Remove one of the split_by columns. Caches should be invalidated.
    test_meta = test_meta[setdiff(colnames(test_meta), "EvenStart")]
    mg$replace_region_metadata(test_meta)    
    checkIdentical(mg$get_params()[["split_by"]], "region_name")
    
    # Try using a data-frame without enough rows. We should get an error message.
    obs <- tryCatch(mg$replace_region_metadata(test_meta[1:4,]),
                    error = conditionMessage)
    exp <- "region_metadata must have one row per region."
    checkIdentical(obs, exp)    
    
}

test.metagene_strand_specific <- function() {
    plus_file = file.path(system.file("extdata", package = "metagene2"), "fake_align1.bam")
    minus_file = file.path(system.file("extdata", package = "metagene2"), "fake_align1_minus.bam")
    
    test_region_list=GRangesList(Plus=GRanges("chr1:1000001-1000100:+"),
                                 Minus=GRanges("chr1:1001001-1001100:-"))
    
    mg = metagene2$new(bam_files=c(Plus=plus_file, Minus=minus_file),
                       regions=test_region_list,
                       strand_specific=TRUE)
    
    # Make sure raw coverages are strand-specific.
    checkTrue(all(mg$get_raw_coverages()[["+"]][[plus_file]][test_region_list[["Plus"]]] == 1))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[plus_file]][test_region_list[["Plus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["+"]][[minus_file]][test_region_list[["Plus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[minus_file]][test_region_list[["Plus"]]] == 0))
    
    checkTrue(all(mg$get_raw_coverages()[["+"]][[plus_file]][test_region_list[["Minus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[plus_file]][test_region_list[["Minus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["+"]][[minus_file]][test_region_list[["Minus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[minus_file]][test_region_list[["Minus"]]] == 1))

    # Make sure split coverages are strand-specific
    split_matrices = mg$split_coverages_by_regions()
    checkTrue(all(split_matrices[["Plus"]][["Plus"]] == 1))
    checkTrue(all(split_matrices[["Plus"]][["Minus"]] == 0))
    checkTrue(all(split_matrices[["Minus"]][["Plus"]] == 0))
    checkTrue(all(split_matrices[["Minus"]][["Minus"]] == 1))
    
    # Now test with inverted strand.
    mg = metagene2$new(bam_files=c(Plus=plus_file, Minus=minus_file),
                       regions=test_region_list, 
                       strand_specific=TRUE, invert_strand=TRUE)
    
    # Make sure raw coverages are strand-specific.
    checkTrue(all(mg$get_raw_coverages()[["+"]][[plus_file]][test_region_list[["Plus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[plus_file]][test_region_list[["Plus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["+"]][[minus_file]][test_region_list[["Plus"]]] == 1))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[minus_file]][test_region_list[["Plus"]]] == 0))
    
    checkTrue(all(mg$get_raw_coverages()[["+"]][[plus_file]][test_region_list[["Minus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[plus_file]][test_region_list[["Minus"]]] == 1))
    checkTrue(all(mg$get_raw_coverages()[["+"]][[minus_file]][test_region_list[["Minus"]]] == 0))
    checkTrue(all(mg$get_raw_coverages()[["-"]][[minus_file]][test_region_list[["Minus"]]] == 0))

    # Make sure split coverages are strand-specific
    split_matrices = mg$split_coverages_by_regions()
    checkTrue(all(split_matrices[["Plus"]][["Plus"]] == 0))
    checkTrue(all(split_matrices[["Plus"]][["Minus"]] == 1))
    checkTrue(all(split_matrices[["Minus"]][["Plus"]] == 1))
    checkTrue(all(split_matrices[["Minus"]][["Minus"]] == 0))    
}

test.metagene_extend_reads_basic <- function() {
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, 
                        regions=fake_bam_region_1_test_region_unique, 
                        design=fake_bam_design,
                        extend_reads=100)
    
    # No reads overlap in fake_align1. The reads are all 50nt, 
    # on the + strand, and the coverage is 1 everywhere.
    # So we should have a coverage of 1 for the first 50 nucleotides,
    # of 2 for the rest. Coverage should extend on the + stran
    checkTrue(all(mg$get_raw_coverages()[[fake_bam_design$BAM[1]]][["chr1"]][1000000:1000049] == 1))
    checkTrue(all(mg$get_raw_coverages()[[fake_bam_design$BAM[1]]][["chr1"]][1000050:1002500] == 2))
}

test.metagene_extend_reads_over_edge <- function() {    
    # Some reads might not overlap the specified regions, yet extend into them.
    # Those must be calculated into the coverage.
    edge_region = flank(fake_bam_region_1_test_region_unique, width=50, start=FALSE)
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, 
                        regions=edge_region, 
                        design=fake_bam_design,
                        extend_reads=100, bin_count=25)
    checkTrue(all(mg$get_raw_coverages()[[fake_bam_design$BAM[1]]][["chr1"]][start(edge_region):end(edge_region)] == 1))
}

test.metagene_extend_reads_directionality <- function() {                          
    # Reads should be extended only in their specified direction.
    # So if we go upstream from the first read in the bam, coverage should be 0.
    edge_region = flank(fake_bam_region_1_test_region_unique, width=50, start=TRUE)
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, 
                        regions=edge_region, 
                        design=fake_bam_design,
                        extend_reads=100)
    checkTrue(all(mg$get_raw_coverages()[[fake_bam_design$BAM[1]]][["chr1"]][start(edge_region):end(edge_region)] == 0))
}

test.metagene_extend_reads_no_shrink <- function() {     
    # Reads should not be shrunk when expand_reads is smaller than read size.
    # mg <- metagene2$new(bam_files=fake_bam_design$BAM, 
    #                     regions=fake_bam_region_1_test_region_unique, 
    #                     design=fake_bam_design,
    #                     extend_reads=35)
    # checkTrue(all(mg$get_raw_coverages()[[fake_bam_design$BAM[1]]][["chr1"]][1000000:1002500] == 1))  
    TRUE
}

test.metagene_extend_reads_beyond_chr_end <- function() {     
    # Expanding past the end of chromosomes (Or past the start) should
    # not raise any errors.
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, 
                    regions=fake_bam_region_1_test_region_unique, 
                    design=fake_bam_design,
                    extend_reads=195471971, bin_count=100)
    mg$produce_metagene()
}

test.metagene_discontiguous <- function() {
    fake_genes = GRangesList(
        Single=        GRanges(c("chr1:1000001-1000100:+")),
        Multiple=      GRanges(c("chr1:1000001-1000100:+",
                                 "chr1:1000201-1000300:+")),
        MultipleValues=GRanges(c("chr1:1000001-1000100:+",
                                 "chr1:1002501-1002600:+")))
                               
    mg <- metagene2$new(bam_files=fake_bam_design$BAM, 
                    regions=fake_genes, 
                    design=fake_bam_design,
                    region_mode="stitch")
    mg$produce_metagene()
    
    # Make sure each design has binned coverages for each gene.
    split_coverages = mg$split_coverages_by_regions()
    for(design_coverage in split_coverages) {
        checkTrue(length(design_coverage) == length(fake_genes))
    }
    
    # Check numerical values of binned coverages.
    checkTrue(all(split_coverages[["fake_align1"]][["Single"]] == 1))
    checkTrue(all(split_coverages[["fake_align1"]][["Multiple"]] == 1))
    checkTrue(all(split_coverages[["fake_align1"]][["MultipleValues"]] == 1))

    checkTrue(all(split_coverages[["fake_align3"]][["Single"]] == 4))
    checkTrue(all(split_coverages[["fake_align3"]][["Multiple"]] == 4))
    checkTrue(all(split_coverages[["fake_align3"]][["MultipleValues"]][1,1:50] == 4))
    checkTrue(all(split_coverages[["fake_align3"]][["MultipleValues"]][1,51:100] == 8))
    
    # Try grouping genes using region_metadata.
    region_metadata = data.frame(CustomGroup=c("A", "A", "B"))
    mg$replace_region_metadata(region_metadata)
    mg$produce_metagene(split_by="CustomGroup")
    
    split_coverages = mg$split_coverages_by_regions()
    checkTrue(all(split_coverages[["fake_align1"]][["A"]] == 1))
    checkTrue(all(split_coverages[["fake_align1"]][["B"]] == 1))
    
    checkTrue(all(split_coverages[["fake_align3"]][["A"]] == 4))
    checkTrue(all(split_coverages[["fake_align3"]][["B"]][1,1:50] == 4))
    checkTrue(all(split_coverages[["fake_align3"]][["B"]][1,51:100] == 8))    
}

test.metagene_bin_count_larger_than_regions <- function() {
    mg <- metagene2$new(bam_files=get_demo_bam_files(),
                        regions=get_demo_regions(),
                        bin_count=1000000)
                    
    # Test for contiguous regions.                    
    obs = tryCatch(mg$bin_coverages(), error = conditionMessage)
    smallest_region = min(width(unlist(get_demo_regions())))
    exp = paste0("The specified bin_count (", 1000000, "nt) ",
              "is larger than the smallest effective region (", smallest_region, "nt).\n",
              "Please make sure bin_count is >= ", smallest_region, "\n")
    checkTrue(obs==exp)
   
    # Test for discontiguous regions.
    region_test = GRangesList(
        Test1=GRanges(c("chr1:1000001-1000100:+",
                        "chr1:1000301-1000400:+")),
        Test2=GRanges(c("chr1:1000001-1000200:+",
                        "chr1:1000301-1000500:+")))
    mg <- metagene2$new(bam_files=get_demo_bam_files(),
                        regions=region_test,
                        bin_count=200, region_mode="stitch")   
    
    # Individual region minimum is 100, but the sum of stitched regions is 200,
    # so a 200 bin_count should not fail.
    mg$bin_coverages()
    
    # However, 1000000 should fail.
    obs = tryCatch(mg$bin_coverages(bin_count=1000000), error = conditionMessage)
    smallest_region = min(unlist(lapply(region_test, function(x) { sum(width(x)) })))
    exp = paste0("The specified bin_count (", 1000000, "nt) ",
             "is larger than the smallest effective region (", smallest_region, "nt).\n",
             "Please make sure bin_count is >= ", smallest_region, "\n")
    checkTrue(obs==exp)            
}

test.metagene_invalid_extend_paired_end <- function() {
    obs = tryCatch(metagene2$new(bam_files=get_demo_bam_files(),
                                 regions=get_demo_regions(),
                                 paired_end=TRUE, extend_reads=200),
                   error=conditionMessage)
    exp = "extend_reads and paired_end cannot both be set at the same time."
    checkTrue(obs==exp)

}