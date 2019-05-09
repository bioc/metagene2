#' A class to manage metagene analysis.
#'
#' This metagene2 class encapsulates all of the steps necessary to perform
#' metagene analyses, which are aggregations of coverages over multiple regions 
#' (genes) to reveal patterns that might not be apparent from looking at
#' individual regions. It will allow to load, convert and normalize bam 
#' alignments and regions files/data. Once the data is ready, the user can then
#' choose to produce metagene plots on the data or some subset of it.
#'
#' Most metagene analyses are a two-step affair:
#' \enumerate{
#'  \item Initialize the object using \code{mg = metagene2$new()}, specifying 
#'        which regions and bam files should be used for the analysis.
#'  \item Generate a metagene plot using \code{mg$produce_metagene}, specifying 
#'        any additional parameter (Number of bins, facetting variables, etc).
#' }
#'
#' The \code{metagene2} object will then internally chain all 6 required  
#' processing steps, updating its internal caches along the way:
#' \enumerate{
#'  \item Coverages are inferred from bam files (\code{metagene2$new}).
#'  \item Coverages from multiple bam files are grouped and normalized 
#'       (\code{mg$group_coverages}).
#'  \item Coverages are binned together (\code{mg$bin_coverages})
#'  \item Binned coverages are split according to the type of region they belong 
#'        to (\code{mg$split_coverages_by_regions}).
#'  \item Coverage means and confidence intervals are calculated for each 
#'        region * group combination (\code{mg$calculate_ci}).
#'  \item Metadata is added to the calculated coverages ( 
#'        \code{mg$add_metadata}).
#'  \item The metagene is plotted (\code{mg$plot}).
#' }
#'
#' Each of these steps has an associated function, which takes as input certain 
#' parameters of the metagene analysis and returns an intermediary structure
#' of interest (coverages, binned coverages, long-form data frame of confidence 
#' intervals, etc). Those are described below, in the "Processing methods" 
#' section.
#'
#' All processing methods automatically call previous processing steps if
#' those have not already been run. For example, there is no need to call
#' \code{mg$group_coverages()} before calling \code{mg$bin_coverages()}: the
#' metagene2 object will automatically detect that certain prerequisite steps 
#' have not yet been performed, and run them.
#' 
#' Additionally, when calling produce_metagene a second time to change certain
#' analysis parameters after generating an initial metagene plot, only the 
#' required caches are reset: all non-impacted aspects of the analysis are left 
#' untouched, decreasing processing time.
#'
#' For further examples, see the metagene2 vignette.
#'
#' @section Constructor:
#'
#' \strong{Usage:}
#'
#'  \code{mg <- metagene2$new(regions, bam_files, padding_size = 0,
#'                            cores = SerialParam(), verbose = FALSE,
#'                            force_seqlevels = FALSE, paired_end = FALSE,
#'                            assay = 'chipseq', strand_specific=FALSE,
#'                            paired_end_strand_mode=2,
#'                            region_mode="auto", region_metadata=NULL, 
#'                            extend_reads=0, invert_strand=FALSE, ...)}
#'
#' \strong{Description:}
#'
#' This method returns a new \code{metagene2} object. Upon initialization, a  
#' \code{metagene2} object calculates coverages over all given regions in the 
#' provided bam files. Any and all parameter associated with any of the 
#' processing steps can be initialized upon object construction. All analysis  
#' parameters that are not explicitly specified in the constructor call are  
#' initialized to sensible defaults.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{regions}{A description of all regions over which metagenes will be calculated.
#'
#'                   When \code{region_mode} is "separate", those should be provided using
#'                   a \code{GRanges} object representing all individual, contiguous regions
#'                   to be examined.
#'
#'                   When \code{region_mode} is "stitch", those should be provided using a 
#'                   \code{GRangesList} object where each individual \code{GRanges} element
#'                   represents a set of regions to be stitched together.
#'
#'                   As a convenience, in "separate" mode, \code{metagene2} will convert any
#'                   passed in \code{GRangesList} into an unlisted \code{GRanges} with an 
#'                   additional \code{region_name} metadata column containing 
#'                   the name of the \code{GRangesList} element it was extracted 
#'                   from.
#'
#'                   Also as a convenience, \code{regions} can also be a \code{character} 
#'                   \code{vector} of filenames, which are then imported into a GRangesList.
#'                   Supported file formats are BED, narrowPeak, broadPeak, gff and gtf.}
#'    \item{bam_files}{A \code{vector} of BAM filenames. The BAM files must be
#'                    indexed. i.e.: if a file is named file.bam, there must
#'                    be a file named file.bam.bai or file.bai in the same 
#'                    directory. If \code{bam_files} is a named vector, then the provided names
#'                    can be used downstream to refer to those bam files. If no
#'                    names are provided, \code{metagene2} will try to infer appropriate ones.}
#'    \item{assay}{\code{'chipseq'}, \code{'rnaseq'} or NULL. If non-NULL, metagene will
#'                 set other parameters, such as region_mode and strand_specific, to logical
#'                 values for the given assay. Default: \code{'chipseq'}}
#'    \item{region_mode}{Set the way the \code{regions} parameter is interpreted. Can be 
#'                       \code{'separate'}, \code{'stitch'} or \code{'auto'}. In separate mode,
#'                       \code{regions} is expected to be a GRanges defining individual, contiguous regions. In
#'                       \code{'stitch'} mode, \code{regions} is expected to be a GRangesList 
#'                       where each \code{GRanges} element represents a set of regions to be 
#'                       stitched together and treated as a single logical region. If \code{'auto'}
#'                       then a logical value is inferred from the \code{assay} parameter.
#'                       Default: \code{'auto'}}
#'    \item{region_metadata}{A data-frame of metadata to be associated with the elements of \code{regions}.
#'                           It must contain has many rows as there are elements in \code{regions}. If 
#'                           \code{region_metadata} is NULL but \code{regions} has an mcols element,
#'                           then it is used.}
#'    \item{padding_size}{The provided \code{regions} will be extended on each side by the
#'                        value of this parameter. The padding_size must be a
#'                        non-negative integer. Default = 0.}
#'    \item{cores}{The number of cores available to parallelize the analysis.
#'                Either a positive integer or a \code{BiocParallelParam}.
#'                Default: \code{SerialParam()}.}
#'    \item{verbose}{Print progression of the analysis. A logical constant.
#'                    Default: \code{FALSE}.}
#'    \item{force_seqlevels}{If \code{TRUE}, remove regions that are not found
#'                in bam file header. Default: \code{FALSE}. \code{TRUE} and \code{FALSE}
#'                respectively correspond to pruning.mode = "coarse"
#'                and "error" in ?seqinfo.}
#'    \item{paired_end}{Set this to \code{TRUE} if the provided bam files
#'                      describe paired-end reads. If \code{FALSE}, single-ended
#'                      data are expected. Default: \code{FALSE}}
#'    \item{strand_specific}{If \code{TRUE}, only reads which align to the same 
#'                           strand as those specified in \code{regions} will
#'                           count toward coverage for that region. Useful for RNA-seq
#'                           profiles generated from strand-specific libraries, such
#'                           as Illumina TruSeq. Default: \code{'FALSE'}}
#'    \item{paired_end_strand_mode}{\code{'1'} or \code{'2'}. In paired-end mode,
#'                                  indicates which read in a pair sets the pair's strand.
#'                                  If \code{1}, this is the first read (This should be used
#'                                  with directional protocols such as Directional Illumina 
#'                                  (Ligation) or Standard SOLiD).
#'                                  If \code{2}, this is the second read (This should be used
#'                                  with directional protocols such as dUTP, NSR, NNSR, 
#'                                  or Illumina stranded TruSeq PE).
#'                                  Ignored if either paired_end or strand_specific is FALSE.
#'                                  Default: \code{'2'}}
#'    \item{extend_reads}{Extend individual reads to have a minimum length equal to this parameter.
#'                        When set to 0, no read extension occurs. This is useful for single-end
#'                        chip-seq experiments, where the length of the captured fragment is usually
#'                        longer than the sequenced read.}
#'    \item{invert_strand}{If \code{TRUE}, coverages for the given regions will be inferred 
#'                         from the coverage on the strand opposite theirs. Useful
#'                         for single-end stranded experiments which use cDNA.
#'                         This parameter is ignored if strand-specific is \code{FALSE}.}
#'    \item{...}{Additional parameters for the metagene analysis. See \code{produce_metagene}
#'               for a list of possible parameters.}
#' }
#'
#'    \code{metagene2$new} returns a \code{metagene2} object that contains the
#'        coverages for every BAM files in the regions from the \code{regions}
#'        parameter.
#'
#' @return
#' \code{metagene2$new} returns a \code{metagene2} object which contains the
#' normalized coverage values for every regions in all specified BAM files.
#'
#' @section produce_metagene():
#'
#' \strong{Usage:}
#'
#' \code{mg$produce_metagene(...)}
#'
#' \strong{Description:}
#'
#' \code{produce_metagene} is the workhorse method of the metagene2 object.
#' This method performs all of the necessary analysis steps for the production
#' of the metagene plot, and returns that plot. Any and all parameters of the 
#' metagene analysis, as documented in the individual processing steps, can be 
#' passed to \code{produce_metagene}. The metagene2 object will then determines
#' which intermediate caches would be affected by changes to those parameters,
#' invalidate them, and rerun all steps up to the plotting. This makes 
#' \code{produce_metagene} ideal for fast, iterative takes on the data.
#'
#' Below we present those parameters and a brief description of their usage.
#' Please refer to the affected processing step for a more in-depth explanation
#' of each parameter.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{design}{A \code{data.frame} that describes the grouping of the bam files
#'            into design groups. By default, each bam file is its own design group.
#'            See \code{group_coverages}.}
#'    \item{normalization}{The algorithm to use to normalize coverages,
#'                        \code{NULL} (no normalization), "RPM" or "NCIS". By default,
#'                        no normalization occurs. See \code{group_coverages}.}
#'    \item{design_filter}{Indices indicating which subset of design groups should 
#'                         be included in the analysis. By default, all design  
#'                         groups/bam files are included. See 
#'                         \code{group_coverages}.}
#'    \item{bin_count}{The number of bins regions should be split into. Defaults
#'                     to 100. See \code{bin_coverages}.}
#'    \item{region_filter}{The subset of regions to be kept for the analysis.
#'                         By default, all regions are kept. See \code{bin_coverages}}
#'    \item{split_by}{Which metadata columns should we use to split the set of 
#'                    regions into subset of interests? Defaults to "region_name",
#'                    an automatically added column. 
#'                    See \code{split_coverages_by_regions}.}
#'    \item{alpha}{The alpha level of the confidence interval estimates.
#'                 Defaults to 0.05. See \code{calculate_ci}.}
#'    \item{sample_count}{The number of draws to perform in the bootstrap
#'                        calculations used to calculate the confidence inteval.  
#'                        Defaults to 1000. See \code{calculate_ci}}
#'    \item{resampling_strategy}{The resampling strategy to be used when performing the
#'                               bootstrap analysis, which can be either \code{'profile'}
#'                               or \code{'bin'}. Defaults to \code{'bin'}. See
#'                               \code{calculate_ci}.}
#'    \item{design_metadata}{A data-frame containing metadata for the design groups.
#'                           By default, no metadata is associated. See 
#'                           \code{add_metadata}.}
#'    \item{title}{A title to add to the graph. See \code{plot}.}
#'    \item{x_label}{X-axis label for the metagene plot. See \code{plot}.}
#'    \item{facet_by}{A formula to be used for facetting the metagene plot. 
#'                    By default, no facetting occurs. See \code{plot}.}
#'    \item{group_by}{The metadata column used to build the color scale. By 
#'                    default, the combination of design and region name is 
#'                    used. See \code{plot}.}
#' }
#'
#' @section Processing methods:
#' 
#' Each of the following methods perform one step of metagene processing.
#' Most do not need to be called explicitly. Instead, you can simply call 
#' \code{produce_metagene}. However, you can use them to access intermediary
#' results: grouped coverages, binned coverages, split coverages, and long-form
#' data-frame of coverages with confidence intervals.
#'
#' @section group_coverages:
#'
#' \strong{Usage:}
#'
#'  \code{mg$group_coverages(design=NA, normalization=NA, 
#'            design_filter=NA, simplify=FALSE)}
#'
#' \strong{Description:}
#'
#'  This method normalizes genome-wide coverages, then groups
#'  them according to the specified design groups. It returns
#'  a list of possible read orientations (+, -, *), each element
#'  of which is either NULL (depending on the value of the 
#'  strand_specific parameter) or a list of possible design groups.
#'  In turn, the lists of design groups contain lists of \code{Rle}
#'  objects representing coverage over a specific chromosome or sequence.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{design}{A \code{data.frame} that describes the grouping of the bam files
#'            into design groups. The first column of the design should contain the 
#'            names of bam_files passed on initialization. Each subsequent columns
#'            represents a design group, that is to say a combination of bam files
#'            whose coverages should be grouped together into a logical unit.
#'            These columns should contain integer values indicating whether the 
#'            bam files on that row should be excluded (0), included as an`
#'            "input" (1) or included as a "control" (2) within the specified 
#'            design group.
#'            \code{NA} can be used keep previous design value. Default: \code{NA}.}
#'    \item{normalization}{The algorithm to use to normalize coverages. Possible
#'                        values are \code{NULL} (no normalization), "RPM", "log2_ratio"
#'                        and "NCIS". "RPM" transforms raw counts into Reads-Per-Million.
#'                        "log2_ratio" uses the formula log2((input RPM + 1) / (control RPM + 1))
#'                        to calculate a log-ratio between input and control. NCIS
#'                        attempts to subtract control from input. See
#'                        Liand and Keles 2012 for the NCIS algorithm. \code{NA} can 
#'                        be used keep the previous value. Default: \code{NA}}
#'    \item{design_filter}{A logical vector specifying which of the design groups specified
#'                         within the \code{design} parameter should be included in the metagene.
#'                         Useful for quickly reprocessing a subset of samples.
#'                         \code{NA} can be used keep previous design value. Default: \code{NA}}
#'    \item{simplify}{In single strand mode, set \code{simplify} to \code{TRUE} to return 
#'                    only the '*' coverage and omit the empty '+' and '-' components.
#'                    Default: \code{FALSE}}
#' }
#' @section bin_coverages:
#'
#' \strong{Usage:}
#'
#'  \code{mg$bin_coverages(bin_count=NA, region_filter=NA)}
#'
#' \strong{Description:}
#'
#' This method summarizes the coverage over regions of interests
#' into a specified number of bins. For each design group, it 
#' produces a matrix of binned coverages where each row represents a region,
#' and each column represents a bin. Those are returned in a named list where
#' each element contains the resulting matrix for a specific design group.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{bin_count}{The number of bins regions should be split into. The specified 
#'                     bin_count must always be equal or higher than the minimum size of
#'                     the specified regions. \code{NA} can be used to keep the previous
#'                     value. Default: \code{NA}.}
#'    \item{region_filter}{This parameter defines the subset of regions within the \code{regions}
#'                         parameter passed on initialization on which the metagene
#'                         should be generated. \code{region_filter} can be (1) a quosure, to be evaluated
#'                         in the context of the \code{region_metadata} data-frame, (2) a character
#'                         vector containing the names of the regions to be used or (3) a logical or numeric
#'                         vector to be used for subsetting. \code{NA} can be used to keep the previous
#'                         value. Default: \code{NA}}
#' }
#' @section split_coverages_by_regions:
#'
#' \strong{Usage:}
#'
#'  \code{mg$split_coverages_by_regions(split_by=NA)}
#'
#' \strong{Description:}
#'
#' This methods splits the matrices generated by mg$bin_coverages
#' into groups of regions where the values of the metadata columns
#' specified by \code{split_by} are homogeneous. It returns a list
#' where each element represents a design group: each of those
#' element is in turn a list representing groups of regions for which
#' all metadata values specified by "split_by" are equal. The leaf elements
#' of this list hierarchy are coverage matrices where each row represents a
#' region, and each column represents a bin.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{split_by}{A vector of column names from the region_metadata
#'                    parameter, as specified on metagene initialization. The 
#'                    selected columns must allow conversion into a factor.
#'                    By default, this is set to region_name, a metadata column
#'                    which is automatically generated by metagene. \code{NA} can 
#'                    be used to keep the previous value. Default: \code{NA}}
#' }
#' @section calculate_ci:
#'
#' \strong{Usage:}
#'
#'  \code{mg$calculate_ci(alpha = NA, sample_count = NA, resampling_strategy=NA)}
#'
#' \strong{Description:}
#'
#' This method calculates coverage means and confidence intervals for all 
#' design_group * region * bin combination. These are returned as a long-form
#' data-frame.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{alpha}{The alpha level of the confidence interval estimate.
#'                \code{NA} can be used to keep the previous value. 
#'                Default: \code{NA}}
#'    \item{sample_count}{The number of draws to perform in the bootstrap
#'                        calculations used to calculate the confidence inteval.  
#'                        \code{NA} can be used to keep the previous value. Default: \code{NA}}
#'    \item{resampling_strategy}{The resampling strategy to be used when performing the
#'                               bootstrap analysis, which can be either \code{'profile'}
#'                               or \code{'bin'}. In \code{'profile'} mode, whole profiles
#'                               across all bins are resampled. In \code{'bin'} mode,
#'                               each bin is resampled individually and independantly from
#'                               all others. \code{NA} can be used to keep the previous value.
#'                               Default: \code{NA}}
#' }
#' @section add_metadata:
#'
#' \strong{Usage:}
#'
#'  \code{mg$add_metadata(design_metadata=NA)}
#'
#' \strong{Description:}
#'
#' This method adds design group and region metadata to the data-frame
#' produced by \code{mg$calculate_ci} for easier plotting.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{design_metadata}{A data-frame containing metadata for the design groups.
#'                           It must contain as many rows as there are design groups,
#'                           and must contain at least one column named 'design'
#'                           which is used to match the rows to design groups.}
#' }
#' @section plot:
#'
#' \strong{Usage:}
#'
#'  \code{mg$plot(region_names = NULL, design_names = NULL,
#'                title = NA, x_label = NA, facet_by=NA, group_by=NA)}
#'
#' \strong{Description:}
#'
#' This method produces a ggplot object giving a graphical representation
#' of the metagene analysis.
#'
#' \strong{Parameters:}
#' \describe{
#'    \item{region_names}{The names of the regions to be plotted. If \code{NULL},
#'                        all the regions are plotted. Default: \code{NULL}.}
#'    \item{design_names}{The names of the design groups to be plotted. If 
#'                        \code{NULL}, all the design groups are
#'                        plotted. Default: \code{NULL}.}
#'    \item{title}{A title to add to the graph. \code{NA} can be used to keep 
#'                 the previous value. Default: \code{NA}}
#'    \item{x_label}{X-axis label for the metagene plot. \code{NA} can be 
#'                   used to keep the previous value. Default: \code{NA}.}
#'    \item{facet_by}{A formula to be used for facetting the metagene plot. This
#'                    formula can include any design metadata, or region_metadata 
#.                    columns that were part of the \code{split_by} argument.
#'                    \code{NA} can be used to keep the previous value.
#'                    Default: \code{NA}.}
#'    \item{group_by}{A string representing a single column from design_metadata or region_metadata
#'                    which will be used to group observations together into lines and which will
#'                    be used to generate the color scale.
#'                    \code{NA} can be used to keep the previous value.
#'                    Default: \code{NA}.}
#' }
#' @section Getter methods:
#' The following methods return various informations about the metagene object.
#'
#' \strong{mg$get_params()}
#' \describe{
#'    \item{}{Returns a list of all parameters used to perform this metagene analysis.}
#' }
#'
#' \strong{mg$get_design()}
#' \describe{
#'    \item{}{Returns the design used to perform this metagene analysis.}
#' }
#'
#' \strong{mg$get_regions()}
#' \describe{
#'    \item{}{Returns the regions used for this metagene analysis.}
#' }
#'
#' \strong{mg$get_data_frame(region_names = NULL, design_names = NULL)}
#' \describe{
#'    \item{}{Returns full data-frame of results.}
#'    \item{region_names}{The names of the regions to extract. If \code{NULL},
#'                        all the regions are returned. Default: \code{NULL}.}
#'    \item{design_names}{The names of the design groups to extract. If \code{NULL},
#'                        design groups are returned. Default: \code{NULL}.}
#' }
#'
#' \strong{mg$get_plot()}
#' \describe{
#'    \item{}{Returns the ggplot object generated by the \code{metagene2$plot} function.}
#' }
#'
#' \strong{mg$get_raw_coverages()}
#' \describe{
#'    \item{}{Returns raw coverages over the regions specified on initialization.}
#'}
#'
#' \strong{mg$get_normalized_coverages()}
#' \describe{
#'    \item{}{Returns normalized coverages over the regions specified on initialization.}
#'}
#'
#' @examples
#' mg <- metagene2$new(regions = get_demo_regions(), bam_files = get_demo_bam_files())
#' \dontrun{
#'    mg$plot()
#' }
#'
#' @import GenomicRanges
#' @import BiocParallel
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom R6 R6Class
#' @importFrom data.table data.table
#' @importFrom tools file_path_sans_ext
#' @importFrom rtracklayer import
#' @importFrom purrr map
#' @importFrom magrittr "%>%"
#' @export
#' @format A metagene experiment manager
metagene2 <- R6Class("metagene2",
    public = list(
    # Methods
        initialize = function(regions, bam_files, padding_size = 0,
                              cores = SerialParam(), verbose = FALSE,
                              force_seqlevels = FALSE, paired_end = FALSE,
                              assay = 'chipseq', strand_specific=FALSE,
                              paired_end_strand_mode=2,
                              region_mode="auto", region_metadata=NULL, 
                              extend_reads=0, invert_strand=FALSE, ...) {

            # Validate the format of bam_files, since it is used to preprocess certain
            # parameters before initialization.
            validate_bam_files_format(bam_files)
            
            if(region_mode=="auto") {
                region_mode = ifelse(assay=='rnaseq', "stitch", "separate")
            }
            
            # Define default design.
            default_design = private$get_complete_design(bam_files)
            design_metadata = data.frame(design=as.character(default_design[,1]), stringsAsFactors=FALSE)
            
            # Initialize parameter handler.
            private$ph <- parameter_manager$new(
                param_values=list(
                    design=default_design,
                    design_metadata=design_metadata,
                    bam_files=private$name_from_path(bam_files),
                    padding_size=padding_size,
                    verbose=verbose,
                    force_seqlevels=force_seqlevels,
                    paired_end=paired_end,
                    assay=tolower(assay),
                    strand_specific=strand_specific,
                    paired_end_strand_mode=paired_end_strand_mode,
                    normalization=NULL,
                    avoid_gaps=FALSE,
                    gaps_threshold=0,
                    bin_count=100,
                    alpha=0.05,
                    sample_count=1000,
                    region_mode=region_mode,
                    split_by="region_name",
                    region_filter=TRUE,
                    design_filter=TRUE,
                    resampling_strategy="bin",
                    facet_by=NULL,
                    group_by=NULL,
                    title=NULL,
                    x_label=NULL,
                    extend_reads=extend_reads,
                    invert_strand=invert_strand),
                param_validations=list(
                    design=private$validate_design,
                    bam_files=validate_bam_files,
                    padding_size=validate_padding_size,
                    verbose=validate_verbose,
                    force_seqlevels=validate_force_seqlevels,
                    assay=validate_assay,
                    normalization=validate_normalization,
                    avoid_gaps=validate_avoid_gaps,
                    gaps_threshold=validate_gaps_threshold,
                    bin_count=validate_bin_count,
                    alpha=validate_alpha,
                    sample_count=validate_sample_count,
                    region_mode=validate_region_mode,
                    extend_reads=validate_extend_reads,
                    split_by=validate_split_by),
                overall_validation=validate_combination)
                
            # Update any other parameter passed as a ... argument.
            # Since the parameter manager is locked, any non-existant
            # parameter will cause an error.
            private$ph$update_params(...)
                
            # Prepare objects for parralel processing.
            self$set_cores(cores)
            
            # Prepare bam files
            bm = private$start_bm("Prepare bam files")
            private$bam_handler <- Bam_Handler$new(private$ph$get("bam_files") , cores = cores,
                                        paired_end = paired_end,
                                        strand_specific=strand_specific,
                                        paired_end_strand_mode=paired_end_strand_mode,
                                        extend_reads=private$ph$get("extend_reads"),
                                        invert_strand=private$ph$get("invert_strand"))
            private$stop_bm(bm)

            # Prepare regions
            private$print_verbose("Prepare regions...")
            private$regions <- private$prepare_regions(regions, private$ph$get("region_mode"), region_metadata)

            # Parse bam files
            bm = private$start_bm("Producing coverage")
            private$coverages <- private$produce_coverages()
            private$stop_bm(bm)            
        },
        get_bam_count = function(filename) {
            # Parameters validation are done by Bam_Handler object
            private$bam_handler$get_aligned_count(filename)
        },
        get_params = function() {
            private$ph$get_all()
        },
        get_design = function() {
            private$ph$get("design")
        },
        get_design_group_names = function() {
            private$get_design_names_internal(private$ph$get('design'))
        },        
        get_regions = function() {
            return(private$regions)
        },
        get_regions_metadata = function() {
            return(private$region_metadata)
        },
        get_split_regions = function() {
            return(split_regions(private$regions,
                                 private$region_metadata, 
                                 private$ph$get("split_by")))
        },
        get_data_frame = function(region_names = NULL, design_names = NULL) {
            if(!is.null(private$ci_meta_df)) {
                results = private$ci_meta_df
            } else if(!is.null(private$ci_df)) {
                results = private$ci_df
            } else {
                return(NULL)
            }
            
            region_subset = TRUE
            if (!is.null(region_names)) {
                region_subset = results$region %in% region_names
            }
            
            design_subset = TRUE
            if (!is.null(design_names)) {
                design_subset = results$design %in% design_names
            }
            return(results[region_subset & design_subset,,drop=FALSE])
        },
        get_plot = function() {
            private$graph
        },
        get_raw_coverages = function() {
            if(!private$ph$get('strand_specific')) {
                return(private$get_raw_coverages_internal()[['*']])
            } else {
                return(private$get_raw_coverages_internal())
            }        
        },
        get_normalized_coverages = function() {
            if(!private$ph$get('strand_specific')) {
                return(private$get_normalized_coverages_internal()[['*']])
            } else {
                return(private$get_normalized_coverages_internal())
            }
        },
        set_cores = function(cores) {
            validate_cores(cores)
            if(is.null(private$parallel_job)) {
                private$parallel_job <- Parallel_Job$new(cores)
            } else {
                private$parallel_job$set_core_count(cores)
            }
        },
        group_coverages = function(design=NA, normalization=NA, design_filter=NA, simplify=TRUE) {
            # Clean up the design so it'll have the expected format.
            design = private$clean_design(design, private$ph$get("bam_files"))

            if(private$ph$have_params_changed(design)) {
                # design has changed.
                # Generate a new design metadata.
                private$ph$set("design_metadata", data.frame(design=design[,1]))
            }
            
            private$update_params_and_invalidate_caches(design, normalization, design_filter)
            
            if(is.null(private$grouped_coverages)) {
                bm <- private$start_bm("Grouping and normalizing coverages")
                
                # Identify design subset.
                design_col_to_keep = rep_len(private$ph$get("design_filter"), ncol(private$ph$get("design")) - 1)
                design_col_to_keep = 1 + which(design_col_to_keep)
                design_subset = private$ph$get("design")[, c(1, design_col_to_keep)]
                
                # Determine how data is to be merged
                if(is.null(private$ph$get("normalization"))) {
                    merge_operation = "+"
                } else if(private$ph$get("normalization")=="RPM") {
                    merge_operation = "mean"
                } else if(private$ph$get("normalization")=="NCIS") {
                    merge_operation = "NCIS"
                } else if(private$ph$get("normalization")=="log2_ratio") {
                    merge_operation = "log2_ratio"
                } else {
                    stop("Unsupported normalization value.")
                }
                
                private$grouped_coverages = lapply(private$get_coverages_internal(),
                                                   group_coverages_s,
                                                   design=design_subset,
                                                   bam_handler=private$bam_handler,
                                                   merge_operation=merge_operation)
                private$stop_bm(bm)                                                   
            }

            if(!private$ph$get("strand_specific") && simplify) {
                invisible(private$grouped_coverages[["*"]])
            } else {
                invisible(private$grouped_coverages)
            }
        },
        bin_coverages = function(bin_count=NA, region_filter=NA) {
            # Make sure the previous step has been performed.
            self$group_coverages()
            private$update_params_and_invalidate_caches(bin_count, region_filter)
            
            if(is.null(private$binned_coverages)) {
                bm = private$start_bm("Binning coverages")
                private$binned_coverages = bin_coverages_s(private$grouped_coverages,
                                                           regions=private$select_regions(private$ph$get("region_filter")),
                                                           bin_count=private$ph$get("bin_count"))
                private$stop_bm(bm)
            }
           
            invisible(private$binned_coverages)
        },
        split_coverages_by_regions = function(split_by=NA) {
            # Make sure the previous step has been performed.
            self$bin_coverages()
            private$update_params_and_invalidate_caches(split_by)
            
            if(is.null(private$split_coverages)) {
                bm = private$start_bm("Splitting coverages by region type")
                split_res = split_matrices(private$binned_coverages,
                                           private$select_region_metadata(private$ph$get("region_filter")),
                                           private$ph$get('split_by'))
                private$split_coverages = split_res$Matrices
                private$split_metadata_cache = split_res$Metadata
                private$stop_bm(bm)                                         
            }
            
            invisible(private$split_coverages)
        },        
        calculate_ci = function(alpha=NA, sample_count=NA, resampling_strategy=NA) {
            # Make sure the previous steps have been completed.
            self$split_coverages_by_regions()
            private$update_params_and_invalidate_caches(alpha, sample_count, resampling_strategy)

            if(is.null(private$ci_df)) {
                bm = private$start_bm("Producing data-frame")
                private$ci_df = calculate_matrices_ci(private$split_coverages,
                                                      private$ph$get('sample_count'), 
                                                      private$ph$get('alpha'),
                                                      private$ph$get('resampling_strategy'),
                                                      private$parallel_job)
                private$stop_bm(bm)
            }
            
            invisible(private$ci_df)
        },
        add_metadata = function(design_metadata=NA) {
            # Make sure the previous steps have been completed.
            self$calculate_ci()
            private$update_params_and_invalidate_caches(design_metadata)
                                                    
            if(is.null(private$ci_meta_df)) {
                filtered_design = private$ph$get("design_metadata")[private$ph$get("design_filter"),, drop=FALSE]
                private$ci_meta_df = add_metadata_to_ci(private$ci_df,
                                                        private$split_metadata_cache,
                                                        filtered_design)
            }
            
            invisible(private$ci_meta_df)
        },
        plot = function(region_names = NULL, design_names = NULL, title = NA,
                        x_label = NA, facet_by=NA, group_by=NA) {
            # 1. Get the correctly formatted table
            self$add_metadata()

            private$update_params_and_invalidate_caches(title, x_label, facet_by, group_by)
            
            plot_df <- self$get_data_frame(region_names, design_names)
            
            # 3. Produce the graph
            if (is.null(title)) {
                title <- paste(unique(plot_df[["group"]]), collapse=" vs ")
            }
            private$graph <- private$plot_graphic(df = plot_df, 
                                        title = private$ph$get("title"), 
                                        x_label = private$ph$get("x_label"),
                                        facet_by=private$ph$get("facet_by"),
                                        group_by=private$ph$get("group_by"))
            private$graph
        },
        produce_metagene = function(...) {
            private$update_params_and_invalidate_caches(...)
            self$plot()
        },
        plot_single_region = function(region, facet_by=NA, group_by=NA,
                                      no_binning=FALSE) {
            # Clone the mg object
            single_mg = self$clone(deep=TRUE)
            
            # Select the one region to plot, and make sure it ends up as a single GRange object.
            single_region = private$select_regions(region)
            stopifnot(length(single_region)==1)
            if(is(self$get_regions(), "GRangesList")) {
                single_region = single_region[[1]]
            }
            
            # If we are skipping binning, make the number of bins equal
            # the length in nucleotide.
            if(no_binning) {
                bin_count = sum(width(single_region))
            } else {
                bin_count = private$ph$get("bin_count")
            }
            
            # Re-bin with the new bin_count and the new filter.
            single_mg$bin_coverages(bin_count=bin_count, region_filter=region)
            
            # With a single region, splitting cannot be applied, so
            # we'll pass in a single row-name that we know to be valid.
            single_mg$split_coverages_by_regions(split_by="region_name")
            
            # Produce the new base single plot.
            out_plot = single_mg$plot(facet_by=facet_by, group_by=group_by)
            
            # Adjust the plot if binning was skipped.
            if(no_binning) {
                # New x label.
                out_plot <- out_plot + labs(x="Distance in nucleotides") 
                
                # If dealing with a stitched region, display "exon" boundaries as.
                # vertical lines.
                if(length(single_region) > 1) {
                    cumulative_width = 0
                    for(i in width(single_region)) {
                        cumulative_width = cumulative_width + i
                        out_plot <- out_plot + geom_vline(xintercept=cumulative_width)
                    }
                }
            }

            # Return the new plot.
            out_plot
        },
        replace_region_metadata = function(region_metadata) {
            # Validate that the old and new metadata have the same number of rows.
            if(nrow(region_metadata)!=nrow(private$region_metadata)) {
                stop("region_metadata must have one row per region.")
            }
            
            # Make sure the region_name column is still present.
            if(is.null(region_metadata$region_name)) {
                message("region_name is missing from the new metadata. Recreating it.")
                region_metadata$region_name = private$region_metadata$region_name
            }
            
            # Make sure that the split_by columns are all present and did not change.
            # If they did, invalidate everything after split_by
            for(split_column in private$ph$get("split_by")) {
                if(is.null(region_metadata[[split_column]]) ||
                   !all(region_metadata[[split_column]]==private$region_metadata[[split_column]])) {
                    # The split_by columns were removed or changed. Restore the original parameter value.
                    message("Replace region_metadata with metadata which would result in a different ",
                            "region split. All caches at the 'split_regions' step will be invalidated, ",
                            "and split_by will be reset to its default value.")
                    private$update_params_and_invalidate_caches(split_by="region_name")
                }
            }
            
            # Everything checks out, replace the metadata.
            private$region_metadata = region_metadata
        }
    ),
    private = list(
        # Region information. Both are kept separate from the parameter
        # handler since both can be very large, and making comparisons
        # to see if they've changed would be onerous.
        regions = GRangesList(),
        region_metadata = NULL,
        
        # Caches of intermediary step results.
        coverages = list(),
        grouped_coverages = NULL,
        binned_coverages = NULL,
        split_coverages = NULL,
        ci_df=NULL,
        ci_meta_df=NULL,
        graph = NULL,

        # Internal caches.
        split_metadata_cache = NULL,

        # Objects for handling bams, parameters and parallel jobs.
        bam_handler = "",
        parallel_job = NULL,
        ph=NULL,
        
        print_verbose = function(to_print) {
            if (private$ph$get("verbose")) {
                cat(paste0(to_print, "\n"))
            }
        },
        prepare_regions = function(regions, region_mode, region_metadata) {
            if (is(regions, "character")) {
                # Validation specific to regions as a vector
                if (!all(vapply(regions, file.exists, TRUE))) {
                    stop("regions is a list of files, but some of those files do not exist.")
                }
                regions = private$import_regions_from_disk(regions)
            }

            if (is(regions, "GRanges")) {
                if(region_mode=="stitch") {
                    if(length(unique(seqnames(regions))) > 1) {
                        stop(paste0("In stitch region_mode, such as in rnaseq assays, regions should be a ",
                                    "GRangesList of transcripts, or a GRanges ",
                                    " object representing a single transcript. ",
                                    "Here regions spans several seqnames, indicating ",
                                    "it might include many transcripts."))
                    }
                }
                regions <- GRangesList("regions" = regions)                
            } else if (is(regions, "list")) {
                regions <- GRangesList(regions)
            } else if (!is(regions, "GRangesList")) {
                stop(paste0("regions must be either a vector of BED ",
                            "filenames, a GRanges object or a GrangesList object"))            
            }
            
            # If regions do not have names, generate generic names for them.
            if (is.null(names(regions))) {
                names(regions) <- paste0("region_", seq_along(regions))
            }
            
            # In stitch mode, make sure disjoint regions part of the same
            # group do not overlap each other.
            if (region_mode=="stitch"){
                # If some "exons" overlap, then the total size will be smaller than the reduced size.
                total_size = sum(width(regions))
                reduced_size = sum(width(GenomicRanges::reduce(regions)))
                if(!all(total_size==reduced_size)) {
                    stop("In stitch region_mode, no overlap should exist between the individual ",
                         "GRanges making up the elements of the GRangesList")
                }
            }
            # TODO: Check if there is a id column in the mcols of every ranges.
            #    If not, add one by merging seqnames, start and end.

            # Apply padding and sortSeqLevels
            pad_regions = function(x, padding_size) {
                start(x) <- pmax(start(x) - padding_size, 1)
                end(x) <- end(x) + padding_size
                # Clean seqlevels
                x <- sortSeqlevels(x)
                x
            }
            regions = GRangesList(lapply(regions, pad_regions, padding_size = private$ph$get("padding_size")))
            
            # Add a region column to all GRangesList elements.
            regions_with_extra_col = list()
            for(region_name in names(regions)) {
                regions_with_extra_col[[region_name]] = regions[[region_name]]
                mcols(regions_with_extra_col[[region_name]])$region_name = region_name
            }
            regions = GRangesList(regions_with_extra_col)
            
            # In separate mode, simplify regions into a single GRanges object.
            if(region_mode=="separate") {
                regions = unlist(regions, use.names=FALSE)
            }

            # Build metadata from the mcols of the given regions.
            if(region_mode=="separate") {
                mcol_metadata = mcols(regions)
            } else {
                first_or_null = function(x) {
                    if(length(x)>0) {
                        return(mcols(x)[1,, drop=FALSE])
                    } else {
                        return(NULL)
                    }
                }
                mcol_metadata = do.call(rbind, lapply(regions, first_or_null))
                rownames(mcol_metadata) = names(regions)
            }

            # Merge the passed metadata object with the mcol metadata.
            if(is.null(region_metadata)) {
                private$region_metadata = mcol_metadata
            } else {
                stopifnot(nrow(region_metadata)==length(regions))
                non_duplicate_columns = setdiff(colnames(mcol_metadata), colnames(region_metadata))
                if(length(non_duplicate_columns) > 0) {
                    private$region_metadata = cbind(region_metadata, mcol_metadata[,non_duplicate_columns, drop=FALSE])
                }
            }
            
            if(!is.null(names(regions)) && is.null(rownames(private$region_metadata))) {
                rownames(private$region_metadata) = names(regions)
            }            
            
            return(regions)
        },
        produce_coverages = function() {
            if(private$ph$get("region_mode")=="stitch") {
                regions = BiocGenerics::unlist(private$regions)
            } else {
                regions = private$regions
            }
            
            regions <- GenomicRanges::reduce(regions)
            bam_files = private$ph$get("bam_files")

            res <- private$parallel_job$launch_job(
                        data = bam_files,
                        FUN = private$bam_handler$get_coverage,
                        regions = regions,
                        force_seqlevels= private$ph$get("force_seqlevels"),
                        simplify=FALSE)

            names(res) <- bam_files
            
            # Turn res inside out so that strand is at the top level,
            # and bam files on the second.
            res = list('+'=purrr::map(res, '+'),
                       '-'=purrr::map(res, '-'),
                       '*'=purrr::map(res, '*'))
            replace_nulls = function(x) { 
                if(all(purrr::map_lgl(x, is.null))) {
                    return(NULL)
                } else {
                    return(x)
                }
            }
            res = lapply(res, replace_nulls)

                       
            sortseq_or_null <- function(x) {
                if(is.null(x)) {
                    return(x)
                } else {
                    return(lapply(x, GenomeInfoDb::sortSeqlevels))
                }
            }
            
            lapply(res, sortseq_or_null)
        },
        plot_graphic = function(df, title, x_label, facet_by, group_by) {
            # Prepare x label
            assay = private$ph$get("assay")
            if (is.null(x_label)) {
                if (assay == "chipseq") {
                    x_label <- "Distance in bins"
                } else if (assay == "rnaseq") {
                    x_label = ifelse(is.null(private$ph$get("bin_count")),
                                     "Distance in nucleotides",
                                     "Distance in bins")
                }
            }

            # Prepare y label
            y_label <- "Mean coverage"
            if (is.null(private$ph$get("normalization"))) {
                y_label <- paste(y_label, "(raw)")
            } else if(private$ph$get("normalization") == "log2_ratio") {
                y_label <- "log2((Treatment RPM + 1) / (Control RPM + 1))"
            } else {
                y_label <- paste(y_label, "(RPM)")
            }
            
            # Produce plot
            p <- plot_metagene(df, facet_by=facet_by, group_by=group_by) +
                ylab(y_label) +
                xlab(x_label) +
                ggtitle(title)
            p
        },
        get_bam_names = function(filenames) {
            if (all(filenames %in% colnames(private$ph$get("design"))[-1])) {
                filenames
            } else {
                stopifnot(private$check_bam_files(filenames))
                vapply(filenames,
                    private$bam_handler$get_bam_name,
                    character(1))
            }
        },
        check_bam_files = function(bam_files) {
            all(vapply(bam_files,
            function(x) {
                !is.null((private$bam_handler$get_bam_name(x)))
            },
            logical(1)))
        },
        get_design_names_internal = function(design) {
            if(is.null(design)) {
                return(NULL)
            } else {
                return(colnames(design)[-1])
            }
        },
        get_design_number = function(design) {
            if(is.null(design)) {
                return(NULL)
            } else {
                return(ncol(design) - 1)
            }        
        },
        get_bam_in_design = function(design, design_name) {
            private$get_x_in_design(design, design_name, 1)
        },
        get_control_in_design = function(design, design_name) {
            private$get_x_in_design(design, design_name, 2)        
        },
        get_x_in_design = function(design, design_name, value) {
            if(is.null(design)) {
                return(NULL)
            } else {
                return(design$Samples[design[[design_name]] == value])
            }        
        },
        get_bam_by_design = function(design) {
            map(private$get_design_names_internal(design), ~private$get_bam_in_design(design, .x))
        },
        get_coverage_names = function(coverages) {
            stopifnot(length(setdiff(names(coverages), c("+", "-", "*")))==0)
            if(!is.null(coverages[['+']])) {
                return(names(coverages[['+']]))
            } else if(!is.null(coverages[['-']])) {
                return(names(coverages[['-']]))
            } else if(!is.null(coverages[['*']])) {
                return(names(coverages[['*']]))
            }
        },
        get_raw_coverages_internal = function() {
            private$coverages
        },
        get_normalized_coverages_internal = function() {
            # Define a function which will normalize coverage for a single
            # BAM file.
            normalize_coverage <- function(work_item) {
                weight <- 1 / (work_item$Count / 1000000)
                return(list(Strand=work_item$Strand, 
                            BamFile=work_item$BamFile, 
                            Coverage=work_item$Coverage * weight))
            }
            
            # Get the raw coverages.
            coverages <- private$get_raw_coverages_internal()
            
            # Serialize the workload
            work_items = list()
            for(strand in c("+", "-", "*")) {
                if(!is.null(coverages[[strand]])) {
                    for(bam_file in names(coverages[[strand]])) {
                        work_items[[length(work_items) + 1]] = list(Coverage=coverages[[strand]][[bam_file]],
                                                                  Strand=strand,
                                                                  BamFile=bam_file,
                                                                  Count=private$bam_handler$get_aligned_count(bam_file))
                    }
                }
            }
            
            serialized_coverages <- private$parallel_job$launch_job(data = work_items,
                                                                    FUN = normalize_coverage)
                                                                    
            # Deserialize the coverages.
            norm_coverages = list("+"=NULL, "-"=NULL, "*"=NULL)
            for(i in serialized_coverages) {
                if(is.null(norm_coverages[[i$Strand]])) {
                    norm_coverages[[i$Strand]] = list()
                }
                
                norm_coverages[[i$Strand]][[i$BamFile]] = i$Coverage
            }

            norm_coverages
        },
        get_coverages_internal = function() {
            if (!is.null(private$ph$get("normalization"))) {
                coverages <- private$get_normalized_coverages_internal()    
            } else {
                coverages <- private$get_raw_coverages_internal()
            }
            
            return(coverages)
        },
        validate_design = function(design) {
            validate_design_format(design)
            private$validate_design_values(design)
        },
        validate_design_values = function(design) {
            # At least one file must be used in the design
            if (sum(rowSums(design[ , -1, drop=FALSE]) > 0) == 0) {
                stop("At least one BAM file must be used in the design.")
            }
        
            # Check if used bam files exist.
            non_empty_rows = rowSums(design[, -1, drop=FALSE]) > 0
            if (!all(file.exists(as.character(design[non_empty_rows,1])))) {
                warning("At least one BAM file does not exist")
            }
        },        
        get_complete_design = function(bam_files) {
            bam_files = private$name_from_path(bam_files)
            
            # Concatenate the bam names and the identity matrix, then
            # rename all columns but the first.
            design <- cbind(data.frame(bam_files = bam_files, stringsAsFactors=FALSE), diag(length(bam_files)))
            colnames(design)[-1] = names(bam_files)
            
            design
        },
        name_from_path = function(file_paths) {
            alt_names = tools::file_path_sans_ext(basename(file_paths))
            if(is.null(names(file_paths))) {
                names(file_paths) = alt_names
            } else {
                names(file_paths) = ifelse(names(file_paths)=="", alt_names, names(file_paths))
            }
            return(file_paths)
        },
        import_regions_from_disk = function(file_names) {
            file_names = private$name_from_path(file_names)
            import_file <- function(region) {
                ext <- tolower(tools::file_ext(region))
                if (ext == "narrowpeak") {
                    extraCols <- c(signalValue = "numeric",
                                    pValue = "numeric", qValue = "numeric",
                                    peak = "integer")
                    rtracklayer::import(region, format = "BED",
                                        extraCols = extraCols)
                } else if (ext == "broadpeak") {
                    extraCols <- c(signalValue = "numeric",
                                    pValue = "numeric", qValue = "numeric")
                    rtracklayer::import(region, format = "BED",
                                        extraCols = extraCols)
                } else if (ext == "gtf" | ext == "gff") {
                    gxf_regions = rtracklayer::import(region)
                    if("gene_id" %in% colnames(mcols(gxf_regions))) {
                        return(split(gxf_regions, gxf_regions$gene_id))
                    } else {
                        return(gxf_regions)
                    }
                } else {
                    rtracklayer::import(region)
                }
            }
            regions <- private$parallel_job$launch_job(data = file_names,
                                                    FUN = import_file)
            names(regions) <- names(file_names)
            return(regions)
        },
        deep_clone = function(name, value) {
            # With x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
            # each field, with the name and value.
            if (name == "ph") {
                 # `a` is an environment, so use this quick way of copying
                 value$clone()
            } else {
                # For all other fields, just return the value
                value
            }
        },
        clean_design = function(design, bam_files) {
            # If no design is provided, use the default one.
            if(is.null(design)) {
                design = private$get_complete_design(bam_files)
            }
            
            # NA will be overwritten with the previous design later on.
            if(!is.data.frame(design) && is.na(design)) {
                return(NA)
            }
            
            # Make sure the design is in the correct format (data-frame
            # with at least 2 columns) before we try improving it.
            validate_design_format(design)
            
            # Standardize names used in first column to match those
            # used in bam_files.
            if(!all(as.character(design[,1]) %in% bam_files)) {
                stop("Design contains samples absent from the list of bam files provided on initialization.")
            }

            return(design)
        },
        start_bm = function(msg) {
            private$print_verbose(paste0(msg, "..."))
            #return(list(Message=msg, Time=Sys.time(), Memory=pryr::mem_used()))
            return(list(Message=msg, Time=Sys.time()))
        },
        stop_bm = function(bm_obj) {
            bm_after_time = Sys.time()
            #bm_after_mem = pryr::mem_used()            
            bm_time = difftime(bm_after_time, bm_obj$Time, unit="secs")
            #bm_mem = structure(bm_after_mem - bm_obj$Memory, class="bytes")
            
            private$print_verbose(paste0("BENCHMARK-TIME-", bm_obj$Message, ":", bm_time))
            #private$print_verbose(paste0("BENCHMARK-MEMORY-", bm_obj$Message, ":", bm_mem))        
        },
        update_params_and_invalidate_caches = function(...) {
            # This prologue makes it possible to infer parameter names from the
            # name of the variable it is passed in. This allows us to avoid
            # design=design, bin_count=bin_count repetitive code.
            #
            # It cannot be factorized into a function, since in any further call,
            # the argument list will deparse as "list(...)".
            param_names_alt = unlist(lapply( substitute(list(...)), deparse)[-1])
            arg_list = list(...)
            if(is.null(names(arg_list))) {
                names(arg_list) = param_names_alt
            } else {
                names(arg_list) = ifelse(names(arg_list)=="", param_names_alt, names(arg_list))
            }
        
            # Associate each parameter witht he step it is used in.
            param_step_map = c(design="group_coverages",
                               normalization="group_coverages", 
                               noise_removal="group_coverages", 
                               design_filter="group_coverages",
                               bin_count="bin_coverages", 
                               region_filter="bin_coverages",
                               split_by="split_coverages",
                               region_filter="split_coverages",
                               alpha="calculate_ci", 
                               sample_count="calculate_ci", 
                               resampling_strategy="calculate_ci",
                               design_metadata="add_metadata")
                               
            # Associate each step with the cache it generates,
            # in reverse order, so we can easily determine which
            # caches to invalidate when a particular step needs to be re-run.
            step_cache_map = c(add_metadata="ci_meta_df",
                               calculate_ci="ci_df",
                               split_coverages="split_coverages",
                               bin_coverages="binned_coverages",
                               group_coverages="grouped_coverages")
                               
            cache_invalidated=FALSE
            
            # Loop over all passed-in parameters.
            if(length(arg_list) > 0) {
                for(arg_index in seq_along(arg_list)) {
                    arg_name=names(arg_list)[arg_index]
                    
                    # Determine if the parameter has changed from its last value.
                    if(do.call(private$ph$have_params_changed, arg_list[arg_index])) {
                        private$print_verbose(paste0(arg_name, " has changed.\n"))
                        
                        # Determine which step the parameter belongs to.
                        invalidated_step = param_step_map[names(arg_list)[arg_index]]
                        if(!is.na(invalidated_step)) {
                            # Invalidate all caches for the step the parameter belonged to,
                            # as well as all caches for downsteam steps.
                            invalidated_caches = step_cache_map[seq_len(which(names(step_cache_map)==invalidated_step))]
                            private$print_verbose(paste0(paste0(invalidated_caches, collapse=", "), " will be invalidated.\n"))
                            for(cache in invalidated_caches) {
                                private[[cache]] = NULL
                                cache_invalidated = TRUE
                            }
                        }
                    }
                }
            
                do.call(private$ph$update_params, arg_list)
            }
            
            return(cache_invalidated)
        },
        select_region_indices = function(selector) {
            if("quosure" %in% class(selector)) {
                # Using dplyr for this because I'm not comfortable enough with
                # quosures.
                selected_indices = as.data.frame(private$region_metadata) %>%
                                       dplyr::mutate(METAGENE_IDX=seq_len(n())) %>%
                                       dplyr::filter(!! selector) %>%
                                       dplyr::pull(METAGENE_IDX)
            } else if(is(selector, "character")) {
                if(!is.null(names(self$get_regions()))) {
                    selected_indices = selector
                } else {
                    selected_indices = self$get_regions()$region_name %in% selector
                }
            } else if(is.numeric(selector) || is.logical(selector) || is.numeric(selector)) {
                selected_indices = selector
            }
            
            selected_indices
        },
        select_regions = function(selector) {
            self$get_regions()[private$select_region_indices(selector)]
        },
        select_region_metadata = function(selector) {
            private$region_metadata[private$select_region_indices(selector),, drop=FALSE]
        }
    )
)