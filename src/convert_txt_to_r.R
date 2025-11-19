## loading libraries
required_libraries <- c("optparse",
                        "data.table")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}

## Options
options(stringsAsFactors = FALSE)

## --------------- ##
## Parse arguments ##
## --------------- ##

option_list <- list(
    optparse::make_option(c("-i", "--input"), type = "character", default = NULL,
                        help = "Path to the input file"),
    optparse::make_option(c("-s", "--samples"), type = "character", default = NULL,
                        help = "Path to the file with sample names"),
    optparse::make_option(c("-e", "--edges"), type = "character", default = NULL,
                        help = "Path to the file with edge names"),
    optparse::make_option(c("-o", "--output"), type = "character", default = NULL,
                        help = "Path to the output file")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# write function to convert txt to Rdata
convert_txt_to_rdata <- function(lioness_network, lionesssamples, edges, rdata_file) {
    data <- fread(lioness_network)
    lioness <- as.data.frame(data)

    # read sample names
    samples_lioness <- readLines(lioness_samples)
    if(length(samples_lioness) != ncol(lioness)){
        stop("Number of samples does not match number of columns in lioness network")
    }

    #read edge names
    edge_names <- readLines(edges)
    if(length(edge_names) != nrow(lioness)){
        stop("Number of edges does not match number of rows in lioness network")
    }

    # add column names
    colnames(lioness) <- samples_lioness

    # add row names
    rownames(lioness) <- edge_names
    save(lioness, file = rdata_file)
}

# Assign arguments to variables
input <- opt$input
lioness_samples <- opt$samples
output <- opt$output
edges <- opt$edges

# Call the function to convert txt to Rdata
convert_txt_to_rdata(input, lioness_samples, edges, output)