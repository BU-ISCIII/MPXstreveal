################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(plyr)
library(lattice)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
  make_option(c("-i", "--input_table"   ), type="character", default=NULL             , metavar="path"   , help="Path to input table with supporting reads and lineages" ),
  make_option(c("-c", "--column_list"   ), type="character", default=NULL             , metavar="string" , help="List of comma separated column names to be plot."       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$input_table)){
  print_help(opt_parser)
  stop("Please provide a LCRsearch parsed long table input file.", call.=FALSE)
}
if (is.null(opt$input_table)) {
  print_help(opt_parser)
  stop("Please provide a list of columns to use for plot", call.=FALSE)
}

print(opt$column_list)

################################################
################################################
## FUNCTIONS                                  ##
################################################
################################################

create_svg_plot <- function(file_name, columns_name, xlabel, plot_title){
  svg(file_name)
  print(
    ggplot(data=STR_sra, aes(x=get(columns_name), fill = Lineage)) +
    geom_bar() +
    labs(x = xlabel, y = "Number of samples", fill = "Lineage") +
    ggtitle(plot_title) +
    theme_bw()+
    scale_x_continuous(breaks = round(seq(0, max(get(columns_name,STR_sra), na.rm = TRUE), by = 5),1)) +
    expand_limits(x = 0, y = 0)
    )
  dev.off()
}

################################################
################################################
## LOAD DATA                                  ##
################################################
################################################

columns_list <-  unlist(strsplit(opt$column_list, ","))

STR_sra <- read.table(file = opt$input_table, header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = '-')

################################################
################################################
## CREATE PLOTS                               ##
################################################
################################################

for (n in 1:length(columns_list)) {
  columns_name <- columns_list[n]
  STR_id <- unlist(strsplit(x = columns_name, split = '\\.'))[1]
  file_name <- paste(STR_id, "_SRA.svg", sep = "")
  xlabel <- paste("Number of ", paste(STR_id, " repeats", sep = ""), sep = "")
  plot_title <- paste("Distribution of ", paste(STR_id, " repeats in SRA publicly available Nanopre Sequencing data", sep = ""), sep = "")
  create_svg_plot(file_name, columns_name, xlabel, plot_title)
}
