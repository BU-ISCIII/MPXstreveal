################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(stringr)
library(gtools)
library(plyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
  make_option(c("-i", "--long_table"    ), type="character", default=NULL             , metavar="path"   , help="Path to LCRsearch parsed long table input file."   ),
  make_option(c("-p", "--parse_lcr"     ), type="character", default='LCR2'           , metavar="string" , help="LCR name to parse when pattern is 2."              ),
  make_option(c("-s", "--lcr_structure" ), type="character", default='AT'             , metavar="string" , help="LCR 2 character pattern structure"                 ),
  make_option(c("-l", "--left_flanking" ), type="character", default='TACATGTGTTTTAGA', metavar="string" , help="Left flanking pattern for LCR whose pattern is 2." ),
  make_option(c("-r", "--right_flanking"), type="character", default='GGGCAAAACATATAA', metavar="string" , help="Right flanking pattern for LCR whose pattern is 2.")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$long_table)){
  print_help(opt_parser)
  stop("Please provide a LCRsearch parsed long table input file.", call.=FALSE)
}

################################################
################################################
## FUNCTIONS                                  ##
################################################
################################################

STR_DP_plot <- function(input_table, font_size){
  ggplot(input_table, aes(
    x = STR_mark ,
    y=Sample_name,
    fill = Supporting_reads
  )) +
    geom_tile()+
    scale_fill_gradient(low="cornflowerblue", high="black") +
    facet_grid(~STR_mark, scales = "free", switch = "x") +
    ggtitle("Depth") +
    labs(x = "LCR", y = "Sample", fill = "Supporting\nReads") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank()
    )
}

Allele_AF_HM_plot <- function(input_table, text_size){
  ggplot(input_table, aes(
    x = STR_structure,
    y=Sample_name,
    fill = AlleleFrequency
  )) +
    geom_tile()+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    scale_fill_gradient(low="cornflowerblue", high="black") +
    facet_grid(~STR_mark, scales = "free_x", space = "free_x", switch = "x") +
    ggtitle("Allele Frequency") +
    labs(x = "Allele", y = "Sample", fill = "Allele\nFrequency") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1, size = text_size),
    )
}

################################################
################################################
## READ IN LONG TABLE                         ##
################################################
################################################

STR_table <- read.table(file = opt$long_table, header = T, sep = "\t", stringsAsFactors = FALSE, colClasses = c("character", "character", "character","numeric", "numeric", "numeric", "character"))

################################################
################################################
## PARSE LCR OF PATTERN SIZE 2               ##
################################################
################################################

if (!is.null(opt$parse_lcr)) {
  lcr_long_structure <- paste(opt$lcr_structure, opt$lcr_structure, sep = "")
  for (line in 1:nrow(STR_table)) {
    STR_name <- as.character(STR_table[line,2])
    if (STR_name == opt$parse_lcr) {
      STR_structure <- as.character(STR_table[line,3])
      STR_Sequence <- as.character(STR_table[line,7])

      #Remove flanking regions
      STR_Sequence <-gsub(opt$left_flank,"",STR_Sequence)
      STR_Sequence <-gsub(opt$right_flank,"",STR_Sequence)

      #Remove everything between brackets and numbers to get strings outside the repeated pattern
      left_pattern <- gsub('[[:digit:]]+', '', gsub(" ", "", sub("\\[.*", "", STR_structure), fixed = TRUE))
      right_pattern <- gsub('[[:digit:]]+', '', gsub(" ", "", sub(".*\\]", "", STR_structure), fixed = TRUE))
      ##If something else in the Allele that is not the pattern, remove it before counting occurrences of short pattern.
      if (left_pattern != "" | right_pattern != "") {
        STR_Sequence <- gsub(left_pattern, "",STR_Sequence)
        STR_Sequence <- gsub(right_pattern, "",STR_Sequence)
        counts <- str_count(STR_Sequence, opt$lcr_structure)
      } else {
        counts <- str_count(STR_Sequence, opt$lcr_structure)
      }
      gsub_long <- paste("\\[", "\\]", sep = lcr_long_structure)
      gsub_short <- paste("\\[", "\\]", sep = opt$lcr_structure)
      new_STR_structure <- gsub('[[:digit:]]+', counts,STR_structure) #Replace old digit count with new digit
      new_STR_structure<- sub(gsub_long, gsub_short, new_STR_structure) #Replace ATAT with AT for example for LCR2
      STR_table$STR_structure[line]<- new_STR_structure
    }
  }
}

## Put LCR as factors to force order
STR_list_full <- mixedsort(unique(STR_table$STR_mark))
STR_table$STR_mark <- factor(STR_table$STR_mark,levels = STR_list_full)

## Filter data for Depth higher than 10X and Allele Frequency higher than 3%
STR_filtered <- subset(STR_table, Supporting_reads >= 10 & AlleleFrequency >= 0.03)


####### DEPTH PLOT
#Create table with the sum of depths for each LCR for all the alleles
STR_table_sum <- ddply(STR_table, .(STR_mark, Sample_name), summarize, Supporting_reads=sum(Supporting_reads))

#Create plot
svg("STR_nofilter_Depth.svg", width = 14, height = 8)
STR_DP_plot(STR_table_sum)
dev.off()

####### ALLELE FREQUENCY PLOT
# Select LCRs with higher variability and split 353_R and 349_R in different plots
STR_353_349 <- subset(STR_filtered, (Sample_name == "MPXV_353_R_SP_2022_NovaSeq" | Sample_name == "MPXV_349_R_SP_2022" | Sample_name == "MPXV_353_R_SP_2022_MiSeq") & (STR_mark == "LCR2" | STR_mark == "LCR5" | STR_mark == "LCR10" | STR_mark == "LCR11"))
STR_rest_samples <- subset(STR_filtered, Sample_name != "MPXV_353_R_SP_2022_MiSeq" & Sample_name != "MPXV_353_R_SP_2022_NovaSeq" & Sample_name != "MPXV_349_R_SP_2022" & (STR_mark == "LCR7" | STR_mark == "LCR10" | STR_mark == "LCR11" | STR_mark == "LCR12" | STR_mark == "LCR13" | STR_mark == "LCR14" | STR_mark == "LCR19" | STR_mark == "LCR20" | STR_mark == "LCR21"))

#Create plot
svg("353R_349R_filtered_AF_HM.svg", width = 12, height = 3)
Allele_AF_HM_plot(STR_353_349, 8)
dev.off()

#Create plot
svg("Rest_filtered_AF_HM.svg",  width = 15, height = 8)
Allele_AF_HM_plot(STR_rest_samples, 5)
dev.off()
