#!/share/apps/R/bin/Rscript
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset input file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}


df <- read.table(opt$file, header=TRUE)
filtered_df <- df[which(rowMeans(df[, -1]) > 200), ] # remove probes with mean<200

write.table(filtered_df, file=opt$out,  quote=FALSE, sep="\t", row.names=FALSE)
