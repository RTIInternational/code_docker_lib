# Modified from https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R
# Original Author: Francois Aguet
# Modified by: Bryan Quach <bquach@rti.org>

library(peer, quietly = T)
library(argparser, quietly = T)

WriteTable <- function(data, filename, index.name) {
    datafile <- file(filename, open = "wt")
    on.exit(close(datafile))
    header <- c(index.name, colnames(data))
    writeLines(paste0(header, collapse = "\t"), con = datafile, sep = "\n")
    write.table(data, datafile, sep = "\t", col.names = F, quote = F)
}

# Generate usage doc and retrieve command line args
p <- arg_parser("Run PEER using the R interface. Probabilistic Estimation of Expression Residuals (PEER) is a method designed to estimate surrogate variables/latent factors/hidden factors that contribute to gene expression variability, but it can be applied to other data types as well. This usage guide describes input arguments in relation to gene expression data, but the descriptions generalize to other genomic data types. For more information, please refer to https://doi.org/10.1038/nprot.2011.457")
p <- add_argument(p, "expression_file", help = "A tab delimited file containing N+1 rows and G+1 columns, where N is the number of samples, and G is the number of genes. The first row and column must contain sample IDs and gene IDs respectively. Gene expression values should be normalized across samples and variance stabilized.")
p <- add_argument(p, "prefix", help = "")
p <- add_argument(p, "n", help = "Number of hidden confounders to estimate")
p <- add_argument(p, "--covariates", help = "A tab delimited file containing a matrix of size N × C, where C is the number of known covariates to be included in association test regression models of downstream analyses. Examples of common covariates include sex, age, batch variables, and quality metrics. Categorical variables (e.g., batch number) have to be encoded as indicator variables, with a different binary variable for each category. For the indicator variables, a value of 1 signifies membership in the category and a value of 0 indicates otherwise.")
p <- add_argument(p, "--alphaprior_a", help = "", default = 0.001)
p <- add_argument(p, "--alphaprior_b", help = "", default = 0.01)
p <- add_argument(p, "--epsprior_a", help = "", default = 0.1)
p <- add_argument(p, "--epsprior_b", help = "", default = 10)
p <- add_argument(p, "--max_iter", help = "", default = 1000)
p <- add_argument(p, "--output_dir", short = "-o", help = "Output directory", default = ".")
argv <- parse_args(p)


# Load omic data
cat(paste0("Loading data from ", argv$expression_file, " ..."))
if(! file.exists(argv$expression_file)){ stop(paste0("Error: ",
    argv$expression_file, " not found. Check your file path and name.")) }
omic.data <- read.table(argv$expression_file, sep = "\t", header = T, 
    check.names = F, comment.char = "", row.names = 1)
omic.data <- as.matrix(omic.data)
cat("Done.\n")

# run PEER
cat(paste0("PEER: estimating hidden confounders (", argv$n, ")\n"))
model <- PEER()
invisible(PEER_setNk(model, argv$n))
invisible(PEER_setPhenoMean(model, omic.data))
invisible(PEER_setPriorAlpha(model, argv$alphaprior_a, argv$alphaprior_b))
invisible(PEER_setPriorEps(model, argv$epsprior_a, argv$epsprior_b))
invisible(PEER_setNmax_iterations(model, argv$max_iter))
if (!is.null(argv$covariates) && !is.na(argv$covariates)) {
    covar.df <- read.table(argv$covariates, sep = "\t", header = T, row.names = 1, as.is = T)
    covar.df <- sapply(covar.df, as.numeric)
    cat(paste0("  * including ", dim(covar.df)[2], " covariates", "\n"))
    invisible(PEER_setCovariates(model, as.matrix(covar.df)))  # samples x covariates
}
time <- system.time(PEER_update(model))

X <- PEER_getX(model)  # samples x PEER factors
A <- PEER_getAlpha(model)  # PEER factors x 1
R <- t(PEER_getResiduals(model))  # genes x samples

# add relevant row/column names
c <- paste0("InferredCov",1:ncol(X))
rownames(X) <- rownames(omic.data)
colnames(X) <- c
rownames(A) <- c
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha
rownames(R) <- colnames(omic.data)
colnames(R) <- rownames(omic.data)

# write results
cat("PEER: writing results ... ")
WriteTable(t(X), file.path(argv$output_dir, paste0(argv$prefix, ".PEER_covariates.txt")), "ID")  # format(X, digits = 6)
WriteTable(A, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_alpha.txt")), "ID")
WriteTable(R, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_residuals.txt")), "ID")
cat("done.\n")
