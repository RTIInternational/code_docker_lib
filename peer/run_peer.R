# Modified from https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R
# Original Author: Francois Aguet
# Modified by: Bryan Quach <bquach@rti.org>

library(peer, quietly = T)
library(argparser, quietly = T)
peer.version <- 1.3 #software version

WriteTable <- function(data, filename, index.name) {
    datafile <- file(filename, open = "wt")
    on.exit(close(datafile))
    header <- c(index.name, colnames(data))
    writeLines(paste0(header, collapse = "\t"), con = datafile, sep = "\n")
    write.table(data, datafile, sep = "\t", col.names = F, quote = F)
}

# Generate usage doc and retrieve command line args
p <- arg_parser("Run PEER using the R interface. Probabilistic Estimation of Expression Residuals (PEER) is a method designed to estimate surrogate variables/latent factors/hidden factors that contribute to gene expression variability, but it can be applied to other data types as well. For more information, please refer to https://doi.org/10.1038/nprot.2011.457.",
    name = "Probabilistic Estimation of Expression Residuals (PEER)")
p <- add_argument(p, arg = "omic_data_file", 
    help = "A tab delimited file containing N+1 rows and G+1 columns, where N is the number of samples, and G is the number of features (genes, methylation sites, chromatin accessibility windows, etc.). The first row and column must contain sample IDs and feature IDs respectively. Feature values should be normalized across samples and variance stabilized.")
p <- add_argument(p, arg = "output_prefix", 
    help = "File name prefix for output files. To specify an output directory as well, use --output_dir.")
p <- add_argument(p, arg = "num_factors", type = "numeric",
    help = "Number of hidden factors to estimate.")
p <- add_argument(p, arg = "--covariates", help = "A tab delimited file containing a matrix of size N × C, where C is the number of known covariates to be included in association test regression models of downstream analyses. Examples of common covariates include sex, age, batch variables, and quality metrics. Categorical variables (e.g., batch number) have to be encoded as indicator variables, with a different binary variable for each category. For the indicator variables, a value of 1 signifies membership in the category and a value of 0 indicates otherwise.")
p <- add_argument(p, arg = "--alphaprior_a", default = 0.001,
    help = "Shape parameter of the gamma distribution prior of the model noise distribution.")
p <- add_argument(p, arg = "--alphaprior_b", default = 0.01,
    help = "Scale parameter of the gamma distribution prior of the model noise distribution.")
p <- add_argument(p, arg = "--epsprior_a", help = "", default = 0.1,
    help = "Shape parameter of the gamma distribution prior of the model weight distribution.")
p <- add_argument(p, arg = "--epsprior_b", help = "", default = 10,
    help = "Scale parameter of the gamma distribution prior of the model weight distribution.")
p <- add_argument(p, arg = "--tol", type = "numeric", default = 0.001,
    help = "Threshold for the increase in model evidence when optimizing hidden factor values. Estimation completes for a hidden factor when the increase in model evidence exceeds this value.")
p <- add_argument(p, arg = "--var_tol", type = "numeric", default = 0.00001,
    help = "Threshold for the variance of model residuals when optimizing hidden factor values. Estimation completes for a hidden factor when the variance of residuals is smaller than this value.")
p <- add_argument(p, arg = "--max_iter", type = "numeric", default = 1000,
    help = "Max number of iterations for updating values of each hidden factor.")
p <- add_argument(p, arg = "--output_dir", short = "-o", default = ".",
    help = "Directory in which to save outputs.")
p <- add_argument(p, arg = "--version", short = "-v", flag = T, 
    help = "Print PEER version number.")
argv <- parse_args(p)

# Quick execution for printing version number
if(argv$version){
    cat(paste0("PEER v", peer.version))
    quit(save = "no")
}

# Check validity of argument inputs
if(is.na(argv$omic_data_file)){ 
    stop(paste0("Error: Please provide an 'omic_data_file'. Use --help for more details.")) 
}
if(! file.exists(argv$omic_data_file)){ 
    stop(paste0("Error: ", argv$omic_data_file, 
        " not found. Check your file path and name.")) 
}
if(is.na(argv$output_prefix)){ 
    stop(paste0("Error: Please provide an 'output_prefix'. Use --help for more details.")) 
}
if(is.na(argv$num_factors) | argv$num_factors <= 0 | 
   !is.finite(argv$num_factors) | argv$num_factors != as.integer(argv$num_factors)){ 
    stop(paste0("Error: Please provide a valid number for 'num_factors'. Use --help for more details."))
}
if(argv$max_iter <= 0 | !is.finite(argv$max_iter) | 
   argv$max_iter != as.integer(argv$max_iter)){ 
    stop(paste0("Error: Please provide a valid number for --max_iter. Use --help for more details."))
}
if(argv$tol <= 0 | !is.finite(argv$tol)){ 
    stop(paste0("Error: Please provide a valid value for --tol. Use --help for more details."))
}
if(argv$var_tol <= 0 | !is.finite(argv$var_tol)){ 
    stop(paste0("Error: Please provide a valid value for --var_tol. Use --help for more details."))
}
if(argv$alphaprior_a < 0 | !is.finite(argv$alphaprior_a)){ 
    stop(paste0("Error: Please provide a valid value for --alphaprior_a. Use --help for more details."))
}
if(argv$alphaprior_b < 0 | !is.finite(argv$alphaprior_b)){ 
    stop(paste0("Error: Please provide a valid value for --alphaprior_b. Use --help for more details."))
}
if(argv$epsprior_a < 0 | !is.finite(argv$epsprior_a)){ 
    stop(paste0("Error: Please provide a valid value for --epsprior_a. Use --help for more details."))
}
if(argv$epsprior_b < 0 | !is.finite(argv$epsprior_b)){ 
    stop(paste0("Error: Please provide a valid value for --epsprior_b. Use --help for more details."))
}

# Create output directory if needed
dir.create(argv$output_dir, showWarnings = F)

# Load omic data
cat(paste0("Loading data from ", argv$omic_data_file, " ..."))
omic.data <- read.table(argv$omic_data_file, sep = "\t", header = T, 
    check.names = F, comment.char = "", row.names = 1)
omic.data <- as.matrix(omic.data)
n.samples <- nrows(omic.data)
n.features <- ncols(omic.data)
cat("Done.\n")
cat(paste0("Loaded data matrix with ", n.samples, " rows and ", 
    n.features, " columns.\n"))

# Set method parameters
cat(paste0("Setting initialization parameters ..."))
model <- PEER()
invisible(PEER_setNk(model, argv$num_factors))
invisible(PEER_setPhenoMean(model, omic.data))
invisible(PEER_setPriorAlpha(model, argv$alphaprior_a, argv$alphaprior_b))
invisible(PEER_setPriorEps(model, argv$epsprior_a, argv$epsprior_b))
invisible(PEER_setTolerance(model, argv$tol))
invisible(PEER_setVarTolerance(model, argv$var_tol))
invisible(PEER_setNmax_iterations(model, argv$max_iter))
if (!is.null(argv$covariates) && !is.na(argv$covariates)) {
    covar.df <- read.table(argv$covariates, sep = "\t", header = T, row.names = 1, as.is = T)
    covar.df <- sapply(covar.df, as.numeric)
    cat(paste0("  * including ", dim(covar.df)[2], " covariates", "\n"))
    invisible(PEER_setCovariates(model, as.matrix(covar.df)))  # samples x covariates
}
cat("Done.\n")

# Begin inference routine
cat(paste0("Estimating ", argv$num_factors, " hidden confounders ...\n"))
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
WriteTable(t(X), file.path(argv$output_dir, paste0(argv$output_prefix, ".PEER_covariates.txt")), "ID")  # format(X, digits = 6)
WriteTable(A, file.path(argv$output_dir, paste0(argv$output_prefix, ".PEER_alpha.txt")), "ID")
WriteTable(R, file.path(argv$output_dir, paste0(argv$output_prefix, ".PEER_residuals.txt")), "ID")
cat("done.\n")
