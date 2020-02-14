# Title     : TODO
# Objective : TODO
# Created by: awaldrop
# Created on: 2020-02-07



############## Create CL parser
library(optparse)
parser = OptionParser(usage = "\n%prog [options] --salmon_dir <salmon_dir> --output_basename <output_basename>",
                      description = "Program for merging salmon quant transcript-level output files for downstream RNAseq analysis",
                      prog="Rscript merge_salmon_quant.R")
parser = add_option(object=parser, opt_str=c("--salmon_dir"), default=NULL, type="character",
                    help="[REQUIRED] Path to directory containing salmon output files")
parser = add_option(object=parser, opt_str=c("--output_basename"), default=NULL, type="character",
                    help="[REQUIRED] Basename to prefix output files")
############## Parse command line
argv = parse_args(parser)

# parse cmdline params
salmon_dir = argv$salmon_dir
output_basename = argv$output_basename

# Options
options(stringsAsFactors = F)
setwd(salmon_dir)

# Install packages
library(tximport)
library(readr)
library(GenomicFeatures)

# Make tx-to-gene mapping
gencode.gtf <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"
txdb <- makeTxDbFromGFF(file = gencode.gtf, format = "gtf",
                        dataSource = "GencodeV28 Main Annotation File",
                        organism = "Homo sapiens")
txdb.keys <- keys(txdb, keytype = "TXNAME")
tx2gene.v28 <- select(txdb, txdb.keys, "GENEID", "TXNAME")


# Import transcript level data from Salmon
salmon.files <- list.files(".", pattern = "quant.sf", full.names = F)
names(salmon.files) <- gsub("_quant.sf", "", salmon.files)
txi.dte.tx <- tximport(files = salmon.files, type = "salmon", txIn = T, txOut = T,
    countsFromAbundance = "no", tx2gene = tx2gene.v28,
    ignoreAfterBar = T)

# Rename transcripts
tx.rownames <- gsub(x = rownames(txi.dte.tx$counts), "\\|.+", "")
rownames(txi.dte.tx$counts) <- tx.rownames
rownames(txi.dte.tx$abundance) <- tx.rownames
rownames(txi.dte.tx$length) <- tx.rownames

# Summarize to gene
txi.gene.v28 <- summarizeToGene(txi.dte.tx, tx2gene = tx2gene.v28, countsFromAbundance = "no")

# Save data
dte_file = paste0(output_basename, "_salmon_dte_tx_data_gencode28.rds")
gene_file = paste0(output_basename, "_salmon_gene_data_gencode28.rds")
saveRDS(txi.dte.tx, dte_file)
saveRDS(txi.gene.v28, gene_file)