#!/usr/bin/env bash

show_usage () {

    # Print error message if one provided
    if [ -n "$1" ]; then
        echo "$1"
    fi

    echo "Tool for extracting specific columns from a GWAS metadata file. Outputs to stdout"
    echo "usage: subset_gxg_cols <gwas_metadata_file> <cols>"
    echo "Arguments:"
    echo $'\t' "gwas_metadata_file:     Path to gwas metadata input file"
    echo $'\t' "col:                    One or more space-delimited column (1-based) numbers to select"

    # Exit with error status
    exit 1
}

# Check to see if someone wants to see the help menus
if [ "$1" == "-h" ]; then
    show_usage
elif [ "$1" == "--help" ]; then
    show_usage
# Check number of arguments
elif [ $# -lt 2 ]; then
    show_usage "Incorrect number of arguments!"
fi

# Get input file
INPUT=$1

# Shift to be able to pass columns to cut function
shift

# Cut out the columns specified on the command line and replace certain words
cat $INPUT | cut -d$'\t' -f "$*" \
    | sed 's/Allele1/A1/g' \
    | sed 's/Allele2/A2/g' \
	| sed 's/Effect/BETA/g' \
	| sed 's/P.value/P/g'