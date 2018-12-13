#!/usr/bin/env bash

show_usage () {

    # Print error message if one provided
    if [ -n "$1" ]; then
        echo "$1"
    fi

    echo "Tool for extracting gtex variants from a single chromosome. Outputs to stdout"
    echo "usage: split_gtex_by_chr <gtex_input_file> <chr>"
    echo "Arguments:"
    echo $'\t' "gtex_input_file:   Path to gtex input file"
    echo $'\t' "chr:               Name of chr to extract"

    # Exit with error status
    exit 1
}

# Check to see if someone wants to see the help menus
if [ "$1" == "-h" ]; then
    show_usage
elif [ "$1" == "--help" ]; then
    show_usage
# Check number of arguments
elif [ $# -ne 2 ]; then
    show_usage "Incorrect number of arguments!"
fi

# Input file from cmd line
INPUT=$1

# Chromosome
CHROM=$2

# Send header to stdout
grep "chromosome" $INPUT

# Loop through and send lines matching chromosome to stdout
perl -lane 'if ($F[0] == '$CHROM') { print; }' $INPUT