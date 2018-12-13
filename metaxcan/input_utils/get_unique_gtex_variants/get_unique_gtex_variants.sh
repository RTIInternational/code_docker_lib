#!/usr/bin/env bash

show_usage () {

    # Print error message if one provided
    if [ -n "$1" ]; then
        echo "$1"
    fi

    echo "Tool for extracting unique variants froma gtex file. Outputs to stdout"
    echo "usage: get_unique_gtex_variants <gtex_input_file>"
    echo "Arguments:"
    echo $'\t' "gtex_input_file:   Path to gtex input file"

    # Exit with error status
    exit 1
}

# Check to see if someone wants to see the help menus
if [ "$1" == "-h" ]; then
    show_usage
elif [ "$1" == "--help" ]; then
    show_usage
# Check number of arguments
elif [ $# -ne 1 ]; then
    show_usage "Incorrect number of arguments!"
fi

# Get input file
INPUT=$1

# Echo header to stdout
echo "variant_id"

# Output unique entries to stdout
grep -v "chromosome" $INPUT | sort | uniq