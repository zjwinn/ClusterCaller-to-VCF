#!/bin/bash

# Function to display script usage
usage() {
    echo
    echo "###############################################"
    echo "#                                             #"
    echo "#               Help Information              #"
    echo "#                                             #"
    echo "###############################################"
    echo
    echo "Usage:"
    echo -e "\t$0 [OPTIONS] ARGUMENT"
    echo
    echo "Description:"
    echo -e "\tThis script will take KASP assay calls made through ClusterCaller and"
    echo -e "\tformat them into a variant calling format (VCF) file. The output of this script"
    echo -e "\tmay then be taken and used in downstream process"
    echo 
    echo "Options:"
    echo -e "\t-v, --verbose           Enable verbose mode"
    echo -e "\t-h, --help              Display this help and exit"
    echo
    echo "Arguments:"
    echo -e "\t-f, --genome-file            Input reference genome file"
    echo -e "\t-k, --key-file               a key file for interpretation of marker calls (tab delimited)"
    echo -e "\t-c, --clustercaller-file     ClusterCaller output file (tab delimited)"
    echo
    echo "Examples:"
    echo -e "\tbash $0 -f 'example.fa' -p 200 -l 10 -c 'Chr1A' -a 'A' -r 'T'"
    exit 1
}

# Default values
verbose=false

# Parse command line options
while getopts ":f:k:c:vh" opt; do
    case ${opt} in
        f | --genome-file )
            genome_file="$OPTARG"
            ;;
        k | --key-file )
            key_file="$OPTARG"
            ;;
        c | --clustercaller-file )
            cc_file="$OPTARG"
            ;;
        v | --verbose )
            verbose=true
            ;;
        h | --help )
            usage
            ;;
        \? )
            echo "Error: Invalid option -$OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Error: Option -$OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done
shift $((OPTIND -1))

if [ "$verbose" = true ]; then
    # Print header
    echo
    echo "###############################################"
    echo "#                                             #"
    echo "#          ClusterCaller to VCF v1.0          #"
    echo "#                                             #"
    echo "###############################################"
    echo
    echo "Written by: Zachary J. Winn PhD"
    echo "Contact information:"
    echo -e "\tGovernment Email: zachary.winn@usda.gov"
    echo -e "\tPersonal Email: zwinn@outlook.com"
    echo
    echo "###############################################"
    echo "# WARNING: This program is not under warranty #"
    echo "#          Use at your own discretion!        #"
    echo "###############################################"
    echo
fi

# Function to replace genotype calls with nucleotides
replace_genotypes() {
    local marker=$1
    local genotype=$2
    case $genotype in
        "X:X") echo "${genotype_data[$marker,1]}" ;;
        "Y:Y") echo "${genotype_data[$marker,2]}" ;;
        "X:Y") echo "${genotype_data[$marker,3]}" ;;
        *) echo "$genotype" ;;
    esac
}

# Read the first dataframe into an array
mapfile -t dataframe1 < dataframe1.csv

# Read the second dataframe into an associative array
declare -A genotype_data
while IFS=$'\t' read -r marker x y no_call; do
    genotype_data["$marker,1"]=$x
    genotype_data["$marker,2"]=$y
    genotype_data["$marker,3"]=$no_call
done < dataframe2.csv

# Process and replace genotype calls in the first dataframe
for ((i = 0; i < ${#dataframe1[@]}; i++)); do
    line="${dataframe1[i]}"
    IFS=$'\t' read -r -a fields <<< "$line"
    for ((j = 1; j < ${#fields[@]}; j++)); do
        fields[$j]=$(replace_genotypes "${fields[0]}" "${fields[j]}")
    done
    echo -e "${fields[*]}"
done
with