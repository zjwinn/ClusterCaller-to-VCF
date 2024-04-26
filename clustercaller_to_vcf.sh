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
while getopts ":k:c:vh" opt; do
    case ${opt} in
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

# # Check if required options are provided
if [ -z "$key_file" ] || [ -z "$cc_file" ]; then
    echo "Error: Required options are missing. Please provide genome_file, key_file, and clustercaller_file."
    usage
fi

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo "Error: R is not installed! Please install R staitistical language before running this script. See link or more details: https://cran.r-project.org/ "
    exit 1
fi


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



#run in R
Rscript - "$key_file" "$cc_file" <<EOF

### R code goes here ###

# Get command-line arguments from Bash
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
key_file <- args[1]
cc_file <- args[2]

# Get message
print_string="### Initiating conversion of ClusterCaller data to VCF in R! ###"

# Measure string
n=nchar(print_string)

# Print message
print(paste(rep("#", n), collapse = ""))
print(print_string)
print(paste(rep("#", n), collapse = ""))

# Read in data
kasp_data<-read.table(cc_file)
key_file<-read.table(key_file)

# Display head
print("### ClusterCaller head")
print(head(kasp_data))
print("### Keyfile head")
print(head(key_file))

EOF