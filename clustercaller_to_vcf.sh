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
    echo -e "\tmay then be taken and used in downstream process."
    echo 
    echo "Options:"
    echo -e "\t-v, --verbose           Enable verbose mode"
    echo -e "\t-h, --help              Display this help and exit"
    echo
    echo "Arguments:"
    echo -e "\t-k, --key-file               A key file for interpretation of marker calls (tab delimited)"
    echo -e "\t-c, --clustercaller-file     ClusterCaller output file (tab delimited)"
    echo -e "\t-o, --out-file               Name of the vcf to write out (string with no spaces and no following .vcf)"
    echo
    echo "Examples:"
    echo -e "\tbash $0 -k keyfile_example.txt -c clustercaller_example.txt -o output_example -v"
    exit 1
}

# Default values
verbose=false

# Parse command line options
while getopts ":k:c:o:vh" opt; do
    case ${opt} in
        k | --key-file )
            key_file="$OPTARG"
            ;;
        c | --clustercaller-file )
            cc_file="$OPTARG"
            ;;
        o | --out-file )
            out_file="$OPTARG"
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
if [ -z "$key_file" ] || [ -z "$cc_file" ] || [ -z "$out_file" ]; then
    echo "Error: Required options are missing. Please provide key, clustercaller, and output file."
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

# Get realpath of files
key_file_realpath=$(realpath "$key_file")
cc_file_realpath=$(realpath "$cc_file")

# Check if the real paths exist
if [ ! -e "$key_file_realpath" ]; then
    echo "Error: The realpath of key file '$key_file_realpath' does not exist."
    exit 1
fi

if [ ! -e "$cc_file_realpath" ]; then
    echo "Error: The realpath of cluster caller file '$cc_file_realpath' does not exist."
    exit 1
fi

# Run R script
Rscript - "$key_file_realpath" "$cc_file_realpath" "$verbose" "$out_file" <<EOF

### R code goes here ###

# Get command-line arguments from Bash
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
key_file <- args[1]
cc_file <- args[2]
verbose <- ifelse(args[3]=="true", TRUE, FALSE)
out_file <- args[4]

# Check for verbose
if(verbose==TRUE){
    # Get message
    print_string="### Initiating conversion of ClusterCaller data to VCF in R! ###"

    # Measure string
    n=nchar(print_string)

    # Print message
    print(paste(rep("#", n), collapse = ""))
    print(print_string)
    print(paste(rep("#", n), collapse = ""))
}


# Read in data
kasp_data <- read.table(cc_file, header=TRUE, sep="\t", na.strings=c("", "NA"), check.names=FALSE)
key_file <- read.table(key_file, header=TRUE, sep="\t", na.strings=c("", "NA"), check.names=FALSE)

# Check for verbose
if(verbose==TRUE){
    # Display head
    print("")
    print("### ClusterCaller head")
    print(kasp_data[1:5,1:3])
    print("")
    print("### Keyfile head")
    print(key_file[1:5,1:5])
    print("")
}

# Pull list of markers
markers_kasp_data <- colnames(kasp_data)[2:ncol(kasp_data)]
markers_key_file <- key_file[,"marker"]
markers_key_file <- markers_key_file[markers_key_file %in% markers_kasp_data]

# Check if all markers are found in files
if(length(markers_key_file)==length(markers_kasp_data)){
  # Make an empty vector for the vcf output
  vcf<-c()
  
  # Run for loop
  for(i in markers_kasp_data){
    # Print
    if(verbose==TRUE){print(paste("### Formatting marker =", i))}
    
    # Get colnames
    temp1 <- c(colnames(kasp_data)[1], i)
    
    # Pull data
    temp1 <- as.data.frame(kasp_data[,colnames(kasp_data) %in% temp1])
    
    # Make into rownames
    rownames(temp1) <- temp1[,1]
    
    # Remove column
    temp1 <- data.frame(temp1[,2],
                        row.names = rownames(temp1),
                        check.names = FALSE)
    
    # Pull key
    temp2 <- key_file[key_file[,"marker"]==i,]
    
    # if
    if(as.character(temp2[1,"X"])=="N"){
      # Replace things in temp1
      temp1[,1] <- suppressWarnings(gsub("X", temp2[1,"X"], temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub("Y", temp2[1,"Y"], temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub("No Call", 
                                         paste(temp2[,"No Call"], 
                                               ":", 
                                               temp2[,"No Call"],
                                               sep = ""), 
                                         temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub(":", "/", temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub(temp2[1,"Y"], 0, temp1[,1]))
      #temp1[,1] <- suppressWarnings(gsub(temp2[1,"X"], 1, temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub("N", ".", temp1[,1]))
      temp1.1 <- as.vector(temp1[,1])
      names(temp1.1) <- rownames(temp1)
      temp1 <- temp1.1
      remove(temp1.1)
      
      # Pull info
      temp3 <- suppressWarnings(as.character(temp2[1,"chr"]))
      temp4 <- suppressWarnings(as.numeric(temp2[1,"position"]))
      temp5 <- suppressWarnings(as.character(temp2[1,"marker"]))
      temp6 <- suppressWarnings(as.character(temp2[1,"X"]))
      temp7 <- suppressWarnings(as.character(temp2[1,"Y"]))
      
      # Make vcf line
      temp3 <- data.frame("#CHROM" = temp3,
                          POS = ifelse(is.na(temp4), ".", temp4),
                          ID = temp5,
                          REF = ifelse(temp7=="N", ".", temp7),
                          ALT = ifelse(temp6=="N", ".", temp6),
                          QUAL = ".",
                          FILTER = "PASS",
                          INFO = ".",
                          FORMAT = "GT",
                          t(temp1),
                          check.names = FALSE)
      rownames(temp3) <- NULL
      
      # Rbind into vcf
      vcf <- rbind(vcf, temp3)
      
      # Remove
      remove(temp1, temp2, temp3, temp4, temp5, temp6, temp7)        
    }else{
      # Replace things in temp1
      temp1[,1] <- suppressWarnings(gsub("X", temp2[1,"X"], temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub("Y", temp2[1,"Y"], temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub("No Call", 
                                         paste(temp2[,"No Call"], 
                                               ":", 
                                               temp2[,"No Call"],
                                               sep = ""), 
                                         temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub(":", "/", temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub(temp2[1,"Y"], 0, temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub(temp2[1,"X"], 1, temp1[,1]))
      temp1[,1] <- suppressWarnings(gsub("N", ".", temp1[,1]))
      temp1.1 <- as.vector(temp1[,1])
      names(temp1.1) <- rownames(temp1)
      temp1 <- temp1.1
      remove(temp1.1)
      
      # Pull info
      temp3 <- suppressWarnings(as.character(temp2[1,"chr"]))
      temp4 <- suppressWarnings(as.numeric(temp2[1,"position"]))
      temp5 <- suppressWarnings(as.character(temp2[1,"marker"]))
      temp6 <- suppressWarnings(as.character(temp2[1,"X"]))
      temp7 <- suppressWarnings(as.character(temp2[1,"Y"]))
      
      # Make vcf line
      temp3 <- data.frame("#CHROM" = temp3,
                          POS = ifelse(is.na(temp4), ".", temp4),
                          ID = temp5,
                          REF = ifelse(temp6=="N", ".", temp6),
                          ALT = ifelse(temp7=="N", ".", temp7),
                          QUAL = ".",
                          FILTER = "PASS",
                          INFO = ".",
                          FORMAT = "GT",
                          t(temp1),
                          check.names = FALSE)
      rownames(temp3) <- NULL
      
      # Rbind into vcf
      vcf <- rbind(vcf, temp3)
      
      # Remove
      remove(temp1, temp2, temp3, temp4, temp5, temp6, temp7)           
    }
  }
  
  # Add metadata header
  header <- matrix(nrow = 3, ncol = ncol(vcf))
  
  # Turn NA into nothing
  header[is.na(header)] <- ""
  
  # Make header text
  header[,1] <- c("##fileformat=VCFv4.2",
                  '##clustercaller_to_vcf=<ID=GenotypeTable,Version=1.0,Description="KASP assays converted to VCF format. Missing positions reported as sudo-positions starting at 1.">',
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
  

  # Turn missing positions into a new position
  vcf_mod <- c()
  
  # Chromosome names
  chrs <- unique(vcf[,"#CHROM"])

  # For loop
  for (i in chrs){
    #print
    print(i)

    # Pull markers
    temp1 <- vcf[vcf[,"#CHROM"]==i & vcf[,"POS"]!=".",]
    
    # Separate out markers with "." for those positions
    temp2 <- vcf[vcf[,"#CHROM"]==i & vcf[,"POS"]==".",]
    
    # Check
    if (nrow(temp2)==0){
      # Place in the modified VCF
      vcf_mod<-rbind(vcf_mod, temp1)
      # Remove
      remove(temp1, temp2)
    }else if (nrow(temp1)>0 & nrow(temp2)>0){
      # print
      print(2) 
      # Pull markers with positions and markers without position
      temp3 <- vcf[vcf[,"POS"]!=".",]
      # Assign number for position
      temp2[,"POS"] <- seq(1:nrow(temp2))
      # Rbind
      temp3 <- rbind(temp2, temp3)
      # Report
      vcf_mod<-rbind(vcf_mod, temp3)
      # Remove
      remove(temp1, temp2, temp3)
    }else if (nrow(temp1)==0){
      # Assign number for position
      temp2[,"POS"] <- seq(1:nrow(temp2))
      # Place in the modified VCF
      vcf_mod<-rbind(vcf_mod, temp2)
      # Remove
      remove(temp1, temp2)
    }else{
      # Stop and exit with error
      stop("Something is wrong with markers on CHR = ",i)
      exit(status = 0)
    }
  }
  # Replace
  vcf <- vcf_mod

  # Make numeric
  vcf[,"POS"] <- as.numeric(vcf[,"POS"])
  
  # Order the vcf
  vcf <- vcf[order(vcf[,"#CHROM"], vcf[,"POS"]),]
  
  # Add header
  vcf <- rbind(colnames(vcf), vcf)
  
  # Make matrix
  vcf <- as.matrix(vcf)
  
  # Get rid of column names
  colnames(vcf) <- NULL
  
  # Bind header and body together
  vcf <- rbind(header, vcf)
  
  # Write vcf
  write.table(vcf, 
              paste(out_file, ".vcf", sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE,
              sep = "\t")
}else{
  # Throw error
  if(verbose==TRUE){
    print("################################################")
    print("### Printing marker names in files for debug ###")
    print("################################################")
    print("")
    print("######################################")
    print("Markers in ClusterCaller file:")
    print("------------------------------")
    print(markers_kasp_data[order(markers_kasp_data)])
    print("######################################")
    print("")
    print("######################################")
    print("Markers in Key file:")
    print("--------------------")
    print(key_file[order(key_file[,"marker"]),"marker"])
    print("######################################")
  }
  stop("Not all marker names are found in both files. Check case, presence, and spelling of names listed above!")
  quit(status = 0)
}

# Displaying warnings
if(verbose==TRUE){warnings()}
EOF

# Display step
if [ "$verbose" = true ]; then
echo
echo "###########################################"
echo "### Compressing and Indexing VCF output ###"
echo "###########################################"
fi

# Bgzip the file
bgzip -f $out_file.vcf
bcftools index $out_file.vcf.gz

# Display end
if [ "$verbose" = true ]; then
echo
echo "#############"
echo "### Done! ###"
echo "#############"
fi

# Exit without error
exit 1