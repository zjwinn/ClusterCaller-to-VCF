### FOR DEBUG ###
keyfile <- "klustercaller_keyfile_example.txt"
kc_file <- "klustercaller_file_example.txt"
verbose <- TRUE
out_file <- "klustercaller_to_vcf_example_output"

# Get command-line arguments from Bash
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
# keyfile <- args[1]
# kc_file <- args[2]
# verbose <- ifelse(args[3]=="true", TRUE, FALSE)
# out_file <- args[4]

# Check for verbose
if(verbose==TRUE){
  
  # Get message
  print_string="### Initiating conversion of KlusterCaller data to VCF in R! ###"
  
  # Measure string
  n=nchar(print_string)
  
  # Print message
  print(paste(rep("#", n), collapse = ""))
  print(print_string)
  print(paste(rep("#", n), collapse = ""))
  
}


# Read in data
kasp_data <- read.table(kc_file, header=TRUE, sep="\t", na.strings=c("", "NA"), check.names=FALSE)
keyfile <- read.table(keyfile, header=TRUE, sep="\t", na.strings=c("", "NA"), check.names=FALSE)

# Check for verbose
if(verbose==TRUE){
  
  # Display head of KlusterCaller data
  print("### KlusterCaller head")
  print(kasp_data[1:5,1:3])
  
  # Display head of KlusterCaller keyfile
  print("### Keyfile head")
  print(keyfile[1:5,1:5])
  
}

# Pull list of markers in KlusterCaller data
markers_kasp_data <- colnames(kasp_data)[2:ncol(kasp_data)]

# Pull list of markers in keyfile
markers_keyfile <- keyfile[,"marker"]

# Pull all markers in the keyfile found in the kasp data
markers_keyfile <- markers_keyfile[markers_keyfile %in% markers_kasp_data]

# Check for duplicated markers in either file
if(any(duplicated(markers_kasp_data))){
  
  # If there is an error, print the following
  print("Error: There are duplicated markers in the ClusterCaller file. Check inputs and resubmit with unique marker names!")
  stop("There are duplicated markers in the ClusterCaller file. Check inputs and resubmit with unique marker names!")
  
}else if(any(duplicated(markers_keyfile))){
  
  # If there is an error, print the following
  print("Error: There are duplicated markers in the key file. Check inputs and resubmit with unique marker names!")
  stop("There are duplicated markers in the key file. Check inputs and resubmit with unique marker names!")
  
}

# make a pattern for illegal characters
illegal_pattern <- "[^XY:No Call]"

# Check if all markers are found in both files
if(length(markers_keyfile)==length(markers_kasp_data)){

  # Make an empty vector for the vcf output
  vcf<-c()
  
  # Run for loop
  for(i in unique(markers_kasp_data)){

    # Print
    if(verbose==TRUE){print(paste("### Formatting marker =", i))}
    
    # Get colnames
    temp1 <- c(colnames(kasp_data)[1], i)
    
    # Pull data
    temp1 <- as.data.frame(kasp_data[,colnames(kasp_data) %in% temp1])
    
    # Check illegal character is found
    if (any(grepl(illegal_pattern, temp1[,2]))) {

      # If illegal characters are found, throw an error
      print(paste("Error: Illegal character found in column ", i, " of ClusterCaller file!", sep = ""))
      stop("Illegal character found in column ", i, " of ClusterCaller file!")
      
    }
    
    # Make into rownames
    rownames(temp1) <- temp1[,1]
    
    # Remove column
    temp1 <- data.frame(temp1[,2],
                        row.names = rownames(temp1),
                        check.names = FALSE)
    
    # Pull key
    temp2 <- keyfile[keyfile[,"marker"]==i,]
    
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
    temp1[,1] <- suppressWarnings(gsub(temp2[1,"X"], 0, temp1[,1]))    
    temp1[,1] <- suppressWarnings(gsub(temp2[1,"Y"], 1, temp1[,1]))
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
                        POS = temp4,
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
  
  # Make header text
  header <- c("##fileformat=VCFv4.2",
              '##clustercaller_to_vcf=<ID=GenotypeTable,Version=1.0,Description="KASP assays converted to VCF format. Missing positions reported as sudo-positions starting at 1.">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
  
  
  # Turn missing positions into a new position
  vcf_mod <- c()
  
  # Chromosome names
  chrs <- unique(vcf[,"#CHROM"])
  chrs <- chrs[order(chrs)]
  
  # For loop
  for (i in chrs){
    # Pull markers
    temp1 <- vcf[vcf[,"#CHROM"]==i & is.na(vcf[,"POS"])==FALSE,]
    
    # Separate out markers with "." for those positions
    temp2 <- vcf[vcf[,"#CHROM"]==i & is.na(vcf[,"POS"])==TRUE,]
    
    # Check
    if (nrow(temp2)==0){
      # Place in the modified VCF
      vcf_mod<-rbind(vcf_mod, temp1)
      # Add row to header
      header <- c(header, paste("##contig=<ID=",
                                i,
                                ",length=",
                                max(temp1[,"POS"]),
                                ">",
                                sep=""))
      # Remove
      remove(temp1, temp2)
    }else if (nrow(temp1)>0 & nrow(temp2)>0){
      # Pull markers with positions and markers without position
      temp3 <- temp1[is.na(temp1[,"POS"])==FALSE,]
      # Assign number for position
      temp2[,"POS"] <- seq(1:nrow(temp2))
      # Rbind
      temp3 <- rbind(temp2, temp3)
      # Report
      vcf_mod<-rbind(vcf_mod, temp3)
      # Add row to header
      header <- c(header, paste("##contig=<ID=",
                                i,
                                ",length=",
                                max(temp3[,"POS"]),
                                ">",
                                sep="")) 
      # Remove
      remove(temp1, temp2, temp3)
    }else if (nrow(temp1)==0){
      # Assign number for position
      temp2[,"POS"] <- seq(1:nrow(temp2))
      # Place in the modified VCF
      vcf_mod<-rbind(vcf_mod, temp2)
      # Add row to header
      header <- c(header, paste("##contig=<ID=",
                                i,
                                ",length=",
                                max(temp2[,"POS"]),
                                ">",
                                sep=""))
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
  
  # Add metadata header
  header_alt <- matrix(nrow = length(header), ncol = ncol(vcf))
  header_alt[,1] <- as.character(header)
  
  # Turn NA into nothing
  header_alt[is.na(header_alt)] <- ""
  
  # Bind header and body together
  vcf <- rbind(header_alt, vcf)
  
  # Write vcf
  write.table(vcf, 
              paste(out_file, ".vcf", sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE,
              sep = "\t")
  
# Else if there is not a complete list of markers in both files  
}else{
  
  # Throw error
  if(verbose==TRUE){
    
    # Print error for read out
    print("ERROR: Not all marker names are found in both files. Check case, presence, and spelling of markers listed below!")
    print("################################################")
    print("### Printing marker names in files for debug ###")
    print("################################################")
    print("Markers in KlusterCaller file yet not in key file:")
    print("------------------------------")
    string <- ifelse(length(markers_kasp_data[!markers_kasp_data %in% markers_keyfile])==0, 
                     "No Missing Markers!",
                     markers_kasp_data[!markers_kasp_data %in% markers_keyfile])
    print(string)
    print("######################################")
    
  }
  
  # Exit with error
  stop("Not all marker names are found in both files. Check case, presence, and spelling of marker names!")
  
}

# Displaying warnings
if(verbose==TRUE){warnings()}