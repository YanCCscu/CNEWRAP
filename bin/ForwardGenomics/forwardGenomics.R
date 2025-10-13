#!/usr/bin/env Rscript
# Xavier Prudent, 2015

#######################################################
## Available approaches:
## - GLS on global percentID (between the taxa and the global ancestor) and phenotypes
## - branch method with the local percentID computed for every branch
## - Check for a perfect match between global percentID and phenotype
#######################################################


# get the dir of the forwardGenomics.R script
args <- commandArgs(trailingOnly = F)
scriptDir <- dirname(sub("--file=","",args[grep("--file",args)]))
 
## Functions for the initialization of the analysis
source(paste(scriptDir, "/forwardGenomics_initialization.R", sep=""))
       
## Functions for the analysis
source(paste(scriptDir, "/forwardGenomics_fullAnalysis.R", sep=""))

## Input arguments
readArguments()

logfile<-paste(out_data,"log",sep=".")
if (file.exists(logfile)) file.remove(logfile)
logcat <- function(..., file = logfile, append = TRUE) {
    cat(..., file = file, append = append)
  }
## Open the input files
inputFiles()
cat(" complete input files \n")
## Prepare the output file
outputFiles()
cat(" build output files \n")

## Run forward genomics 
loopElements()
cat(" compelete analysis $\n")

