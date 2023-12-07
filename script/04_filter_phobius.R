# filter proteins with signalp and  without transmembrane domain identified with phobius

# Author: MVAP

# Packages ----

library('tidyverse')
library('here')
library('janitor')
library('skimr')
library('dplyr')
library("Biostrings")
library("GenomicRanges")
library("GenomicFeatures")
library("Rsamtools")
library("rtracklayer")
library("data.table")


# Read files ----

# get phobius files path 
pho_dire <- list.files(path = "phobius/", recursive = TRUE, full.names = T)
pho_file <- str_split(pho_dire, "/", simplify = T)[,2]
pho_file <- gsub("_phobius.out", "", pho_file)
fileList <- tibble(phobius=pho_dire, species=pho_file)

# get fasta file for signalp selected proteins
signalP <- list.files(path = "signalP_output/", pattern = "_signalp.faa", recursive = FALSE, full.names = T)
signalp_file <- str_split(signalP, "/", simplify = T)[,2]
signalPfile <- tibble(signalP = signalP, species = signalp_file)

# File List
fileList2 <- fileList %>% 
  left_join(signalPfile, by = "species") %>% 
  drop_na()


# Code ----

for (line in 1:nrow(fileList2)){

  # read phobius file
  phobius_file <- read.table(fileList2[line,]$phobius, fill = T,
                             sep = " ", header = T) 
                         
  # Filter proteins with signalP without transmembrane domain
  pho_pos <- phobius_file %>% 
    filter(SP == "Y")%>%
    filter(!grepl("i", PREDICTION))
    
  # read fasta file
  pep_file <- Biostrings::readAAStringSet(fileList2[line,]$signalP)
  
  proteins <-data.table(as.data.frame(pep_file), keep.rownames = TRUE) 
  
  proteins2 <- proteins %>% 
    mutate(V1 = gsub("\\|", ":", rn)) # replace "|" with ":" otherwise I get an error downstream
  
  sel_fasta <- subset(proteins2,
                      grepl(paste0(rownames(pho_pos), collapse = "|"), V1, 
                            ignore.case = TRUE))
  
  set12<-apply(sel_fasta[ ,2], 1, paste, collapse="")
  
  # get protein format
  set12 <- AAStringSet(set12)
  
  # add the names of the proteins
  names(set12)<-sel_fasta$rn
  
  # compose output file
  tmpfile <- paste("tmhmm_phobius/", fileList2[line,]$species, "_phobius.faa", sep="")
  print (tmpfile)
  
  # write fasta file with the longest proteins per proteome
  writeXStringSet(set12, tmpfile)
  
}
