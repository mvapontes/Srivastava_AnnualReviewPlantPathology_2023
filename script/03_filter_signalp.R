# Filter out proteins without signalP

# Author: MVAP

# Packages ----
# Remember to insert %>% is ctrl + Shift +M

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

# . get fasta file path ---- 
prot_dire <- list.files(path = "speciesSelected/", pattern = "fasta|faa",  full.names = T)
prot_file <- str_split(prot_dire, "/", simplify = T)[,2]
prot_file <- str_replace(prot_file, ".fasta|.faa", "")
fileList <- tibble(original=prot_dire, species=prot_file)


# . get signalp file path ----
signalP <- grep(list.files(path = "signalP_output/", recursive = FALSE, full.names = T), 
                pattern = "fasta|faa", invert = TRUE, value = TRUE)

signalp_file <- str_split(signalP, "/", simplify = T)[,2]
species <- gsub('_summary.signalp5', "", signalp_file)
signalPfile <- tibble(signalP = signalP, species = species)


# . file List ----
fileList2 <- fileList %>% 
  left_join(signalPfile, by = "species") %>% 
  drop_na()


# Code ----

for (line in 1:nrow(fileList2)){
  
  # read signalP file
  sig_file <- read.table(fileList2[line,]$signalP, 
                         fill = T,
                         sep = "\t")
  
  # filter proteins with signalP
  sig_pos <- sig_file %>% 
    filter(str_detect(V2, "SP")) %>% 
    filter(V3 >= 0.7) %>% # prob >= 0.7
    mutate(V1 = gsub("\\|", ":", V1)) # replace "|" with ":" otherwise I get an error downstream
  
  # read fasta file
  pep_file <- Biostrings::readAAStringSet(fileList2[line,]$original)

  proteins <-data.table(as.data.frame(pep_file), keep.rownames = TRUE) 
  
  proteins2 <- proteins %>% 
    mutate(V1 = gsub("\\|", ":", rn)) %>% # replace "|" with ":" otherwise I get an error downstream 
    mutate(x = gsub("\\*", "", x)) # replace "*" at the end of protein with nothing otherwise I get an error downstream
  
  # select proteins with signalP
  sel_fasta <- subset(proteins2, grepl(paste0(sig_pos$V1, collapse = "|"), V1,ignore.case = TRUE))
  
  set12<-apply(sel_fasta[ ,2], 1, paste, collapse = "")
  
  # get protein fomrat
  set12 <- AAStringSet(set12)
  
  # add the names of the proteins
  names(set12)<-sel_fasta$rn
  
  # compose output file
  tmpfile <- paste("signalP_output/", fileList2[line,]$species, "_signalp.faa", sep="")
  print (tmpfile)
  
  # write fasta file with the longest proteins per proteome
  writeXStringSet(set12, tmpfile)
  
}
